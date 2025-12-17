#include "DelaunayTriangulation.hpp"

#include <numeric>
#include <assert.h>

namespace Arns
{

namespace Math
{

//[Predicates]

bool isPointInCircumcircle(const Vector2D& p, const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    float ax = a.x - p.x;
    float ay = a.y - p.y;
    float bx = b.x - p.x;
    float by = b.y - p.y;
    float cx = c.x - p.x;
    float cy = c.y - p.y;

    float det = (ax * ax + ay * ay) * (bx * cy - cx * by)
        - (bx * bx + by * by) * (ax * cy - cx * ay)
        + (cx * cx + cy * cy) * (ax * by - bx * ay);

    return det > 0; // det < 0 for CW
}

bool circumcircleSquared(const Vector2D& a, const Vector2D& b, const Vector2D& c, Vector2D& center, float& radius2)
{
    const Vector2D d = b - a;
    const Vector2D e = c - a;

    const float bl = d.lengthSquared();
    const float cl = e.lengthSquared();
    const float det = d.cross(e);

    if (approximatelyZero(det))
    {
        center = {0.0f, 0.0f};
        radius2 = -1.0f; // invalid
        return false;
    }

    const Vector2D offset(
        (e.y * bl - d.y * cl) * 0.5f / det,
        (d.x * cl - e.x * bl) * 0.5f / det
    );

    center = a + offset;
    radius2 = offset.lengthSquared();
    return true;
}

float circumradiusSquared(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    const Vector2D d = b - a;
    const Vector2D e = c - a;

    const float bl = d.lengthSquared();
    const float cl = e.lengthSquared();
    const float det = d.cross(e);

    if (approximatelyZero(det))
    {
        return (std::numeric_limits<float>::max)();
    }

    const Vector2D radius(
        (e.y * bl - d.y * cl) * 0.5 / det,
        (d.x * cl - e.x * bl) * 0.5 / det);

    return radius.lengthSquared();
}

//Todo: Move this somewhere else
template <class T>
T orient2D(T aX, T aY, T bX, T bY, T cX, T cY)
{
    return (bX - aX) * (cY - aY) - (bY - aY) * (cX - aX);
}

float orient2D(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool isCCW(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return orient2D(a, b, c) > FloatAbsEpsilon;
}

bool isCCW(float aX, float aY, float bX, float bY, float cX, float cY)
{
    return orient2D(aX, aY, bX, bY, cX, cY) > FloatAbsEpsilon;
}

bool isCW(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return orient2D(a, b, c) < -FloatAbsEpsilon;
}

bool isCW(float aX, float aY, float bX, float bY, float cX, float cY)
{
    return orient2D(aX, aY, bX, bY, cX, cY) < -FloatAbsEpsilon;
}

bool isColinear(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return approximatelyZero(orient2D(a,b,c));
}


//[]

bool isDelaunay(const TriangleMesh2D& mesh)
{
    const auto& verts = mesh.m_vertices;
    const auto& tris = mesh.m_triangles;

    for (const auto& tri : tris)
    {
        const Vector2D& A = verts[tri.v0];
        const Vector2D& B = verts[tri.v1];
        const Vector2D& C = verts[tri.v2];

        Vector2D center;
        float radius2;
        circumcircleSquared(A, B, C, center, radius2);

        // Skip degenerate triangles
        if (radius2 < 0)
            continue;

        // Check all other vertices
        for (size_t i = 0; i < verts.size(); ++i)
        {
            if (i == tri.v0 || i == tri.v1 || i == tri.v2)
                continue;

            float d2 = center.distanceSquared(verts[i]);

            if (approximatelyLess(d2, radius2))
                return false;
        }
    }
    return true;
}


/*
  Original Algorithm: Delaunator (JavaScript) by Mapbox, licensed under the ISC License
  Adaptation: Modern C++ translation, based on the above work
  License: ISC (see THIRD_PARTY_LICENSES.txt)
*/

class Delaunator
{
public:
    /* Special value indicating that no valid index exists. */
    static constexpr std::size_t INVALID_INDEX = std::numeric_limits<size_t>::max();

    // --- Core triangulation structures ---

    /**
     * @brief Indices into the input vertex array, grouped by triples (one triple per triangle).
     * triangles[3 * t + 0], triangles[3 * t + 1], triangles[3 * t + 2] are the vertices of triangle t.
     */
    std::vector<std::size_t> triangles;

    /**
     * @brief Halfedge connectivity: halfedges[e] = opposite halfedge index, or INVALID_INDEX if boundary edge.
     */
    std::vector<std::size_t> halfedges;

    // --- Convex hull structures ---

    /* Previous vertex index in the convex hull linked list. */
    std::vector<std::size_t> hull_prev;
    /* Next vertex index in the convex hull linked list. */
    std::vector<std::size_t> hull_next;

    /* One adjacent triangle index for each convex hull vertex. */
    std::vector<std::size_t> hull_tri;

    /* The starting vertex index for iterating the convex hull linked list. */
    std::size_t hull_start;

private:
    std::tuple<std::size_t, std::size_t, std::size_t>
        find_seed_triangle(const std::vector<Vector2D>& points);

    void init_hull(const std::vector<Vector2D>& points,
        std::size_t i0, std::size_t i1, std::size_t i2);

    void insert_point_into_hull(const std::vector<Vector2D>& points,
        size_t i, size_t start);

    void insert_points(const std::vector<Vector2D>& points,
        const std::vector<std::size_t>& sorted_ids,
        std::size_t i0, std::size_t i1, std::size_t i2);

public:
    /**
    @brief Builds a Delaunay triangulation for the given point set.
    @param inputPoints Array of 2D points to triangulate.
    */
    explicit Delaunator(const std::vector<Vector2D>& inputPoints);

    //[Retrieval Method's]

    const std::vector<size_t>& getTriangleIndices()
    {
        return triangles;
    }

    std::vector<TriangleIndices> getTriangleTriplets()
    {
        std::vector<TriangleIndices> triangleTriplets;
        triangleTriplets.reserve(triangles.size() / 3);
        for (size_t i = 0; i < triangles.size(); i += 3)
        {
            triangleTriplets.emplace_back(triangles[i + 0], triangles[i + 1], triangles[i + 2]);
        }
        return triangleTriplets;
    }

    TriangleMesh2D getTriangleMesh(const std::vector<Vector2D>& inputPoints)
    {
        return TriangleMesh2D(inputPoints, getTriangleTriplets());
    }

private:
    Vector2D m_center;
    std::size_t m_hash_size;
    std::vector<std::size_t> m_hash;
    std::vector<std::size_t> m_edge_stack;

    [[nodiscard]] static float pseudo_angle(float dx, float dy) noexcept;
    [[nodiscard]] std::size_t hash_key(const Vector2D& point) const noexcept;

    [[nodiscard]] std::size_t find_hull_start(const Vector2D& point) const noexcept;

    std::size_t legalize(const std::vector<Vector2D>& inputPoints, std::size_t a);
    std::size_t add_triangle(std::size_t i0, std::size_t i1, std::size_t i2,
        std::size_t a, std::size_t b, std::size_t c);
    void link(std::size_t a, std::size_t b);
};

//DBL
// --- Monotonic mapping of vector direction to [0, 1) without atan2
float Delaunator::pseudo_angle(float dx, float dy) noexcept
{
    const float adx = std::abs(dx);
    const float ady = std::abs(dy);

    if (adx + ady == 0.0) return 0.0;

    const float p = dx / (adx + ady);
    return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0;
}

//DBL
// --- Converts point direction to hash bucket index
std::size_t Delaunator::hash_key(const Vector2D& point) const noexcept
{
    const float dx = point.x - m_center.x;
    const float dy = point.y - m_center.y;

    const float angle_fraction = pseudo_angle(dx, dy);
    const std::size_t bucket = static_cast<std::size_t>(
        std::floor(angle_fraction * static_cast<float>(m_hash_size)));

    return bucket % m_hash_size;
}

std::size_t Delaunator::find_hull_start(const Vector2D& point) const noexcept
{
    size_t key = hash_key(point);
    for (size_t j = 0; j < m_hash_size; j++)
    {
        size_t start = m_hash[(key + j) % m_hash_size];
        if (start != INVALID_INDEX && start != hull_next[start])
        {
            return hull_prev[start]; // previous edge
        }    
    }
    return INVALID_INDEX; // fallback
}

std::tuple<std::size_t, std::size_t, std::size_t>
Delaunator::find_seed_triangle(const std::vector<Vector2D>& points)
{
    const std::size_t n = points.size();

    // --- Build BBox and calculate Center
    BBox2D bbox(points[0]);
    for (const auto& p : points) bbox.encapsulate(p);
    m_center = bbox.centroid();

    // --- Find Seed Triangle ---
    constexpr float maxDist = std::numeric_limits<float>::max();
    std::size_t i0 = INVALID_INDEX, i1 = INVALID_INDEX, i2 = INVALID_INDEX;

    // Find i0: closest to center
    float min_dist = maxDist;
    for (std::size_t i = 0; i < n; ++i)
    {
        float d = m_center.distanceSquared(points[i]);
        if (d < min_dist) { i0 = i; min_dist = d; }
    }
    Vector2D p0 = points[i0];

    // Find i1: closest to i0
    min_dist = maxDist;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i == i0) continue;
        float d = p0.distanceSquared(points[i]);
        if (d < min_dist && d > 0.0f) { i1 = i; min_dist = d; }
    }
    if (i1 == INVALID_INDEX) throw std::runtime_error("All points are duplicates.");
    Vector2D p1 = points[i1];

    // Find i2: smallest circumcircle with p0 and p1
    float min_radius = maxDist;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (i == i0 || i == i1) continue;
        float r = circumradiusSquared(p0, p1, points[i]);
        if (r < min_radius) { i2 = i; min_radius = r; }
    }
    if (i2 == INVALID_INDEX) throw std::runtime_error("All points are collinear.");
    Vector2D p2 = points[i2];

    // Ensure CCW ordering
    if (isCCW(p0, p1, p2)) std::swap(i1, i2);

    return { i0, i1, i2 };
}

void Delaunator::init_hull(const std::vector<Vector2D>& points,
    std::size_t i0, std::size_t i1, std::size_t i2)
{
    const std::size_t n = points.size();
    constexpr float GOLDEN_RATIO = 1.618f;
    m_hash_size = static_cast<std::size_t>(GOLDEN_RATIO * std::ceil(std::sqrt(n)));
    m_hash.assign(m_hash_size, INVALID_INDEX);

    hull_prev.resize(n);
    hull_next.resize(n);
    hull_tri.resize(n);

    hull_start = i0;

    hull_next[i0] = hull_prev[i2] = i1;
    hull_next[i1] = hull_prev[i0] = i2;
    hull_next[i2] = hull_prev[i1] = i0;

    hull_tri[i0] = 0;
    hull_tri[i1] = 1;
    hull_tri[i2] = 2;

    m_hash[hash_key(points[i0])] = i0;
    m_hash[hash_key(points[i1])] = i1;
    m_hash[hash_key(points[i2])] = i2;

    std::size_t max_triangles = n < 3 ? 1 : 2 * n - 5;
    triangles.reserve(max_triangles * 3);
    halfedges.reserve(max_triangles * 3);

    add_triangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
}

void Delaunator::insert_point_into_hull(const std::vector<Vector2D>& points,
    size_t i, size_t start)
{
    size_t e = start;
    size_t q;

    // Advance until we find a visible edge
    while (true)
    {
        q = hull_next[e];
        if (points[i] == points[e] || points[i] == points[q]) { e = INVALID_INDEX; break; }
        if (isCCW(points[i], points[e], points[q])) break;
        e = q;
        if (e == start) { e = INVALID_INDEX; break; }
    }
    if (e == INVALID_INDEX) return; //near-duplicate

    // First triangle
    size_t t = add_triangle(e, i, hull_next[e], INVALID_INDEX, INVALID_INDEX, hull_tri[e]);
    hull_tri[i] = legalize(points, t + 2);
    hull_tri[e] = t;

    // Forward walk
    size_t next = hull_next[e];
    while (true)
    {
        q = hull_next[next];
        if (!isCCW(points[i], points[next], points[q])) break;
        t = add_triangle(next, i, q, hull_tri[i], INVALID_INDEX, hull_tri[next]);
        hull_tri[i] = legalize(points, t + 2);
        hull_next[next] = next; // mark as removed
        next = q;
    }

    // Backward walk
    if (e == start) {
        while (true)
        {
            q = hull_prev[e];
            if (!isCCW(points[i], points[q], points[e])) break;
            t = add_triangle(q, i, e, INVALID_INDEX, hull_tri[e], hull_tri[q]);
            legalize(points, t + 2);
            hull_tri[q] = t;
            hull_next[e] = e; // mark as removed
            e = q;
        }
    }

    // Update hull
    hull_prev[i] = e;
    hull_start = e;
    hull_prev[next] = i;
    hull_next[e] = i;
    hull_next[i] = next;

    m_hash[hash_key(points[i])] = i;
    m_hash[hash_key(points[e])] = e;
}

void Delaunator::insert_points(const std::vector<Vector2D>& points,
    const std::vector<std::size_t>& sorted_ids,
    std::size_t i0, std::size_t i1, std::size_t i2)
{
    Vector2D prevPoint(std::numeric_limits<float>::quiet_NaN(),
                       std::numeric_limits<float>::quiet_NaN());

    for (size_t k = 0; k < points.size(); k++)
    {
        const size_t i = sorted_ids[k];
        const Vector2D curPoint = points[i];

        if (k > 0 && prevPoint == curPoint) continue;
        prevPoint = curPoint;

        if (i == i0 || i == i1 || i == i2) continue;

        size_t start = find_hull_start(curPoint);
        assert(start != INVALID_INDEX);
        insert_point_into_hull(points, i, start);
    }
}

Delaunator::Delaunator(const std::vector<Vector2D>& inputPoints)
{
    // --- Validation ---
    const std::size_t n = inputPoints.size();
    assert(n >= 3);

    std::vector<std::size_t> m_ids(n);
    m_ids.resize(inputPoints.size());
    std::iota(m_ids.begin(), m_ids.end(), 0);

    auto [i0, i1, i2] = find_seed_triangle(inputPoints);

    float circumRadius;
    circumcircleSquared(inputPoints[i0], inputPoints[i1], inputPoints[i2], m_center, circumRadius);

    // --- Sort dists ---
    std::vector<float> m_dists;
    m_dists.reserve(inputPoints.size());
    for (const Vector2D& p : inputPoints)
    {
        m_dists.push_back(p.distanceSquared(m_center));
    }

    std::sort(m_ids.begin(), m_ids.end(),
        [&m_dists](std::size_t i, std::size_t j) {
            return m_dists[i] < m_dists[j];
        });

    init_hull(inputPoints, i0, i1, i2);

    insert_points(inputPoints, m_ids, i0, i1, i2);
}

std::size_t Delaunator::legalize(const std::vector<Vector2D>& inputPoints, std::size_t a)
{
    std::size_t stackIndex = 0;
    std::size_t currentEdgeIndex = a;
    std::size_t oppositeRightEdgeIndex = 0;
    m_edge_stack.clear();

    while (true)
    {
        const size_t oppositeEdgeIndex = halfedges[currentEdgeIndex];

        const size_t currentTriangleStart = 3 * (currentEdgeIndex / 3);
        const size_t currentRightEdgeIndex = currentTriangleStart + (currentEdgeIndex + 2) % 3;
        oppositeRightEdgeIndex = currentRightEdgeIndex;

        // If the current edge is on the boundary (no opposite triangle)
        // => Pop from the stack to continue.
        if (oppositeEdgeIndex == INVALID_INDEX)
        {
            if (stackIndex > 0)
            {
                stackIndex--;
                currentEdgeIndex = m_edge_stack[stackIndex];
                continue;
            }
            else
            {
                break;
            }
        }

        // Get the four points involved in a potential flip.
        // p0 and p1 are the vertices of the opposing triangles.
        // pl and pr are the vertices of the shared edge.
        const size_t oppositeTriangleStart = 3 * (oppositeEdgeIndex / 3);
        const size_t currentLeftEdgeIndex  = currentTriangleStart  + (currentEdgeIndex + 1)  % 3;
        const size_t oppositeLeftEdgeIndex = oppositeTriangleStart + (oppositeEdgeIndex + 2) % 3;

        const std::size_t p0_vertex = triangles[currentRightEdgeIndex];
        const std::size_t pr_vertex = triangles[currentEdgeIndex];
        const std::size_t pl_vertex = triangles[currentLeftEdgeIndex];
        const std::size_t p1_vertex = triangles[oppositeLeftEdgeIndex];

        // Check if the Delaunay condition is violated.
        const bool violatesDelaunay = isPointInCircumcircle(
            inputPoints[p0_vertex], inputPoints[p1_vertex], inputPoints[pr_vertex], inputPoints[pl_vertex]
        );

        if (violatesDelaunay)
        {
            // Flip the edge to satisfy the Delaunay condition.
            triangles[currentEdgeIndex] = p1_vertex;
            triangles[oppositeEdgeIndex] = p0_vertex;

            auto oppositeHalfedge = halfedges[oppositeLeftEdgeIndex];

            // An edge on the convex hull was swapped.
            // => We need to fix the hull's halfedge reference.
            if (oppositeHalfedge == INVALID_INDEX)
            {
                std::size_t e = hull_start;
                do {
                    if (hull_tri[e] == oppositeLeftEdgeIndex)
                    {
                        hull_tri[e] = currentEdgeIndex;
                        break;
                    }
                    e = hull_prev[e];
                } while (e != hull_start);
            }

            link(currentEdgeIndex, oppositeHalfedge);
            link(oppositeEdgeIndex, halfedges[currentRightEdgeIndex]);
            link(currentRightEdgeIndex, oppositeLeftEdgeIndex);

            // Push the new edge to be checked onto the stack.
            const size_t oppositeRightEdge = oppositeTriangleStart + (oppositeEdgeIndex + 1) % 3;
            if (stackIndex < m_edge_stack.size())
            {
                m_edge_stack[stackIndex] = oppositeRightEdge;
            }
            else
            {
                m_edge_stack.push_back(oppositeRightEdge);
            }
            stackIndex++;
        }
        else
        {
            // No flip needed. Pop the next edge from the stack to continue.
            if (stackIndex > 0)
            {
                stackIndex--;
                currentEdgeIndex = m_edge_stack[stackIndex];
            }
            else
            {
                break;
            }
        }
    }
    return oppositeRightEdgeIndex;
}

std::size_t Delaunator::add_triangle(
    std::size_t i0, std::size_t i1, std::size_t i2,
    std::size_t a, std::size_t b, std::size_t c)
{
    const std::size_t t = triangles.size();
    triangles.push_back(i0);
    triangles.push_back(i1);
    triangles.push_back(i2);
    link(t    , a);
    link(t + 1, b);
    link(t + 2, c);
    return t;
}

void Delaunator::link(const std::size_t a, const std::size_t b)
{
    auto link_one_way = [&](std::size_t from, std::size_t to) 
    {
        if (from == halfedges.size())
            halfedges.push_back(to);
        else if (from < halfedges.size())
            halfedges[from] = to;
        else
            throw std::runtime_error("Cannot link edge: index out of range");
    };

    link_one_way(a, b);
    if (b != INVALID_INDEX) {
        link_one_way(b, a);
    }
}

TriangleMesh2D fastDelaunayTriangulation(const std::vector<Vector2D>& points)
{
    if (points.size() < 3)
        return TriangleMesh2D(points, {});

    Delaunator dt(points);
    return dt.getTriangleMesh(points);
}


} // namespace Math

} // namespace Arns