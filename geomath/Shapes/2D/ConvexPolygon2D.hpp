/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "Polygon2D.hpp"

#include <math.h>
#include <optional>

namespace Arns
{

namespace geomath
{

class ConvexPolygon2D : public IPolygonalShape2D //: public Polygon2D
{
public:
    // --- Metadata & Interface Queries ---

    ShapeType2D type() const { return SHAPE2D_CONVEX_POLYGON; }
    static constexpr ShapeType2D shapeType = SHAPE2D_CONVEX_POLYGON;    

    // --- Data Accessors & Mutators ---

    size_t vertexCount() const override
    {
        return m_vertices.size();
    }

    Vector2D vertexAt(size_t index) const
    {
        return m_vertices[index];
    }

    const std::span<const Vector2D> vertices() const override
    {
        return m_vertices;
    }

    using iterator = std::vector<Vector2D>::iterator;
    using const_iterator = std::vector<Vector2D>::const_iterator;

    iterator begin() { return m_vertices.begin(); }
    iterator end() { return m_vertices.end(); }

    const_iterator begin() const { return m_vertices.begin(); }
    const_iterator end() const { return m_vertices.end(); }

    // --- Geometric Properties ---

    bool isConvex() const noexcept
    {
        return true;
    }

    real_t area() const override
    {
        double area = 0;
        for (int i = 0; i < m_vertices.size(); i++)
        {
            const Vector2D& a = m_vertices[i];
            const Vector2D& b = m_vertices[(i + 1 == m_vertices.size()) ? 0 : i + 1];
            area += static_cast<double>(a.x) * static_cast<double>(b.y)
              - static_cast<double>(b.x) * static_cast<double>(a.y);
        }

        return (real_t)area / real_t{2};
    }

    real_t perimeter() const override
    {
        real_t perimeter = 0;

        for (int i = 0; i < m_vertices.size(); i++)
        {
            perimeter += (m_vertices[i] - m_vertices[(i + 1) % m_vertices.size()]).length();
        }

        return perimeter;
    }

    Vector2D centroid() const override
    {
        Vector2D centroid;

        for (int i = 0; i < m_vertices.size(); i++)
        {
            centroid += m_vertices[i];
        }

        return centroid / m_vertices.size();
    }

    // --- Transform / modification ---

    ConvexPolygon2D& translate(const Vector2D &translation) override
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i] += translation;
        return *this;
    }

    ConvexPolygon2D& rotate(real_t angle)
    {
        const Vector2D centroid = this->centroid();
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i].rotateAround(angle, centroid);
        return *this;
    }

    ConvexPolygon2D& rotate(real_t angle, const Vector2D& point)
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i].rotateAround(angle, point);
        return *this;
    }

    ConvexPolygon2D& scale(real_t factor, const Vector2D& point)
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i] = point + (m_vertices[i] - point) * factor;
        return *this;
    }

    ConvexPolygon2D& scale(real_t factor)
    {
        return scale(factor, this->centroid());
    }

    ConvexPolygon2D& transform(const Matrix3x3& matrix)
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i] = matrix.transformPoint(m_vertices[i]);
        return *this;
    }

    // --- Lifecycle / Factory Methods ---

    ConvexPolygon2D copy() const
    {
        return *this;
    }

    std::unique_ptr<IFiniteShape2D> clone() const override
    {
        return std::make_unique<ConvexPolygon2D>(*this);
    }

    // --- Spatial Queries ---

    bool contains(const Vector2D &point) const override
    {
        return detail::isPointInsidePolygon(point, std::span<const Vector2D>(m_vertices));
    }

    bool contains(const Segment2D& segment) const
    {
        return isSegmentInsideShape(segment, *this);
    }

    template <typename T> requires FiniteShape2D<T>
    bool contains(const T& innerShape) const
    {
        return isShapeInsideShape(innerShape, *this);
    }

    BBox2D boundingBox() const override
    {
        if (m_vertices.empty())
            return BBox2D();

        Vector2D min = m_vertices[0];
        Vector2D max = m_vertices[0];

        for (int i = 1; i < m_vertices.size(); i++)
        {
            min = Vector2D::min(min, m_vertices[i]);
            max = Vector2D::max(max, m_vertices[i]);
        }

        return BBox2D(min, max);
    }

    // --- Operators ---

    /*
    Vector2D& operator[](size_t index)
    {
        return m_vertices[index];
    }
    */

    const Vector2D& operator[](size_t index) const override
    {
        return m_vertices[index];
    }

public:
    static ConvexPolygon2D fromPoints(const std::vector<Vector2D>& points)
    {
        return ConvexPolygon2D(compute_convex_hull(points));
    }

    // Danger: Assumes whatever is supplied is already convex
    static ConvexPolygon2D createUnchecked(std::vector<Vector2D> vertices)
    {
        return ConvexPolygon2D(std::move(vertices));
    }

    static ConvexPolygon2D fromHull(const std::vector<Vector2D>& vertices)
    {
        return ConvexPolygon2D(compute_convex_hull(vertices));
    }

    static std::optional<ConvexPolygon2D> tryCreate(std::vector<Vector2D> vertices)
    {
        if (!is_convex_sequence(vertices)) return std::nullopt;
        return ConvexPolygon2D(std::move(vertices));
    }

private:
    std::vector<Vector2D> m_vertices;

    // --- Constructors ---

    explicit ConvexPolygon2D(std::vector<Vector2D> vertices) : m_vertices(std::move(vertices)) {}

    static std::vector<Vector2D> compute_convex_hull(const std::vector<Vector2D>& v)
    {
        auto points = v;

        sort(points.begin(), points.end(), [](const Vector2D& a, const Vector2D& b)
            {
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            });

        points.erase(std::unique(points.begin(), points.end(), [](Vector2D a, Vector2D b)
            { return a == b; }), points.end());

        size_t hullSize = 0;
        size_t n = points.size();

        if (n <= 1)
            return points;
        
	    std::vector<Vector2D> hull(n*2);

	    for (size_t i = 0; i < n; ++i)
        {
	    	while (hullSize >= 2 && orient2D(hull[hullSize-2], hull[hullSize-1], points[i]) <= 0)
                hullSize--;
	    	hull[hullSize++] = points[i];
	    }

	    for (size_t i = n-1, t = hullSize+1; i > 0; --i)
        {
	    	while (hullSize >= t && orient2D(hull[hullSize-2], hull[hullSize-1], points[i-1]) <= 0)
                hullSize--;
	    	hull[hullSize++] = points[i-1];
	    }

	    hull.resize(hullSize-1);
	    return hull;
    }

    static bool is_convex_sequence(std::span<const Vector2D> vertices)
    {
        const size_t n = vertices.size();
        if (n < 3) 
            return false;

        int sign = 0;

        for (size_t i = 0; i < n; ++i)
        {
            const Vector2D& a = vertices[i];
            const Vector2D& b = vertices[(i + 1) % n];
            const Vector2D& c = vertices[(i + 2) % n];

            // Compute edge vectors
            const Vector2D ab = b - a;
            const Vector2D bc = c - b;

            real_t cross = ab.cross(bc);

            if (approximatelyZero(cross))
            {
                // Collinear points
                continue;
            }

            int currentSign = (cross > 0) ? 1 : -1;

            if (sign == 0)
            {
                // First turn
                sign = currentSign;
            }
            else if (currentSign != sign)
            {
                // Turning direction changes -> concave
                return false;
            }
        }
        return true;
    }

    // --- Comparison Operators ---

    bool operator==(const ConvexPolygon2D& other) const
    {
        if (this->vertexCount() != other.vertexCount())
            return false;

        for (int i = 0; i < m_vertices.size(); i++)
        {
            if (this->vertexAt(i) != other.vertexAt(i))
                return false;
        }

        return true;
    }

    bool operator!=(const ConvexPolygon2D& other) const
    {
        return !(*this == other);
    }
};

static_assert(FiniteShape2D<Polygon2D>);

} // namespace Math

} // namespace Arns