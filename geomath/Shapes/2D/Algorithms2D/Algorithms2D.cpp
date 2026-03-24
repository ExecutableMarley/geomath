/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms2D.hpp"

#include "CommonMath.hpp"
#include "../Ray2D.hpp"
#include "../Line2D.hpp"
#include "../BBox2D.hpp"
#include "../IShape2D.hpp"
#include "../Circle2D.hpp"
#include "../Triangle2D.hpp"
#include "../Rectangle2D.hpp"
#include "../Polygon2D.hpp"
#include "../ConvexPolygon2D.hpp"

#include <cassert>
#include <span>

namespace Arns
{

namespace Math
{

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd)
{
    Vector2D r = segmentEnd - segmentStart;
    Vector2D s = point - segmentStart;
    real_t cross = r.cross(s);
    if (!approximatelyZero(cross))
        return false;

    real_t dot = r.dot(s);
    if (dot < 0.0 || dot > r.lengthSquared())
        return false;

    return true;
}

bool isPointOnSegment(const Vector2D& point, const Line2D& line)
{
    return isPointOnSegment(point, line.m_start, line.m_end);
}

bool isSegmentOnSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2)
{
    return isPointOnSegment(s1, k1, k2) && isPointOnSegment(s2, k1, k2);
}

bool isSegmentOnSegment(const Segment2D& segment1, const Segment2D& segment2)
{
    return isSegmentOnSegment(segment1.m_start, segment1.m_end, segment2.m_start, segment2.m_end);
}

// Point inside shape

bool isPointInsideBBox(const Vector2D& point, const BBox2D& bbox)
{
    return bbox.contains(point);
}

bool isPointInsideTriangle(const Vector2D& point, const Triangle2D& triangle)
{
    return triangle.contains(point);
}

bool isPointInsideRectangle(const Vector2D& point, const Rectangle2D& rectangle)
{
    return rectangle.contains(point);
}

//Todo: Consider the O(log n) method
bool isPointInsideConvexPolygon(const Vector2D& point, const std::span<const Vector2D> vertices)
{
    const size_t n = vertices.size();
    if (n < 3) return false;

    bool hasPositive = false;
    bool hasNegative = false;

    for (size_t i = 0; i < n; ++i)
    {
        const auto& A = vertices[i];
        const auto& B = vertices[(i + 1) % n];

        real_t val = orient2D(A, B, point);

        if (approximatelyGreaterAbs(val, real_t{0}))   hasPositive = true;
        else if (approximatelyLessAbs(val, real_t{0})) hasNegative = true;

        if (hasPositive && hasNegative) return false;
    }

    return true; 
}

//Non self intersecting polygons
bool isPointInsidePolygon(const Vector2D& point, const std::span<const Vector2D> vertices)
{
    bool contains = false;
    for (size_t i = 0, j = vertices.size() - 1; i < vertices.size(); j = i++)
    {
        if (isPointOnSegment(point, vertices[i], vertices[j]))
            return true;

        if (((vertices[i].y > point.y) != (vertices[j].y > point.y)) &&
            (point.x < (vertices[j].x - vertices[i].x) * (point.y - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x))
        {
            contains = !contains;
        }
    }

    return contains;
}

bool isPointInsideConvexPolygon(const Vector2D& point, const ConvexPolygon2D& polygon)
{
    return isPointInsideConvexPolygon(point, polygon.vertices());
}

bool isPointInsidePolygon(const Vector2D& point, const Polygon2D& polygon)
{
    return polygon.contains(point);
}

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle)
{
    return circle.contains(point);
}

// Segment inside shape

bool isSegmentInsideBBox(const Segment2D& segment, const BBox2D& bbox)
{
    return isPointInsideBBox(segment.m_start, bbox) && isPointInsideBBox(segment.m_end, bbox);
}

bool isSegmentInsideTriangle(const Segment2D& segment, const Triangle2D& triangle)
{
    return isPointInsideTriangle(segment.m_start, triangle) && isPointInsideTriangle(segment.m_end, triangle);
}

bool isSegmentInsideRectangle(const Segment2D& segment, const Rectangle2D& rectangle)
{
    return isPointInsideRectangle(segment.m_start, rectangle) && isPointInsideRectangle(segment.m_end, rectangle);
}

bool isSegmentInsideConvexPolygon(const Segment2D& segment, const ConvexPolygon2D& polygon)
{
    return isPointInsideConvexPolygon(segment.m_start, polygon) && isPointInsideConvexPolygon(segment.m_end, polygon);
}

bool isSegmentInsidePolygon(const Segment2D& segment, const Polygon2D& polygon)
{
    if (isPointInsidePolygon(segment.m_start, polygon) && isPointInsidePolygon(segment.m_end, polygon))
    {
        if (intersect(segment, polygon))
            return false;
        return true;
    }
    return false;
}

bool isSegmentInsideCircle(const Segment2D& segment, const Circle2D& circle)
{
    return isPointInsideCircle(segment.m_start, circle) && isPointInsideCircle(segment.m_end, circle);
}


//Shape inside BBox

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2)
{
    return bbox2.contains(bbox1);
}

bool isPolygonalInsideBBox(const IPolygonalShape2D& polygon, const BBox2D& bbox)
{
    for (const auto& v : polygon.vertices())
        if (!bbox.contains(v))
            return false;
    return true;
}

bool isCircleInsideBBox(const Circle2D& circle, const BBox2D& bbox)
{
    return bbox.contains(circle.m_center) && bbox.minDistanceSquared(circle.m_center) < circle.m_radius * circle.m_radius;
}

//

//Helper function
bool edgeIntersectionPolygonWithPolygon(const std::span<const Vector2D> vertices1, const std::span<const Vector2D> vertices2, HitInfo2D* hitInfo = nullptr)
{
    const size_t countA = vertices1.size();
    const size_t countB = vertices2.size();

    if (countA < 3 || countB < 3)
        return false;

    for (size_t i = 0; i < countA; ++i)
    {
        const Vector2D& a1 = vertices1[i];
        const Vector2D& a2 = vertices1[(i + 1) % countA];

        for (size_t j = 0; j < countB; ++j)
        {
            const Vector2D& b1 = vertices2[j];
            const Vector2D& b2 = vertices2[(j + 1) % countB];

            if (intersectSegmentWithSegment(a1, a2, b1, b2, hitInfo))
                return true;
        }
    }
    return false;
}

//

bool isBBoxInsideConvexPolygon(const BBox2D& bbox, const std::span<const Vector2D> convexPoly)
{
    return isPointInsideConvexPolygon(bbox.bottomLeft(), convexPoly) && isPointInsideConvexPolygon(bbox.topLeft(), convexPoly) && 
        isPointInsideConvexPolygon(bbox.topRight(), convexPoly) && isPointInsideConvexPolygon(bbox.bottomRight(), convexPoly);
}

bool isPolygonInsideConvexPolygon(const std::span<const Vector2D> poly, const std::span<const Vector2D> convexPoly)
{
    for (const auto& point : poly)
    {
        if (!isPointInsideConvexPolygon(point, convexPoly))
            return false;
    }
    return true;
}

bool isCircleInsideConvexPolygon(const Circle2D& circle, const std::span<const Vector2D> convexPoly)
{
    if (!isPointInsideConvexPolygon(circle.m_center, convexPoly))
    {
        return false;
    }

    const size_t n = convexPoly.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Vector2D& v1 = convexPoly[i];
        const Vector2D& v2 = convexPoly[(i + 1) % n];

        real_t dist = distancePointToSegment(circle.m_center, v1, v2);
        if (dist < circle.m_radius)
        {
            return false;
        }
    }
    return true;
}

bool isShapeInsideConvexPolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> convexPoly)
{
    if (isPolygonalShape(shape.type()))
    {
        return isPolygonInsideConvexPolygon(shape.polygonal()->vertices(), convexPoly);
    }
    else if (shape.type() == SHAPE2D_CIRCLE)
    {
        return isCircleInsideConvexPolygon(*shape.shape_cast<Circle2D>(), convexPoly);
    }
    assert(false && "Unreachable code reached!");
    return false;
}

bool isBBoxInsidePolygon(const BBox2D& bbox, const std::span<const Vector2D> poly)
{
    std::vector<Vector2D> bboxPoly = {bbox.bottomLeft(), bbox.bottomRight(), bbox.topRight(), bbox.topLeft()};

    if (!isPointInsidePolygon(bbox.bottomLeft(), poly) || !isPointInsidePolygon(bbox.bottomRight(), poly) ||
        !isPointInsidePolygon(bbox.topRight(), poly)   || !isPointInsidePolygon(bbox.topLeft(), poly))
            return false;

    if (edgeIntersectionPolygonWithPolygon(bboxPoly, poly))
        return false;

    return true;
}

bool isPolygonInsidePolygon(const std::span<const Vector2D> poly1, const std::span<const Vector2D> poly2)
{
    for (const auto& point : poly2)
    {
        if (!isPointInsidePolygon(point, poly2))
        {
            return false;
        }
    }

    if (edgeIntersectionPolygonWithPolygon(poly1, poly2))
        return false;

    return true;
}

//Same as isCircleInsideConvexPolygon. Should be merged to one function
bool isCircleInsidePolygon(const Circle2D& circle, const std::span<const Vector2D> poly)
{
    if (!isPointInsideConvexPolygon(circle.m_center, poly))
    {
        return false;
    }

    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Vector2D& v1 = poly[i];
        const Vector2D& v2 = poly[(i + 1) % n];

        real_t dist = distancePointToSegment(circle.m_center, v1, v2);
        if (dist < circle.m_radius)
        {
            return false;
        }
    }
    return true;
}

bool isShapeInsidePolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> poly)
{
    if (isPolygonalShape(shape.type()))
    {
        return isPolygonInsidePolygon(shape.polygonal()->vertices(), poly);
    }
    else if (shape.type() == SHAPE2D_CONVEX_POLYGON)
    {
        return isCircleInsidePolygon(*shape.shape_cast<Circle2D>(), poly);
    }
    assert(false && "Unreachable code reached!");
    return false;
}

bool isBBoxInsideCircle(const BBox2D& bbox, const Circle2D& circle)
{
    return bbox.maxDistanceSquared(circle.m_center) <= circle.m_radius * circle.m_radius;
}

bool isPolygonInsideCircle(const std::span<const Vector2D> poly, const Circle2D& circle)
{
    for (const auto& point : poly)
    {
        if (!isPointInsideCircle(point, circle))
            return false;
    }
    return true;
}

bool isCircleInsideCircle(const Circle2D& circle1, const Circle2D& circle2)
{
    return circle1.contains(circle2);
}

bool isShapeInsideCircle(const IFiniteShape2D& shape, const Circle2D& circle)
{
    if (isPolygonalShape(shape.type()))
    {
        return isPolygonInsideCircle(shape.polygonal()->vertices(), circle);
    }
    else if (shape.type() == SHAPE2D_CONVEX_POLYGON)
    {
        return isCircleInsideCircle(*shape.shape_cast<Circle2D>(), circle);
    }
    assert(false && "Unreachable code reached!");
    return false;
}


bool isShapeInsideShape(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2)
{
    if (isConvexPolygonal(shape2.type()))
    {
        return isShapeInsideConvexPolygon(shape1, shape2.polygonal()->vertices());
    }
    else if (isPolygonalShape(shape2.type()))
    {
        return isShapeInsidePolygon(shape1, shape2.polygonal()->vertices());
    }
    else if (shape2.type() == SHAPE2D_CIRCLE)
    {
        return isShapeInsideCircle(shape1, *shape2.shape_cast<Circle2D>());
    }
    assert(false && "Unreachable code reached!");
    return false;
}

//Consider throwing an exception after unreachable code asserts




// Distance calculation algorithms

real_t distancePointToSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd, Vector2D* closestPoint)
{
    const Vector2D r = segmentEnd - segmentStart;
    const Vector2D s = point - segmentStart;
    const real_t rDotR = r.dot(r);
    if (rDotR == 0.0f)
    {
        if (closestPoint) *closestPoint = segmentStart;
        return (point - segmentStart).length();
    }
    real_t t = std::clamp(s.dot(r) / rDotR, 0.0f, 1.0f);
    const Vector2D projection = segmentStart + r * t;
    if (closestPoint)
        *closestPoint = projection;
    return (point - projection).length();
}

real_t distance(const Vector2D& point, const Segment2D& segment, Vector2D* closestPoint)
{
    return distancePointToSegment(point, segment.m_start, segment.m_end, closestPoint);
}

real_t distanceSegmentToSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const Vector2D delta1 = s2 - s1;
    const Vector2D delta2 = k2 - k1;
    real_t lengthSquared1 = delta1.lengthSquared();
    real_t lengthSquared2 = delta2.lengthSquared();

    if (lengthSquared1 == 0 || lengthSquared2 == 0)
    {
        if (closestPoint1)
            *closestPoint1 = s1;
        if (closestPoint2)
            *closestPoint2 = k1;
        return (s1 - k1).length();
    }

    const Vector2D startDelta = s1 - k1;
    
    const real_t a = delta1.dot(startDelta);
    const real_t b = delta2.dot(startDelta);
    const real_t c = delta1.dot(delta2);

    const real_t d = lengthSquared1 * lengthSquared2 - c * c;

    real_t t = d != 0 ? (b * c - a * lengthSquared2) / d : 0.0f; 
    real_t u = d != 0 ? (lengthSquared1 * b - c * a) / d : 0.0f;

    t = clamp(t, 0.0f, 1.0f);
    u = clamp(u, 0.0f, 1.0f);

    const Vector2D closestPointOnLine1 = s1 + delta1 * t;
    const Vector2D closestPointOnLine2 = k1 + delta2 * u;

    if (closestPoint1)
        *closestPoint1 = closestPointOnLine1;
    if (closestPoint2)
        *closestPoint2 = closestPointOnLine2;
    
    return (closestPointOnLine1 - closestPointOnLine2).length();
}

real_t distance(const Line2D& segment1, const Line2D& segment2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distanceSegmentToSegment(segment1.m_start, segment1.m_end, segment2.m_start, segment2.m_end, closestPoint1, closestPoint2);
}

//Todo: Distance to shapes

real_t distance(const Vector2D& point, const BBox2D& bbox, Vector2D* closestPoint)
{
    Vector2D c = bbox.closestPoint(point);
    if (closestPoint) *closestPoint = c;
    return point.distance(c);
}

//Todo: Consider renaming slightly

//Helper
real_t distancePointToPolygonVertices(const Vector2D& point, const std::span<const Vector2D> polygon, Vector2D* closestPoint)
{
    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D best;

    for (size_t i = 0; i < polygon.size(); ++i)
    {
        const Vector2D& a = polygon[i];
        const Vector2D& b = polygon[(i + 1) % polygon.size()];

        Vector2D q;
        real_t dist = distancePointToSegment(point, a, b, &q);

        if (dist < minDist)
        {
            minDist = dist;
            best = q;
        }
    }

    if (closestPoint) *closestPoint = best;
    return std::sqrt(minDist);
}

/*
real_t distanceToConvexPolygon(const Vector2D& point, const std::span<const Vector2D> polygon, Vector2D* closestPoint)
{
    if (isPointInsideConvexPolygon(point, polygon))
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }
    return distancePointToPolygonVertices(point, polygon, closestPoint);
}*/ 

real_t distance(const Vector2D& point, const IPolygonalShape2D& polygon, Vector2D* closestPoint)
{
    if (isConvexPolygonal(polygon.type()))
    {
        if (isPointInsideConvexPolygon(point, polygon.vertices()))
        {
            if (closestPoint)
                *closestPoint = point;
            return real_t{0};
        }
        return distancePointToPolygonVertices(point, polygon.vertices(), closestPoint);
    }
    else
    {
        if (isPointInsidePolygon(point, polygon.vertices()))
        {
            if (closestPoint)
                *closestPoint = point;
            return real_t{0};
        }
        return distancePointToPolygonVertices(point, polygon.vertices(), closestPoint);
    }
}

real_t distance(const Vector2D& point, const Circle2D& circle, Vector2D* closestPoint)
{
    Vector2D delta = point - circle.m_center;
    real_t dist = delta.length();

    if (dist <= circle.m_radius)
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }

    //real_t dist = std::sqrt(distSq);
    Vector2D dir = delta / dist;

    if (closestPoint)
        *closestPoint = circle.m_center + dir * circle.m_radius;

    return std::abs(dist - circle.m_radius);
}

real_t distanceSegmentToEdgeList(const Line2D& seg, std::span<const Vector2D> vertices, bool closed, Vector2D* closestSeg, Vector2D* closestEdge)
{
    const size_t n = vertices.size();
    if (n < 2)
        return std::numeric_limits<real_t>::max();

    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D bestSeg;
    Vector2D bestEdge;

    const size_t edgeCount = closed ? n : n - 1;

    for (size_t i = 0; i < edgeCount; ++i)
    {
        const Vector2D& a = vertices[i];
        const Vector2D& b = vertices[(i + 1) % n];

        Vector2D cs, ce;
        real_t d = distanceSegmentToSegment(seg.m_start, seg.m_end, a, b, &cs, &ce);

        if (d < minDist)
        {
            minDist  = d;
            bestSeg  = cs;
            bestEdge = ce;

            // => Intersection
            if (minDist == real_t(0))
                break;
        }
    }

    if (closestSeg)  *closestSeg  = bestSeg;
    if (closestEdge) *closestEdge = bestEdge;

    return minDist;
}

real_t distanceSegmentToBBox(const Segment2D& segment, const BBox2D& bbox, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distanceSegmentToEdgeList(segment, bbox.corners(), true, closestPoint1, closestPoint2);
}

real_t distance(const Segment2D& segment, const IPolygonalShape2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distanceSegmentToEdgeList(segment, polygon.vertices(), true, closestPoint1, closestPoint2);
}

real_t distanceSegmentToCircle(const Segment2D& segment, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    Vector2D p;
    real_t dist = distancePointToSegment(circle.m_center, segment.m_start, segment.m_end, &p);

    //Trivial case
    if (dist <= circle.m_radius)
    {
        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = p;
        return real_t(0);
    }

    //Not intersecting case
    Vector2D dir = p - circle.m_center;
    real_t len = dir.length();

    if (len > real_t{0})
        dir /= len;

    Vector2D q = circle.m_center + dir * circle.m_radius;

    if (closestPoint1) *closestPoint1 = p;
    if (closestPoint2) *closestPoint2 = q;

    return len - circle.m_radius;
}


real_t distance(const BBox2D& a, const BBox2D& b, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    real_t dx = 0;
    real_t dy = 0;

    Vector2D pa;
    Vector2D pb;

    // --- X axis ---
    if (a.m_max.x < b.m_min.x)
    {
        dx = b.m_min.x - a.m_max.x;
        pa.x = a.m_max.x;
        pb.x = b.m_min.x;
    }
    else if (b.m_max.x < a.m_min.x)
    {
        dx = a.m_min.x - b.m_max.x;
        pa.x = a.m_min.x;
        pb.x = b.m_max.x;
    }
    else
    {
        // overlap
        real_t x = std::max(a.m_min.x, b.m_min.x);
        pa.x = pb.x = x;
    }

    // --- Y axis ---
    if (a.m_max.y < b.m_min.y)
    {
        dy = b.m_min.y - a.m_max.y;
        pa.y = a.m_max.y;
        pb.y = b.m_min.y;
    }
    else if (b.m_max.y < a.m_min.y)
    {
        dy = a.m_min.y - b.m_max.y;
        pa.y = a.m_min.y;
        pb.y = b.m_max.y;
    }
    else
    {
        // overlap
        real_t y = std::max(a.m_min.y, b.m_min.y);
        pa.y = pb.y = y;
    }

    if (closestPoint1) *closestPoint1 = pa;
    if (closestPoint2) *closestPoint2 = pb;

    return std::sqrt(dx * dx + dy * dy);
}

//real_t distanceCircleToCircle(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const BBox2D& bbox1, const IPolygonalShape2D& poly2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    //Todo: This can (and should) be optimized

    const auto v = poly2.vertices();
    const size_t n = v.size();

    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D bestP1, bestP2;
    
    auto b = bbox1.corners();

    for (int i = 0; i < 4; ++i)
    {
        Vector2D a1 = b[i];
        Vector2D b1 = b[(i + 1) % 4];

        for (int j = 0; j < n; ++j)
        {
            Vector2D a2 = v[j];
            Vector2D b2 = v[(j + 1) % n];

            Vector2D p1, p2;
            real_t d = distanceSegmentToSegment(a1, b1, a2, b2, &p1, &p2);

            if (d < minDist)
            {
                minDist = d;
                bestP1 = p1;
                bestP2 = p2;

                if (d == real_t(0))
                {
                    if (closestPoint1) *closestPoint1 = bestP1;
                    if (closestPoint2) *closestPoint2 = bestP2;
                    return real_t(0);
                }
            }
        }
    }

    if (closestPoint1) *closestPoint1 = bestP1;
    if (closestPoint2) *closestPoint2 = bestP2;

    return minDist;
}

real_t distance(const BBox2D& bbox, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const Vector2D& c = circle.m_center;
    const real_t r = circle.m_radius;

    Vector2D p;
    p.x = std::clamp(c.x, bbox.m_min.x, bbox.m_max.x);
    p.y = std::clamp(c.y, bbox.m_min.y, bbox.m_max.y);

    Vector2D delta = c - p;
    real_t distSq = delta.lengthSquared();

    if (distSq == real_t(0))
    {
        // Find closest point on box boundary instead
        real_t left   = c.x - bbox.m_min.x;
        real_t right  = bbox.m_max.x - c.x;
        real_t bottom = c.y - bbox.m_min.y;
        real_t top    = bbox.m_max.y - c.y;

        // Pick smallest penetration direction
        real_t minDist = left;
        p = { bbox.m_min.x, c.y };

        if (right < minDist) { minDist = right;  p = { bbox.m_max.x, c.y }; }
        if (bottom < minDist){ minDist = bottom; p = { c.x, bbox.m_min.y }; }
        if (top < minDist)   { minDist = top;    p = { c.x, bbox.m_max.y }; }

        delta = c - p;
        real_t dist = delta.length();

        Vector2D dir = (dist > real_t(0)) ? (delta / dist) : Vector2D{1, 0};

        if (minDist <= r)
        {
            Vector2D cp = c - dir * r;

            if (closestPoint1) *closestPoint1 = cp;
            if (closestPoint2) *closestPoint2 = cp;
            return real_t(0);
        }

        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = c - dir * r;

        return minDist - r;
    }

    real_t dist = std::sqrt(distSq);
    Vector2D dir = delta / dist;

    // Overlapping or touching
    if (dist <= r)
    {
        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = p;
        return real_t(0);
    }

    // Separated
    if (closestPoint1) *closestPoint1 = p;
    if (closestPoint2) *closestPoint2 = c - dir * r;

    return dist - r;
}

real_t distance(const IPolygonalShape2D& poly1, const IPolygonalShape2D& poly2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const auto v1 = poly1.vertices();
    const auto v2 = poly2.vertices();

    const int n1 = (int)v1.size();
    const int n2 = (int)v2.size();

    real_t minDist = std::numeric_limits<real_t>::max();

    Vector2D bestP1, bestP2;

    for (int i = 0; i < n1; ++i)
    {
        Vector2D a1 = v1[i];
        Vector2D b1 = v1[(i + 1) % n1];

        for (int j = 0; j < n2; ++j)
        {
            Vector2D a2 = v2[j];
            Vector2D b2 = v2[(j + 1) % n2];

            Vector2D p1, p2;
            real_t d = distanceSegmentToSegment(a1, b1, a2, b2, &p1, &p2);

            if (d < minDist)
            {
                minDist = d;
                bestP1 = p1;
                bestP2 = p2;

                if (minDist == real_t(0))
                {
                    if (closestPoint1) *closestPoint1 = bestP1;
                    if (closestPoint2) *closestPoint2 = bestP2;
                    return real_t(0);
                }
            }
        }
    }

    if (closestPoint1) *closestPoint1 = bestP1;
    if (closestPoint2) *closestPoint2 = bestP2;

    return minDist;
}

real_t distance(const IPolygonalShape2D& poly1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const auto v = poly1.vertices();
    const int n = (int)v.size();

    const Vector2D& center = circle2.m_center;
    const real_t radius = circle2.m_radius;

    real_t minDist = std::numeric_limits<real_t>::max();

    Vector2D bestP1, bestP2;

    // --- Edge vs Circle ---
    for (int i = 0; i < n; ++i)
    {
        Vector2D a = v[i];
        Vector2D b = v[(i + 1) % n];

        Vector2D p1, p2;
        real_t d = distanceSegmentToCircle(Segment2D(a, b), circle2, &p1, &p2);

        if (d < minDist)
        {
            minDist = d;
            bestP1 = p1;
            bestP2 = p2;

            if (minDist == real_t(0))
            {
                if (closestPoint1) *closestPoint1 = bestP1;
                if (closestPoint2) *closestPoint2 = bestP2;
                return real_t(0);
            }
        }
    }

    if (closestPoint1) *closestPoint1 = bestP1;
    if (closestPoint2) *closestPoint2 = bestP2;
    return minDist;
}

real_t distance(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    Vector2D delta = circle2.m_center - circle1.m_center;
    real_t dist = delta.length();

    Vector2D dir = (dist > real_t{0}) ? (delta / dist) : Vector2D{1, 0};

    if (dist <= circle1.m_radius + circle2.m_radius)
    {
        Vector2D p = circle1.m_center + dir * circle1.m_radius;
        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = p;

        return real_t{0};
    }

    if (closestPoint1) *closestPoint1 = circle1.m_center + dir * circle1.m_radius;
    if (closestPoint2) *closestPoint2 = circle2.m_center - dir * circle2.m_radius;

    return dist - circle1.m_radius - circle2.m_radius;
}


// Intersection calculation algorithms

//    --- Rays ---

// Internal Helper
bool intersectParams(const Vector2D& p, const Vector2D& r, const Vector2D& q, const Vector2D& s, 
    real_t& t, real_t& u)
{
    const real_t rxs = r.cross(s);
    const Vector2D qp = q - p;

    if (approximatelyZero(rxs))
        return false; // parallel

    t = qp.cross(s) / rxs;
    u = qp.cross(r) / rxs;
    return true;
}

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    real_t t, u;
    if (!intersectParams(ray.m_origin, ray.m_direction, p1, p2 - p1, t, u))
        return false;

    if (t < t_min || t > t_max || u < 0.f || u > 1.f)
        return false;

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = ray.m_origin + ray.m_direction * t;
    }

    return true;
}

bool intersect(const Ray2D& ray, const Segment2D& segment, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    return intersectRayWithSegment(ray, segment.m_start, segment.m_end, t_min, t_max, hitInfo);
}

bool intersect(const Ray2D& ray, const BBox2D& bbox, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    for (int i = 0; i < 2; i++)
    {
        if (approximatelyZero(ray.m_direction[i]))
        {
            if (ray.m_origin[i] < bbox.m_min[i] || ray.m_origin[i] > bbox.m_max[i])
                return false;
            continue;
        }

        real_t invD = 1.0f / ray.m_direction[i];
        real_t t0 = (bbox.m_min[i] - ray.m_origin[i]) * invD;
        real_t t1 = (bbox.m_max[i] - ray.m_origin[i]) * invD;

        if (invD < 0.0f)
            std::swap(t0, t1);
        
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        if (t_max <= t_min)
            return false;
    }

    if (hitInfo)
    {
        hitInfo->t = t_min;
        hitInfo->intersectionPoint = ray.m_origin + ray.m_direction * t_min;
    }

    return true;
}

bool intersectRayWithPolygonOptimized(const Ray2D& ray, const IPolygonalShape2D& polygon, real_t t_min, real_t t_max)
{
    for (int i = 0; i < polygon.vertexCount(); i++)
    {
        const Vector2D& p1 = polygon[i];
        const Vector2D& p2 = polygon[(i + 1) % polygon.vertexCount()];

        if (intersectRayWithSegment(ray, p1, p2, t_min, t_max, NULL))
        {
            return true;
        }
    }
    return false;
}

bool intersect(const Ray2D& ray, const IPolygonalShape2D& polygon, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    if (hitInfo == nullptr)
        return intersectRayWithPolygonOptimized(ray, polygon, t_min, t_max);

    bool hit = false;
    real_t t = t_max;

    for (int i = 0; i < polygon.vertexCount(); i++)
    {
        const Vector2D& p1 = polygon[i];
        const Vector2D& p2 = polygon[(i + 1) % polygon.vertexCount()];

        if (intersectRayWithSegment(ray, p1, p2, t_min, t, hitInfo))
        {
            hit = true;
            t = hitInfo->t;
        }
    }
    return hit;
}

bool intersect(const Ray2D& ray, const Circle2D& circle, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    const Vector2D oc = ray.m_origin - circle.m_center;
    const real_t a = ray.m_direction.dot(ray.m_direction);
    const real_t b = 2.0f * oc.dot(ray.m_direction);
    const real_t c = oc.dot(oc) - circle.m_radius * circle.m_radius;
    const real_t discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
        return false;

    const real_t sqrtDiscriminant = sqrt(discriminant);
    real_t t = (-b - sqrtDiscriminant) / (2.0f * a);
    if (t < t_min || t > t_max)
    {
        t = (-b + sqrtDiscriminant) / (2.0f * a);
        if (t < t_min || t > t_max)
            return false;
    }

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = ray.pointAt(t);
        //Vector2D normal = (hitInfo->intersectionPoint - circle.m_center).normalize();
    }

    return true;
}

bool intersect(const Ray2D& ray, const IBaseShape2D& shape, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    switch(shape.type())
    {
        case ShapeType2D::SHAPE2D_TRIANGLE:
            return intersect(ray, dynamic_cast<const Triangle2D&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_RECTANGLE:
        case ShapeType2D::SHAPE2D_POLYGON:
            return intersect(ray, dynamic_cast<const IPolygonalShape2D&>(shape), t_min, t_max, hitInfo);
         case ShapeType2D::SHAPE2D_CIRCLE:
            return intersect(ray, dynamic_cast<const Circle2D&>(shape), t_min, t_max, hitInfo);
        default:
            //Todo: we might just want to throw and error here
            assert("Unreachable Code reached!");
            return false;
    }
}

bool intersectSegmentWithSegment(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    return intersectRayWithSegment(Ray2D(p1, p2 - p1), q1, q2, 0.f, 1.f, hitInfo);
}

bool intersect(const Segment2D& segment1, const Segment2D& segment2, HitInfo2D* hitInfo)
{
    return intersectRayWithSegment(Ray2D(segment1.m_start, segment1.deltaVector()), segment2.m_start, segment2.m_end, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    real_t t, u;
    if (!intersectParams(p1, p2 - p1, q1, q2 - q1, t, u))
        return false;

    if (!approximatelyGreater(t, 0.f) || !approximatelyLess(t, 1.f) ||
        !approximatelyGreater(u, 0.f) || !approximatelyLess(u, 1.f))
        return false;

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = p1 + (p2 - p1) * t;
    }

    return true;
}

bool intersect(const Line2D& line, const BBox2D& bbox, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(line.m_start, line.deltaVector()), bbox, 0.f, 1.f, hitInfo);
}

bool intersect(const Line2D& line, const Triangle2D& triangle, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(line.m_start, line.deltaVector()), triangle, 0.f, 1.f, hitInfo);
}

bool intersect(const Line2D& line, const IPolygonalShape2D& polygon, HitInfo2D* hitInfo)
{
    return(intersect(Ray2D(line.m_start, line.direction()), polygon, 0.0f, line.length(), hitInfo));
}

bool intersect(const Line2D& line, const Circle2D& circle, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(line.m_start, line.deltaVector()), circle, 0.f, 1.f, hitInfo);
}

bool intersect(const Segment2D& segment, const IBaseShape2D& shape, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(segment.m_start, segment.direction()), shape, 0.0f, segment.length(), hitInfo);
}

// --- Misc Helper ---

inline void projectPointSetOntoAxis(std::span<const Vector2D> pts, const Vector2D& axis,
                                    real_t& outMin, real_t& outMax)
{
    outMin = std::numeric_limits<real_t>::max();
    outMax = std::numeric_limits<real_t>::lowest();

    for (const auto& p : pts)
    {
        real_t v = p.dot(axis);
        if (v < outMin) outMin = v;
        if (v > outMax) outMax = v;
    }
}

inline void projectBBoxOntoAxis(const BBox2D& bbox, const Vector2D& axis,
                                real_t& outMin, real_t& outMax)
{
    Vector2D center  = (bbox.m_min + bbox.m_max) * 0.5f;
    Vector2D extents = (bbox.m_max - bbox.m_min) * 0.5f;

    // The "radius" of the box projection
    real_t radius = std::abs(extents.x * axis.x) + std::abs(extents.y * axis.y);
    real_t centerProj = center.dot(axis);

    outMin = centerProj - radius;
    outMax = centerProj + radius;
    //return { centerProj - radius, centerProj + radius }
}

// --- Bounding Box ---

bool intersect(const BBox2D& b1, const BBox2D& b2)
{
    //Separating Axis Theorem (SAT)
#if 0
    bool separated_x = (b1.m_min.x > b2.m_max.x) || (b2.m_min.x > b1.m_max.x);
    bool separated_y = (b1.m_min.y > b2.m_max.y) || (b2.m_min.y > b1.m_max.y);

    if (separated_x || separated_y)
        return false;

    return true;
#else
    return intervalsOverlap(b1.m_min.x, b1.m_max.x, b2.m_min.x, b2.m_max.x) &&
           intervalsOverlap(b1.m_min.y, b1.m_max.y, b2.m_min.y, b2.m_max.y);
#endif
}

// Generic SAT-based intersection for a convex polygon with a BBox
inline bool intersectBBoxWithConvexShape(const BBox2D& bbox, std::span<const Vector2D> verts)
{
    if (verts.size() < 3) return false; // degenerate

    real_t vMinX = verts[0].x, vMaxX = verts[0].x;
    real_t vMinY = verts[0].y, vMaxY = verts[0].y;

    for (size_t i = 1; i < verts.size(); ++i)
    {
        if (verts[i].x < vMinX) vMinX = verts[i].x;
        if (verts[i].x > vMaxX) vMaxX = verts[i].x;
        if (verts[i].y < vMinY) vMinY = verts[i].y;
        if (verts[i].y > vMaxY) vMaxY = verts[i].y;
    }

    if (vMaxX < bbox.m_min.x || vMinX > bbox.m_max.x) return false;
    if (vMaxY < bbox.m_min.y || vMinY > bbox.m_max.y) return false;

    real_t bmin, bmax, tmin, tmax;

    for (size_t i = 0; i < verts.size(); ++i)
    {
        // Get the edge between current and next vertex (wrapping around)
        const Vector2D& p1 = verts[i];
        const Vector2D& p2 = verts[(i + 1) % verts.size()];
        
        Vector2D edge = p2 - p1;
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges (nearly identical vertices)
        if (axis.isZero()) continue;

        projectBBoxOntoAxis    (bbox , axis, bmin, bmax);
        projectPointSetOntoAxis(verts, axis, tmin, tmax);

        if (!intervalsOverlap(bmin, bmax, tmin, tmax))
            return false;
    }

    return true; // No separating axis found
}

bool intersectBBoxWithPolygon(const BBox2D& bbox, const Polygon2D& polygon)
{
    if (polygon.isConvex())
    {
        const ConvexPolygon2D* convex = static_cast<const ConvexPolygon2D*>(&polygon);
        return intersectBBoxWithConvexShape(bbox, polygon.vertices());
    }
    
    Vector2D box[4] =
    {
        bbox.m_min,
        Vector2D(bbox.m_max.x, bbox.m_min.y),
        bbox.m_max,
        Vector2D(bbox.m_min.x, bbox.m_max.y)
    };

    for (const Vector2D& p : box)
        if (polygon.contains(p))
            return true;

    for (const Vector2D& p : polygon.m_vertices)
        if (bbox.contains(p))
            return true;

    Vector2D boxEdges[4][2] =
    {
        { box[0], box[1] },
        { box[1], box[2] },
        { box[2], box[3] },
        { box[3], box[0] }
    };

    const size_t n = polygon.m_vertices.size();

    for (size_t i = 0; i < n; ++i)
    {
        Vector2D a = polygon.m_vertices[i];
        Vector2D b = polygon.m_vertices[(i + 1) % n];

        for (int j = 0; j < 4; ++j)
            if (intersectSegmentWithSegment(a, b, boxEdges[j][0], boxEdges[j][1]))
                return true;
    }

    return false;
}

bool intersect(const BBox2D& bbox, const Circle2D& circle)
{
    real_t distSq = bbox.minDistanceSquared(circle.m_center);
    real_t rSq = circle.m_radius * circle.m_radius;
    return distSq <= rSq;
}

// --- Polygons --- 

bool intersectConvexPolygonWithConvexPolygon(const std::span<const Vector2D> v1, const std::span<const Vector2D> v2)
{
    const int n1 = (int)v1.size();
    const int n2 = (int)v2.size();

    if (n1 < 3 || n2 < 3)
        return false;

    real_t min1, max1, min2, max2;

    for (int i = 0; i < n1; ++i)
    {
        const Vector2D& a = v1[i];
        const Vector2D& b = v1[(i + 1) % n1];
        Vector2D edge = b - a;

        // Normal
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        projectPointSetOntoAxis(v1, axis, min1, max1);
        projectPointSetOntoAxis(v2, axis, min2, max2);

        if (!intervalsOverlap(min1, max1, min2, max2))
            return false; // separating axis found
    }

    for (int i = 0; i < n2; ++i)
    {
        const Vector2D& a = v2[i];
        const Vector2D& b = v2[(i + 1) % n2];
        Vector2D edge = b - a;

        // Normal
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        projectPointSetOntoAxis(v1, axis, min1, max1);
        projectPointSetOntoAxis(v2, axis, min2, max2);

        if (!intervalsOverlap(min1, max1, min2, max2))
            return false;
    }

    // No separating axis -> intersection
    return true;
}

bool intersectPolygonWithPolygon(const std::span<const Vector2D> v1, const std::span<const Vector2D> v2)
{
    const int n1 = (int)v1.size();
    const int n2 = (int)v2.size();

    if (n1 < 3 || n2 < 3)
        return false;

    for (const auto& v : v1)
    {
        if (isPointInsidePolygon(v, v2))
            return true;
    }

    for (const auto& v : v2)
    {
        if (isPointInsidePolygon(v, v1))
            return true;
    }

    // --- 2. Edge intersection tests ---
    for (size_t i = 0; i < n1; ++i)
    {
        Vector2D a1 = v1[i];
        Vector2D a2 = v1[(i + 1) % n1];

        for (size_t j = 0; j < n2; ++j)
        {
            Vector2D b1 = v2[j];
            Vector2D b2 = v2[(j + 1) % n2];

            if (intersectSegmentWithSegment(a1, a2, b1, b2))
                return true;
        }
    }

   return false;
}

bool intersectPolygonWithCircle(const std::span<const Vector2D> poly, const Circle2D& circle)
{
    if (isPointInsidePolygon(circle.m_center, poly))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    real_t r = circle.m_radius;

    for (int i = 0; i < poly.size(); ++i)
    {
        const Vector2D& a = poly[i];
        const Vector2D& b = poly[(i + 1) % poly.size()];

        if (distancePointToSegment(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

bool intersectPolygonWithCircle(const ConvexPolygon2D& polygon, const Circle2D& circle)
{
    if (polygon.contains(circle.centroid()))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    real_t r = circle.m_radius;

    for (int i = 0; i < polygon.vertexCount(); ++i)
    {
        const Vector2D& a = polygon.wrappedVertexAt(i);
        const Vector2D& b = polygon.wrappedVertexAt(i+1);

        if (distancePointToSegment(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

// --- Polygonal ---

bool intersect(const IPolygonalShape2D& polygon1, const IPolygonalShape2D& polygon2)
{
    if (isPolygonalShape(polygon1.type()) && isPolygonalShape(polygon2.type()))
        return intersectConvexPolygonWithConvexPolygon(polygon1.vertices(), polygon2.vertices());

    return intersectPolygonWithPolygon(polygon1.vertices(), polygon2.vertices());
}

bool intersect(const IPolygonalShape2D& polygon, const Circle2D& circle)
{
    if (polygon.contains(circle.centroid()))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    size_t n = polygon.vertexCount();
    real_t r = circle.m_radius;

    for (int i = 0; i < polygon.vertexCount(); ++i)
    {
        const Vector2D& a = polygon[i];
        const Vector2D& b = polygon[(i + 1) % n];

        if (distancePointToSegment(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

// --- Circles ---

bool intersectCircleWithCircle(const Circle2D& circle1, const Circle2D& circle2)
{
    return circle1.intersects(circle2);
}

// --- Generic ---

bool intersect(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2)
{
    ShapeType2D type1 = shape1.type();
    ShapeType2D type2 = shape2.type();

    if (type1 > type2)
        return intersect(shape2, shape1);

    if (isConvexPolygonal(type1))
    {
        const IPolygonalShape2D* s1 = shape1.polygonal();
        if (isConvexPolygonal(type2))
        {
            const IPolygonalShape2D* s2 = shape2.polygonal();
            return intersectConvexPolygonWithConvexPolygon(s1->vertices(), s2->vertices());
        }
        if (type2 == ShapeType2D::SHAPE2D_CIRCLE)
        {
            const Circle2D* s2 = shape2.shape_cast<Circle2D>();
            return intersectPolygonWithCircle(s1->vertices(), *s2);
        }
    }
    if (isPolygonalShape(type1))
    {
        const IPolygonalShape2D* s1 = shape1.polygonal();
        if (isPolygonalShape(type2))
        {
            const IPolygonalShape2D* s2 = shape2.polygonal();
            return intersectPolygonWithPolygon(s1->vertices(), s2->vertices());
        }
        if (type2 == ShapeType2D::SHAPE2D_CIRCLE)
        {
            const Circle2D* s2 = shape2.shape_cast<Circle2D>();
            return intersectPolygonWithCircle(shape1.polygonal()->vertices(), *shape2.shape_cast<Circle2D>());
        }
    }
    else if (type1 == ShapeType2D::SHAPE2D_CIRCLE)
    {
        const Circle2D* s1 = shape1.shape_cast<Circle2D>();
        if (type2 == ShapeType2D::SHAPE2D_CIRCLE)
        {
            const Circle2D* s2 = shape2.shape_cast<Circle2D>();
            return intersectCircleWithCircle(*s1, *s2);
        }
    }
    return false;
}

} // namespace Math

} // namespace Arns