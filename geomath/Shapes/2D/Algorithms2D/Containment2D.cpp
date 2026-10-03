/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms2D.hpp"
#include "CommonMath.hpp"
#include "../Ray2D.hpp"
#include "../Line2D.hpp"
#include "../BBox2D.hpp"
#include "../Circle2D.hpp"
#include "../Triangle2D.hpp"
#include "../Rectangle2D.hpp"
#include "../Polygon2D.hpp"
#include "../ConvexPolygon2D.hpp"

#include <cassert>
#include <span>

namespace Arns::geomath
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

bool isPointOnSegment(const Vector2D& point, const Segment2D& segment)
{
    return isPointOnSegment(point, segment.m_start, segment.m_end);
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
    const real_t areaABC = triangle.area();
    const auto& a = triangle.a();
    const auto& b = triangle.b();
    const auto& c = triangle.c();
    const real_t areaPBC = real_t{0.5} * std::abs((b.x - point.x) * (c.y - point.y) - (c.x - point.x) * (b.y - point.y));
    const real_t areaPCA = real_t{0.5} * std::abs((c.x - point.x) * (a.y - point.y) - (a.x - point.x) * (c.y - point.y));
    const real_t areaPAB = real_t{0.5} * std::abs((a.x - point.x) * (b.y - point.y) - (b.x - point.x) * (a.y - point.y));
    return approximatelyEqual(areaABC, areaPBC + areaPCA + areaPAB);
}

bool isPointInsideRectangle(const Vector2D& point, const Rectangle2D& rectangle)
{
    const auto& a = rectangle.a();
    const auto& b = rectangle.b();
    const auto& c = rectangle.c();
    const auto& d = rectangle.d();

    return (b - a).cross(point - a) >= 0 &&
           (c - b).cross(point - b) >= 0 &&
           (d - c).cross(point - c) >= 0 &&
           (a - d).cross(point - d) >= 0;
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
    return isPointInsidePolygon(point, polygon.vertices());
}

bool isPointInsidePolygon(const Vector2D& point, const IPolygonalShape2D& polygon)
{
    return isPointInsidePolygon(point, polygon.vertices());
}

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle)
{
    return (point - circle.m_center).lengthSquared() <= circle.m_radius * circle.m_radius;
}

bool isPointInShape(const Vector2D& point, const IBaseShape2D& shape)
{
    return shape.contains(point);
}

// Segment inside shape

bool isSegmentInsideShape(const Segment2D& segment, const IBaseShape2D& shape)
{
    if (!shape.contains(segment.m_start) || !shape.contains(segment.m_end))
        return false;

    if (shape.type() == SHAPE2D_POLYGON)
    {
        return !intersect(segment, shape);
    }

    if (isConvexPolygonal(shape.type()) || shape.type() == SHAPE2D_CIRCLE)
        return true;

    assert(false && "Unreachable code reached!");
    return false;
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

bool isPolygonInsideConvexPolygon(const std::span<const Vector2D> poly, const std::span<const Vector2D> convexPoly)
{
    for (const auto& point : poly)
    {
        if (!isPointInsideConvexPolygon(point, convexPoly))
            return false;
    }
    return true;
}

bool isPolygonInsidePolygon(const std::span<const Vector2D> poly1, const std::span<const Vector2D> poly2)
{
    for (const auto& point : poly1)
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

bool isCircleInsidePolygon(const Circle2D& circle, const std::span<const Vector2D> poly);

bool isCircleInsideConvexPolygon(const Circle2D& circle, const std::span<const Vector2D> convexPoly)
{
    if (!isPointInsideConvexPolygon(circle.m_center, convexPoly))
    {
        return false;
    }

    return isCircleInsidePolygon(circle, convexPoly);
}

bool isCircleInsidePolygon(const Circle2D& circle, const std::span<const Vector2D> poly)
{
    if (!isPointInsidePolygon(circle.m_center, poly))
    {
        return false;
    }

    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i)
    {
        const Vector2D& v1 = poly[i];
        const Vector2D& v2 = poly[(i + 1) % n];

        real_t dist = detail::distancePointToSegment(circle.m_center, v1, v2);
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

bool isShapeInsidePolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> poly)
{
    if (isPolygonalShape(shape.type()))
    {
        return isPolygonInsidePolygon(shape.polygonal()->vertices(), poly);
    }
    else if (shape.type() == SHAPE2D_CIRCLE)
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
    if (circle1.m_radius > circle2.m_radius)
        return false;

    const real_t radiusDelta = circle2.m_radius - circle1.m_radius;
    return (circle2.m_center - circle1.m_center).lengthSquared() <= radiusDelta * radiusDelta;
}

bool isShapeInsideCircle(const IFiniteShape2D& shape, const Circle2D& circle)
{
    if (isPolygonalShape(shape.type()))
    {
        return isPolygonInsideCircle(shape.polygonal()->vertices(), circle);
    }
    else if (shape.type() == SHAPE2D_CIRCLE)
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

}