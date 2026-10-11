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
#include "../Variants/ShapeConcepts.hpp"
#include "ShapeVariantAlgorithms.hpp"

#include <algorithm>
#include <cassert>
#include <span>
#include <vector>

namespace Arns::geomath
{

namespace detail
{

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

bool isPolygonInsideConvexPolygon(const std::span<const Vector2D> poly, const std::span<const Vector2D> convexPoly)
{
    for (const auto& point : poly)
    {
        if (!detail::isPointInsideConvexPolygon(point, convexPoly))
            return false;
    }
    return true;
}

bool isPolygonInsidePolygon(const std::span<const Vector2D> poly1, const std::span<const Vector2D> poly2)
{
    for (const auto& point : poly1)
    {
        if (!detail::isPointInsidePolygon(point, poly2))
        {
            return false;
        }
    }

    if (detail::edgeIntersectionPolygonWithPolygon(poly1, poly2))
        return false;

    return true;
}

}


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

bool isPointInsideConvexPolygon(const Vector2D& point, const ConvexPolygon2D& polygon)
{
    return detail::isPointInsideConvexPolygon(point, polygon.vertices());
}

bool isPointInsidePolygon(const Vector2D& point, const Polygon2D& polygon)
{
    return detail::isPointInsidePolygon(point, polygon.vertices());
}

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle)
{
    return (point - circle.m_center).lengthSquared() <= circle.m_radius * circle.m_radius;
}

// Segment inside shape

bool isSegmentInsideBBox(const Segment2D &segment, const BBox2D& bbox)
{
    return isPointInsideBBox(segment.m_start, bbox) &&
           isPointInsideBBox(segment.m_end, bbox);
}

bool isSegmentInsideTriangle(const Segment2D &segment, const Triangle2D& triangle)
{
    return isPointInsideTriangle(segment.m_start, triangle) &&
           isPointInsideTriangle(segment.m_end, triangle);
}

bool isSegmentInsideRectangle(const Segment2D &segment, const Rectangle2D& rectangle)
{
    return isPointInsideRectangle(segment.m_start, rectangle) &&
           isPointInsideRectangle(segment.m_end, rectangle);
}

bool isSegmentInsideConvexPolygon(const Segment2D &segment, const ConvexPolygon2D& polygon)
{
    return isPointInsideConvexPolygon(segment.m_start, polygon) &&
           isPointInsideConvexPolygon(segment.m_end, polygon);
}

bool isSegmentInsidePolygon(const Segment2D &segment, const Polygon2D& polygon)
{
    if (!isPointInsidePolygon(segment.m_start, polygon) ||
        !isPointInsidePolygon(segment.m_end, polygon))
    {
        return false;
    }

    if (!polygon.isConvex())
    {
        return !intersect(segment, polygon);
    }
    
    return true;
}

bool isSegmentInsideCircle(const Segment2D &segment, const Circle2D& circle)
{
    return isPointInsideCircle(segment.m_start, circle) &&
           isPointInsideCircle(segment.m_end, circle);
}

//Shape inside BBox

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2)
{
    return bbox2.contains(bbox1);
}

bool isCircleInsideBBox(const Circle2D& circle, const BBox2D& bbox)
{
    return bbox.contains(circle.m_center) && bbox.minDistanceSquared(circle.m_center) < circle.m_radius * circle.m_radius;
}

//

//

bool isBBoxInsideConvexPolygon(const BBox2D& bbox, const std::span<const Vector2D> convexPoly)
{
    return detail::isPointInsideConvexPolygon(bbox.bottomLeft(), convexPoly) && detail::isPointInsideConvexPolygon(bbox.topLeft(), convexPoly) && 
        detail::isPointInsideConvexPolygon(bbox.topRight(), convexPoly) && detail::isPointInsideConvexPolygon(bbox.bottomRight(), convexPoly);
}

bool isBBoxInsidePolygon(const BBox2D& bbox, const std::span<const Vector2D> poly)
{
    std::vector<Vector2D> bboxPoly = {bbox.bottomLeft(), bbox.bottomRight(), bbox.topRight(), bbox.topLeft()};

    if (!detail::isPointInsidePolygon(bbox.bottomLeft(), poly) || !detail::isPointInsidePolygon(bbox.bottomRight(), poly) ||
        !detail::isPointInsidePolygon(bbox.topRight(), poly)   || !detail::isPointInsidePolygon(bbox.topLeft(), poly))
            return false;

    if (detail::edgeIntersectionPolygonWithPolygon(bboxPoly, poly))
        return false;

    return true;
}



bool isCircleInsidePolygon(const Circle2D& circle, const std::span<const Vector2D> poly);

bool isCircleInsideConvexPolygon(const Circle2D& circle, const std::span<const Vector2D> convexPoly)
{
    if (!detail::isPointInsideConvexPolygon(circle.m_center, convexPoly))
    {
        return false;
    }

    return isCircleInsidePolygon(circle, convexPoly);
}

bool isCircleInsidePolygon(const Circle2D& circle, const std::span<const Vector2D> poly)
{
    if (!detail::isPointInsidePolygon(circle.m_center, poly))
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

namespace detail
{

    template <PolygonalShape Shape>
    bool isShapeInsideCircle(const Shape& shape, const Circle2D& outerCircle)
    {
        return isPolygonInsideCircle(shape.vertices(), outerCircle);
    }

    bool isShapeInsideCircle(const Circle2D& inner, const Circle2D& outerCircle)
    {
        return isCircleInsideCircle(inner, outerCircle);
    }


    // -------------------------------------------------------------------------
    // Shape -> convex polygon
    // -------------------------------------------------------------------------

    template <PolygonalShape Shape>
    bool isShapeInsideConvexPolygon(
        const Shape& shape,
        std::span<const Vector2D> convexPolygon)
    {
        if (shape.isConvex())
        {
            return isPolygonInsideConvexPolygon(
                shape.vertices(),
                convexPolygon);
        }

        return isPolygonInsideConvexPolygon(
            shape.vertices(),
            convexPolygon);
    }

    inline bool isShapeInsideConvexPolygon(
        const Circle2D& circle,
        std::span<const Vector2D> convexPolygon)
    {
        return isCircleInsideConvexPolygon(
            circle,
            convexPolygon);
    }


    // -------------------------------------------------------------------------
    // Shape -> arbitrary polygon
    // -------------------------------------------------------------------------

    template <PolygonalShape Shape>
    bool isShapeInsidePolygon(
        const Shape& shape,
        std::span<const Vector2D> polygon)
    {
        return isPolygonInsidePolygon(
            shape.vertices(),
            polygon);
    }

    inline bool isShapeInsidePolygon(
        const Circle2D& circle,
        std::span<const Vector2D> polygon)
    {
        return isCircleInsidePolygon(
            circle,
            polygon);
    }


    // -------------------------------------------------------------------------
    // Shape -> Shape
    // -------------------------------------------------------------------------

    template <FiniteShape2D Shape>
    bool isShapeInsideShape(
        const Shape& shape,
        const PolygonalShape auto& shell)
    {
        if (shell.isConvex())
        {
            return isShapeInsideConvexPolygon(
                shape,
                shell.vertices());
        }

        return isShapeInsidePolygon(
            shape,
            shell.vertices());
    }

    template <FiniteShape2D Shape>
    bool isShapeInsideShape(
        const Shape& shape,
        const Circle2D& circle)
    {
        return isShapeInsideCircle(shape, circle);
    }

}

bool isShapeInsideShape(const FiniteShapeVariant2D& inner, const FiniteShapeVariant2D& shell)
{
#if 0
    // Todo: Don't rely on the function we want to remove in the future
    return std::visit([](const auto& innerShape, const auto& shellShape)
    {
        return isShapeInsideShape(
            static_cast<const IFiniteShape2D&>(innerShape),
            static_cast<const IFiniteShape2D&>(shellShape));
    }, inner, shell);
#else
    return std::visit(
        []<FiniteShape2D Inner, FiniteShape2D Shell>(
            const Inner& inner,
            const Shell& shell)
        {
            return detail::isShapeInsideShape(inner, shell);
        }, inner, shell);
#endif
}

/*
bool isPointInsideShape(const Vector2D& point, const FiniteShapeVariant2D& shell)
{
    
}

bool isSegmentInsideShape(const Segment2D& segment, const FiniteShapeVariant2D& shell)
{
    
}
*/

//Consider throwing an exception after unreachable code asserts

}