/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

namespace Arns::geomath
{
    

namespace detail
{

    real_t distancePointToSegment(const Vector2D &point, const Vector2D &segmentStart, const Vector2D &segmentEnd,
                                  Vector2D *closestPoint = nullptr);

    real_t distanceSegmentToSegment(const Vector2D &s1, const Vector2D &s2, const Vector2D &k1, const Vector2D &k2,
                                    Vector2D *closestPoint1 = nullptr, Vector2D *closestPoint2 = nullptr);

    real_t distancePointToPolygonEdges(const Vector2D &point, std::span<const Vector2D> vertices,
                                       Vector2D *closestPoint = nullptr);

    real_t distanceSegmentToEdgeList(const Segment2D &segment, std::span<const Vector2D> vertices,
                                     bool closed, Vector2D *closestSegmentPoint = nullptr, Vector2D *closestEdgePoint = nullptr);

    real_t distanceSegmentToPolygon(const Segment2D &segment, std::span<const Vector2D> vertices,
                                    bool convex, Vector2D *closestSegmentPoint = nullptr, Vector2D *closestPolygonPoint = nullptr);

    real_t distanceBBoxToPolygon(const BBox2D &bbox, std::span<const Vector2D> vertices, bool convex,
                                 Vector2D *closestPoint1, Vector2D *closestPoint2);

    real_t distancePolygonToPolygon(std::span<const Vector2D> vertices1, bool convex1,
                                    std::span<const Vector2D> vertices2, bool convex2, Vector2D *closestPoint1, Vector2D *closestPoint2);

    real_t distancePolygonToCircle(std::span<const Vector2D> vertices, bool convex,
                                   const Circle2D &circle, Vector2D *closestPoint1, Vector2D *closestPoint2);

    // ========================================================================
    // 1. Point vs Shape Queries (1 Output Point)
    // ========================================================================

    real_t distanceImpl(const Vector2D &point, const Segment2D &segment, Vector2D *c1 = nullptr, Vector2D *c2 = nullptr);
    real_t distanceImpl(const Vector2D &point, const BBox2D &bbox, Vector2D *c1 = nullptr, Vector2D *c2 = nullptr);

    template <PolygonalShape Polygon>
    real_t distanceImpl(const Vector2D &point, const Polygon &polygon, Vector2D *c1 = nullptr, Vector2D *c2 = nullptr)
    {
        const auto vertices = polygon.vertices();
        const bool inside = polygon.isConvex()
                                ? detail::isPointInsideConvexPolygon(point, vertices)
                                : detail::isPointInsidePolygon(point, vertices);

        if (c1) *c1 = point;

        if (inside)
        {
            if (c2) *c2 = point;
            return real_t{0};
        }

        return detail::distancePointToPolygonEdges(point, vertices, c2);
    }

    real_t distanceImpl(const Vector2D &point, const Circle2D &circle, Vector2D *c1 = nullptr, Vector2D *c2 = nullptr);

    // ========================================================================
    // 2. Line Segment vs Shape Queries (2 Output Points)
    // ========================================================================

    real_t distanceImpl(const Segment2D &segment1, const Segment2D &segment2, Vector2D *closestPoint1 = nullptr, Vector2D *closestPoint2 = nullptr);
    real_t distanceImpl(const Segment2D &segment, const BBox2D &bbox, Vector2D *closestPoint1 = nullptr, Vector2D *closestPoint2 = nullptr);

    template <PolygonalShape Polygon>
    real_t distanceImpl(const Segment2D &segment, const Polygon &polygon,
                        Vector2D *closestPoint1 = nullptr, Vector2D *closestPoint2 = nullptr)
    {
        return detail::distanceSegmentToPolygon(
            segment, polygon.vertices(), polygon.isConvex(), closestPoint1, closestPoint2);
    }

    real_t distanceImpl(const Segment2D &segment, const Circle2D &circle, Vector2D *closestPoint1 = nullptr, Vector2D *closestPoint2 = nullptr);

    // ========================================================================
    // 3. Shape vs Shape Queries (2 Output Points)
    // ========================================================================

    real_t distanceImpl(const BBox2D &bbox1, const BBox2D &bbox2, Vector2D *closestPoint1, Vector2D *closestPoint2);
    real_t distanceImpl(const BBox2D &bbox, const Circle2D &circle, Vector2D *closestPoint1, Vector2D *closestPoint2);
    real_t distanceImpl(const Circle2D &circle1, const Circle2D &circle2, Vector2D *closestPoint1, Vector2D *closestPoint2);

    template <PolygonalShape Polygon>
    real_t distanceImpl(const BBox2D &bbox, const Polygon &polygon, Vector2D *closestPoint1, Vector2D *closestPoint2)
    {
        return distanceBBoxToPolygon(bbox, polygon.vertices(), polygon.isConvex(), closestPoint1, closestPoint2);
    }

    template <PolygonalShape Polygon1, PolygonalShape Polygon2>
    real_t distanceImpl(const Polygon1 &polygon1, const Polygon2 &polygon2,
                        Vector2D *closestPoint1, Vector2D *closestPoint2)
    {
        return distancePolygonToPolygon(polygon1.vertices(), polygon1.isConvex(),
                                        polygon2.vertices(), polygon2.isConvex(), closestPoint1, closestPoint2);
    }

    template <PolygonalShape Polygon>
    real_t distanceImpl(const Polygon &polygon, const Circle2D &circle,
                        Vector2D *closestPoint1, Vector2D *closestPoint2)
    {
        return distancePolygonToCircle(polygon.vertices(), polygon.isConvex(), circle, closestPoint1, closestPoint2);
    }

    template <typename A, typename Object>
    concept HasDistanceOrder = requires(const A &a, const Object &b, Vector2D *pointA, Vector2D *pointB) {
        { distanceImpl(a, b, pointA, pointB) } -> std::convertible_to<real_t>;
    };

    template <typename A, typename Object>
    concept HasDistance = HasDistanceOrder<A, Object> || HasDistanceOrder<Object, A>;

} // namespace detail

//==================================
// Public Distance dispatch function 
//==================================

template <typename A, typename Object>
    requires (requires(const A& a, const Object& b, Vector2D* p1, Vector2D* p2) {
                  detail::distanceImpl(a, b, p1, p2);
              } || requires(const A& a, const Object& b, Vector2D* p1, Vector2D* p2) {
                  detail::distanceImpl(b, a, p2, p1);
              })
real_t distance(const A& a, const Object& b, Vector2D* p1 = nullptr, Vector2D* p2 = nullptr)
{
    if constexpr (requires { detail::distanceImpl(a, b, p1, p2); })
    {
        return detail::distanceImpl(a, b, p1, p2);
    }
    else
    {
        return detail::distanceImpl(b, a, p2, p1);
    }
}
    
} // namespace Arns::geomath