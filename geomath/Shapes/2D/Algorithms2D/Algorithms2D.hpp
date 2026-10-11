/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "../Line2D.hpp"
#include "../Ray2D.hpp"
#include "../Variants/ShapeConcepts.hpp"

#include <math.h>
#include <vector>
#include <span>

namespace Arns
{

namespace geomath
{

// Forward declarations
class Ray2D;
class Segment2D;
class BBox2D;
class IBaseShape2D;
class IFiniteShape2D;
class IPolygonalShape2D;
class Circle2D;
class Triangle2D;
class Rectangle2D;
class Polygon2D;
class ConvexPolygon2D;


struct HitInfo2D
{
    real_t t;
    Vector2D intersectionPoint;
    Vector2D normal;
    //IShape2D* shape;
};

struct SweepResult2D
{
    bool hit = false;
    real_t toi = 1.0;
    Vector2D normal{};
    Vector2D point{};
};


//

namespace detail
{

bool isPointInsideConvexPolygon(const Vector2D& point, const std::span<const Vector2D> vertices);

bool isPointInsidePolygon(const Vector2D& point, const std::span<const Vector2D> vertices);

bool isPolygonInsideConvexPolygon(const std::span<const Vector2D> poly, const std::span<const Vector2D> convexPoly);

}

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd);

bool isPointOnSegment(const Vector2D& point, const Segment2D& segment);

bool isSegmentOnSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2);

bool isSegmentOnSegment(const Segment2D& segment1, const Segment2D& segment2);

// --- Point-in-Shape containment

bool isPointInsideBBox(const Vector2D& point, const BBox2D& bbox);

bool isPointInsideTriangle(const Vector2D& point, const Triangle2D& triangle);

bool isPointInsideRectangle(const Vector2D& point, const Rectangle2D& rectangle);

bool isPointInsideConvexPolygon(const Vector2D& point, const ConvexPolygon2D& polygon);

bool isPointInsidePolygon(const Vector2D& point, const Polygon2D& polygon);

//bool isPointInsidePolygon(const Vector2D& point, const IPolygonalShape2D& polygon);

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle);

//bool isPointInShape(const Vector2D& point, const IBaseShape2D& shape);

// --- Segment-in-Shape containment

bool isSegmentInsideBBox(const Segment2D &segment, const BBox2D& bbox);

bool isSegmentInsideTriangle(const Segment2D &segment, const Triangle2D& triangle);

bool isSegmentInsideRectangle(const Segment2D &segment, const Rectangle2D& rectangle);

bool isSegmentInsideConvexPolygon(const Segment2D &segment, const ConvexPolygon2D& polygon);

bool isSegmentInsidePolygon(const Segment2D &segment, const Polygon2D& polygon);

bool isSegmentInsideCircle(const Segment2D &segment, const Circle2D& circle);

//bool isSegmentInsideShape(const Segment2D& segment, const IBaseShape2D& shape);

//

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2);

//bool isPolygonalInsideBBox(const IPolygonalShape2D& polygon, const BBox2D& bbox);

bool isCircleInsideBBox(const Circle2D& circle, const BBox2D& bbox);

//

bool isCircleInsideCircle(const Circle2D& circle1, const Circle2D& circle2);



// Distance calculation algorithms

namespace detail
{

real_t distancePointToSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd, 
    Vector2D* closestPoint = nullptr);

real_t distanceSegmentToSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, 
    Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

real_t distancePointToPolygonEdges(const Vector2D& point, std::span<const Vector2D> vertices,
    Vector2D* closestPoint = nullptr);

real_t distanceSegmentToEdgeList(const Segment2D& segment, std::span<const Vector2D> vertices,
    bool closed, Vector2D* closestSegmentPoint = nullptr, Vector2D* closestEdgePoint = nullptr);

real_t distanceSegmentToPolygon(const Segment2D& segment, std::span<const Vector2D> vertices,
    bool convex, Vector2D* closestSegmentPoint = nullptr, Vector2D* closestPolygonPoint = nullptr);

real_t distanceBBoxToPolygon(const BBox2D& bbox, std::span<const Vector2D> vertices, bool convex,
    Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distancePolygonToPolygon(std::span<const Vector2D> vertices1, bool convex1,
    std::span<const Vector2D> vertices2, bool convex2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distancePolygonToCircle(std::span<const Vector2D> vertices, bool convex,
    const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);

// ========================================================================
// 1. Point vs Shape Queries (1 Output Point)
// ========================================================================

real_t distanceImpl(const Vector2D& point, const Segment2D& segment, Vector2D* c1 = nullptr, Vector2D* c2 = nullptr);
real_t distanceImpl(const Vector2D& point, const BBox2D& bbox, Vector2D* c1 = nullptr, Vector2D* c2 = nullptr);

template <PolygonalShape Polygon>
real_t distanceImpl(const Vector2D& point, const Polygon& polygon, Vector2D* c1 = nullptr, Vector2D* c2 = nullptr)
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

real_t distanceImpl(const Vector2D& point, const Circle2D& circle, Vector2D* c1 = nullptr, Vector2D* c2 = nullptr);


// ========================================================================
// 2. Line Segment vs Shape Queries (2 Output Points)
// ========================================================================

real_t distanceImpl(const Segment2D& segment1, const Segment2D& segment2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);
real_t distanceImpl(const Segment2D& segment, const BBox2D& bbox, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

template <PolygonalShape Polygon>
real_t distanceImpl(const Segment2D& segment, const Polygon& polygon,
    Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr)
{
    return detail::distanceSegmentToPolygon(
        segment, polygon.vertices(), polygon.isConvex(), closestPoint1, closestPoint2);
}

real_t distanceImpl(const Segment2D& segment, const Circle2D& circle, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);


// ========================================================================
// 3. Shape vs Shape Queries (2 Output Points)
// ========================================================================

real_t distanceImpl(const BBox2D& bbox1, const BBox2D& bbox2, Vector2D* closestPoint1, Vector2D* closestPoint2);
real_t distanceImpl(const BBox2D& bbox, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);
real_t distanceImpl(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

template <PolygonalShape Polygon>
real_t distanceImpl(const BBox2D& bbox, const Polygon& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distanceBBoxToPolygon(bbox, polygon.vertices(), polygon.isConvex(), closestPoint1, closestPoint2);
}

template <PolygonalShape Polygon1, PolygonalShape Polygon2>
real_t distanceImpl(const Polygon1& polygon1, const Polygon2& polygon2,
    Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distancePolygonToPolygon(polygon1.vertices(), polygon1.isConvex(),
        polygon2.vertices(), polygon2.isConvex(), closestPoint1, closestPoint2);
}

template <PolygonalShape Polygon>
real_t distanceImpl(const Polygon& polygon, const Circle2D& circle,
    Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distancePolygonToCircle(polygon.vertices(), polygon.isConvex(), circle, closestPoint1, closestPoint2);
}

} // namespace detail

//  --- Generic Distance function ---

template <typename A, typename B>
    requires (requires(const A& a, const B& b, Vector2D* p1, Vector2D* p2) {
                  detail::distanceImpl(a, b, p1, p2);
              } || requires(const A& a, const B& b, Vector2D* p1, Vector2D* p2) {
                  detail::distanceImpl(b, a, p2, p1);
              })
real_t distance(const A& a, const B& b, Vector2D* p1 = nullptr, Vector2D* p2 = nullptr)
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



// Intersection calculation algorithms



namespace detail
{

    bool intersectRayWithSegment(const Ray2D &ray, const Vector2D &p1, const Vector2D &p2,
                                 real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    bool rayBoundsInterval(const Ray2D &ray, const BBox2D &bbox,
                           real_t t_min, real_t t_max, real_t &t_enter, real_t &t_exit);

    bool intersectRayWithPolygon(const Ray2D &ray, std::span<const Vector2D> vertices,
                                 real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    bool intersectSegmentWithSegment(const Vector2D &p1, const Vector2D &p2, const Vector2D &q1, const Vector2D &q2, HitInfo2D *hitInfo = nullptr);

    bool intersectSegmentWithSegmentStrict(const Vector2D &p1, const Vector2D &p2, const Vector2D &q1, const Vector2D &q2, HitInfo2D *hitInfo = nullptr);

    //    --- Rays ---

    bool intersectRayImpl(const Ray2D &ray, const Segment2D &segment, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    bool intersectRayImpl(const Ray2D &ray, const BBox2D &bbox, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    template <PolygonalShape Polygon>
    bool intersectRayImpl(const Ray2D &ray, const Polygon &polygon, real_t t_min = 0,
                          real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr)
    {
        return detail::intersectRayWithPolygon(ray, polygon.vertices(), t_min, t_max, hitInfo);
    }

    bool intersectRayImpl(const Ray2D &ray, const Circle2D &circle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    template <typename Object>
    concept HasRayIntersectionImpl =
        requires(
            const Ray2D &ray,
            const Object &object,
            real_t tMin,
            real_t tMax,
            HitInfo2D *hitInfo) {
            {
                intersectRayImpl(
                    ray, object, tMin, tMax, hitInfo)
            } -> std::convertible_to<bool>;
        };

    bool intersectImpl(const BBox2D &bbox1, const BBox2D &bbox2);
    bool intersectImpl(const BBox2D &bbox, const Circle2D &circle);
    bool intersectImpl(const Circle2D &circle1, const Circle2D &circle2);

    bool intersectBBoxWithPolygon(const BBox2D &bbox, std::span<const Vector2D> vertices, bool convex);
    bool intersectPolygonWithPolygon(std::span<const Vector2D> vertices1, bool convex1,
                                     std::span<const Vector2D> vertices2, bool convex2);
    bool intersectPolygonWithCircle(std::span<const Vector2D> vertices, bool convex, const Circle2D &circle);

    template <PolygonalShape Polygon>
    bool intersectImpl(const BBox2D &bbox, const Polygon &polygon)
    {
        return intersectBBoxWithPolygon(bbox, polygon.vertices(), polygon.isConvex());
    }

    template <PolygonalShape Polygon1, PolygonalShape Polygon2>
    bool intersectImpl(const Polygon1 &polygon1, const Polygon2 &polygon2)
    {
        return intersectPolygonWithPolygon(polygon1.vertices(), polygon1.isConvex(),
                                           polygon2.vertices(), polygon2.isConvex());
    }

    template <PolygonalShape Polygon>
    bool intersectImpl(const Polygon &polygon, const Circle2D &circle)
    {
        return intersectPolygonWithCircle(polygon.vertices(), polygon.isConvex(), circle);
    }
}

//    --- Ray intersection function ---

template <typename B> requires (detail::HasRayIntersectionImpl<B>)
bool intersect(
    const Ray2D& first, const B& second, real_t tMin = real_t{0},
        real_t tMax = std::numeric_limits<real_t>::max(),
        HitInfo2D* hitInfo = nullptr)
    {
        return detail::intersectRayImpl(first, second, tMin, tMax, hitInfo);
    }

//    --- Segment intersection function ---

template <typename Object> requires detail::HasRayIntersectionImpl<Object>
bool intersect(const Segment2D& segment, const Object& object, HitInfo2D* hitInfo = nullptr)
{
    const Ray2D ray(segment.m_start, segment.deltaVector());
    return intersect(ray, object, real_t{0}, real_t{1}, hitInfo);
}

//    --- Shape-Shape intersection function ---

template <typename A, typename B>
    requires (requires(const A& a, const B& b) { detail::intersectImpl(a, b); } ||
              requires(const A& a, const B& b) { detail::intersectImpl(b, a); })
bool intersect(const A& a, const B& b)
{
    if constexpr (requires { detail::intersectImpl(a, b); })
        return detail::intersectImpl(a, b);
    else
        return detail::intersectImpl(b, a);
}

//================================
// CCD (sweep)
//================================

/*
SweepResult sweepConservativeAdvancement(
    const BBox2D& a, Vector2D velocityA,
    const BBox2D& b, Vector2D velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15,
    real_t epsilon = real_t(1e-4));

SweepResult sweepConservativeAdvancement(
    const BBox2D& a, Vector2D velocityA,
    const IFiniteShape2D& b, Vector2D velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15,
    real_t epsilon = real_t(1e-4));

SweepResult sweepConservativeAdvancement(
    const IFiniteShape2D& a, Vector2D velocityA,
    const BBox2D& b, Vector2D velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15,
    real_t epsilon = real_t(1e-4));
*/


// ========================================================================
// 1. Point vs Shape Queries
// ========================================================================

/*
SweepResult2D sweep(const Vector2D& point, const Vector2D& velocityA, const Segment2D& segment, const Vector2D& velocityB, 
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const Vector2D& point, const Vector2D& velocityA, const BBox2D& bbox,       const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const Vector2D& point, const Vector2D& velocityA, const IPolygonalShape2D& polygon, const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const Vector2D& point, const Vector2D& velocityA, const Circle2D& circle, const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
*/

// ========================================================================
// 2. Line Segment vs Shape Queries
// ========================================================================

/*
SweepResult2D sweep(const Segment2D& segment1, const Vector2D& velocityA, const Segment2D& segment2,        const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const Segment2D& segment,  const Vector2D& velocityA, const BBox2D& bbox,               const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const Segment2D& segment,  const Vector2D& velocityA, const IPolygonalShape2D& polygon, const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const Segment2D& segment,  const Vector2D& velocityA, const Circle2D& circle,           const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
*/

// ========================================================================
// 3. Shape vs Shape Queries
// ========================================================================

/*
SweepResult2D sweep(const BBox2D& bbox1, const Vector2D& velocityA, const BBox2D& bbox2,           const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const BBox2D& bbox,  const Vector2D& velocityA, const IPolygonalShape2D& poly, const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const BBox2D& bbox,  const Vector2D& velocityA, const Circle2D& circle,        const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));

SweepResult2D sweep(const IPolygonalShape2D& poly1, const Vector2D& velocityA, const IPolygonalShape2D& poly2, const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
SweepResult2D sweep(const IPolygonalShape2D& poly1, const Vector2D& velocityA, const Circle2D& circle2, const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));

SweepResult2D sweep(const Circle2D& circle1, const Vector2D& velocityA, const Circle2D& circle2,      const Vector2D& velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
*/
/*
SweepResult2D sweep(const IFiniteShape2D& a, Vector2D velocityA, const IFiniteShape2D& b, Vector2D velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4));
*/

} // namespace Math

} // namespace Arns