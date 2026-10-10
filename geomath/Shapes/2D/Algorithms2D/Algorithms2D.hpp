/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
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

/*
template <PolygonalShape Shape>
bool isSegmentInsideShape(const Segment2D &segment, const Shape &shape)
{
    if (shape.isConvex())
        return isPointInsideConvexPolygon(segment.m_start, shape.vertices()) &&
            return isPointInsideConvexPolygon(segment.m_end, shape.vertices());
    else
        return isPointInsidePolygon(segment.m_start, shape.vertices()) &&
            return isPointInsidePolygon(segment.m_end, shape.vertices());
}
*/

//bool isSegmentInsideShape(const Segment2D& segment, const IBaseShape2D& shape);

//

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2);

//bool isPolygonalInsideBBox(const IPolygonalShape2D& polygon, const BBox2D& bbox);

bool isCircleInsideBBox(const Circle2D& circle, const BBox2D& bbox);

//



bool isCircleInsideCircle(const Circle2D& circle1, const Circle2D& circle2);

//bool isShapeInsideConvexPolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> convexPoly);

//bool isShapeInsidePolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> poly);

//bool isShapeInsideCircle(const IFiniteShape2D& shape, const Circle2D& circle);

//bool isShapeInsideShape(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2);


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

}

//GJK dist?

// ========================================================================
// 1. Point vs Shape Queries (1 Output Point)
// ========================================================================

real_t distance(const Vector2D& point, const Segment2D& segment, Vector2D* closestPoint = nullptr);
real_t distance(const Vector2D& point, const BBox2D& bbox, Vector2D* closestPoint = nullptr);

template <PolygonalShape Polygon>
real_t distance(const Vector2D& point, const Polygon& polygon, Vector2D* closestPoint = nullptr)
{
    const auto vertices = polygon.vertices();
    const bool inside = polygon.isConvex()
        ? detail::isPointInsideConvexPolygon(point, vertices)
        : detail::isPointInsidePolygon(point, vertices);

    if (inside)
    {
        if (closestPoint)
            *closestPoint = point;
        return real_t{0};
    }

    return detail::distancePointToPolygonEdges(point, vertices, closestPoint);
}

real_t distance(const Vector2D& point, const Circle2D& circle, Vector2D* closestPoint = nullptr);

// ========================================================================
// 2. Line Segment vs Shape Queries (2 Output Points)
// ========================================================================

real_t distance(const Segment2D& segment1, const Segment2D& segment2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);
real_t distance(const Segment2D& segment, const BBox2D& bbox, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

template <PolygonalShape Polygon>
real_t distance(const Segment2D& segment, const Polygon& polygon,
    Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr)
{
    return detail::distanceSegmentToPolygon(
        segment, polygon.vertices(), polygon.isConvex(), closestPoint1, closestPoint2);
}

real_t distance(const Segment2D& segment, const Circle2D& circle, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

// ========================================================================
// 3. Shape vs Shape Queries (2 Output Points)
// ========================================================================

namespace detail
{

real_t distanceImpl(const BBox2D& bbox1, const BBox2D& bbox2, Vector2D* closestPoint1, Vector2D* closestPoint2);
real_t distanceImpl(const BBox2D& bbox, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);
real_t distanceImpl(const Circle2D& circle1, const Circle2D& circle2,
    Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceBBoxToPolygon(const BBox2D& bbox, std::span<const Vector2D> vertices, bool convex,
    Vector2D* closestPoint1, Vector2D* closestPoint2);
real_t distancePolygonToPolygon(std::span<const Vector2D> vertices1, bool convex1,
    std::span<const Vector2D> vertices2, bool convex2, Vector2D* closestPoint1, Vector2D* closestPoint2);
real_t distancePolygonToCircle(std::span<const Vector2D> vertices, bool convex,
    const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);

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

}

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

//    --- Rays ---

namespace detail
{

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, 
    real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool rayBoundsInterval(const Ray2D& ray, const BBox2D& bbox, 
    real_t t_min, real_t t_max, real_t& t_enter, real_t& t_exit);

}

bool intersect(const Ray2D& ray, const Segment2D& segment, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const BBox2D& bbox, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const IPolygonalShape2D& polygon, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const Circle2D& circle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const IBaseShape2D& shape, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

//    --- Segments ---

bool intersectSegmentWithSegment(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment1, const Segment2D& segment2, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const BBox2D& bbox, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const Circle2D& circle, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const Triangle2D& triangle, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const IPolygonalShape2D& polygon, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const IBaseShape2D& shape, HitInfo2D* hitInfo = nullptr);

// --- Bounding Box ---

bool intersect(const BBox2D& b1, const BBox2D& b2);

bool intersect(const BBox2D& bbox, const IPolygonalShape2D& polygon);

bool intersect(const BBox2D& bbox, const Circle2D& circle);

// --- Polygonal and Circle shapes ---

//bool intersect(const IPolygonalShape2D& polygon1, const IPolygonalShape2D& polygon2);

//bool intersect(const IPolygonalShape2D& polygon, const Circle2D& circle);

bool intersect(const Circle2D& circle1, const Circle2D& circle2);

// --- Generic ---

//bool intersect(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2);

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