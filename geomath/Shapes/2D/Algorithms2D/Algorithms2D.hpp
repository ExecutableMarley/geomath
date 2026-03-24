/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"

#include <math.h>
#include <vector>
#include <span>

namespace Arns
{

namespace Math
{

// Forward declarations
class Ray2D;
class Line2D;
class BBox2D;
class IBaseShape2D;
class IFiniteShape2D;
class IPolygonalShape2D;
class Circle2D;
class Triangle2D;
class Rectangle2D;
class Polygon2D;
class ConvexPolygon2D;

//Consider renaming Line2D in the future
//Segment2D is better
using Segment2D = Line2D;

struct HitInfo2D
{
    real_t t;
    Vector2D intersectionPoint;
    //Todo: Collect more information
    //Vector2D normal;
    //IShape2D* shape;
};


//Todo: Containment functions

//

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd);

bool isPointOnSegment(const Vector2D& point, const Line2D& line);

bool isSegmentOnSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2);

bool isSegmentOnSegment(const Segment2D& segment1, const Segment2D& segment2);

// --- Point-in-Shape containment

bool isPointInsideBBox(const Vector2D& point, const BBox2D& bbox);

bool isPointInsideTriangle(const Vector2D& point, const Triangle2D& triangle);

bool isPointInsideRectangle(const Vector2D& point, const Rectangle2D& rectangle);

bool isPointInsideConvexPolygon(const Vector2D& point, const ConvexPolygon2D& polygon);

bool isPointInsidePolygon(const Vector2D& point, const Polygon2D& polygon);

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle);

// --- Segment-in-Shape containment

bool isSegmentInsideBBox(const Segment2D& segment, const BBox2D& bbox);

bool isSegmentInsideTriangle(const Segment2D& segment, const Triangle2D& triangle);

bool isSegmentInsideRectangle(const Segment2D& segment, const Rectangle2D& rectangle);

bool isSegmentInsideConvexPolygon(const Segment2D& segment, const ConvexPolygon2D& polygon);

bool isSegmentInsidePolygon(const Segment2D& segment, const Polygon2D polygon);

bool isSegmentInsideCircle(const Segment2D& segment, const Circle2D& circle);

//

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2);

bool isPolygonalInsideBBox(const IPolygonalShape2D& polygon, const BBox2D& bbox);

bool isCircleInsideBBox(const Circle2D& circle, const BBox2D& bbox);

//

bool isShapeInsideConvexPolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> convexPoly);

bool isShapeInsidePolygon(const IFiniteShape2D& shape, const std::span<const Vector2D> poly);

bool isShapeInsideCircle(const IFiniteShape2D& shape, const Circle2D& circle);

bool isShapeInsideShape(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2);


// Distance calculation algorithms

real_t distancePointToSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd, Vector2D* closestPoint = nullptr);

real_t distanceSegmentToSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

//Todo: Distance to shapes

real_t distance(const Vector2D& point, const Segment2D& segment, Vector2D* closestPoint = nullptr);

real_t distance(const Vector2D& point, const BBox2D& bbox, Vector2D* closestPoint = nullptr);

real_t distance(const Vector2D& point, const IPolygonalShape2D& polygon, Vector2D* closestPoint = nullptr);

real_t distance(const Vector2D& point, const Circle2D& circle, Vector2D* closestPoint = nullptr);

real_t distance(const Segment2D& segment1, const Segment2D& segment2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

real_t distance(const Segment2D& segment, const BBox2D& bbox, Vector2D* closestPoint = nullptr);

real_t distance(const Segment2D& segment, const IPolygonalShape2D& polygon, Vector2D* closestPoint = nullptr);

real_t distance(const Segment2D& segment, const Circle2D& circle, Vector2D* closestPoint = nullptr);

//

real_t distance(const BBox2D& a, const BBox2D& b, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const BBox2D& bbox1, const IPolygonalShape2D& poly2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const BBox2D& bbox, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const IPolygonalShape2D& poly1, const IPolygonalShape2D& poly2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const IPolygonalShape2D& poly1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distance(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

// Intersection calculation algorithms

//    --- Rays ---

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const Segment2D& segment, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const BBox2D& bbox, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const Triangle2D& triangle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const IPolygonalShape2D& polygon, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const Circle2D& circle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersect(const Ray2D& ray, const IBaseShape2D& shape, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

//    --- Segments ---

bool intersectSegmentWithSegment(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersect(const Line2D& line1, const Line2D& line2, HitInfo2D* hitInfo = nullptr);

bool intersect(const Line2D& line, const BBox2D& bbox, HitInfo2D* hitInfo = nullptr);

bool intersect(const Line2D& line, const Circle2D& circle, HitInfo2D* hitInfo = nullptr);

bool intersect(const Line2D& line, const Triangle2D& triangle, HitInfo2D* hitInfo = nullptr);

bool intersect(const Line2D& line, const IPolygonalShape2D& polygon, HitInfo2D* hitInfo = nullptr);

bool intersect(const Segment2D& segment, const IBaseShape2D& shape, HitInfo2D* hitInfo = nullptr);

// --- Bounding Box ---

bool intersect(const BBox2D& b1, const BBox2D& b2);

bool intersect(const BBox2D& bbox, const IPolygonalShape2D& polygon);

bool intersect(const BBox2D& bbox, const Circle2D& circle);

// --- Polygonal and Circle shapes ---

bool intersect(const IPolygonalShape2D& polygon1, const IPolygonalShape2D& polygon2);

bool intersect(const IPolygonalShape2D& polygon, const Circle2D& circle);

bool intersect(const Circle2D& circle1, const Circle2D& circle2);

// --- Generic ---

bool intersect(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2);

} // namespace Math

} // namespace Arns