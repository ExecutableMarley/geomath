/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"

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

// Distance calculation algorithms

real_t distancePointToLine(const Vector2D& point, const Vector2D& lineStart, const Vector2D& lineEnd, Vector2D* closestPoint = nullptr);

real_t distancePointToLine(const Vector2D& point, const Line2D& line, Vector2D* closestPoint = nullptr);

real_t distanceLineToLine(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

real_t distanceLineToLine(const Line2D& line1, const Line2D& line2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

//Todo: Distance to shapes

real_t distancePointToBBox(const Vector2D& point, const BBox2D& bbox, Vector2D* closestPoint = nullptr);

real_t distancePointToTriangle(const Vector2D& point, const Triangle2D& triangle, Vector2D* closestPoint = nullptr);

real_t distancePointToRectangle(const Vector2D& point, const Rectangle2D& rectangle, Vector2D* closestPoint = nullptr);

real_t distancePointToPolygon(const Vector2D& point, const ConvexPolygon2D& polygon, Vector2D* closestPoint = nullptr);

real_t distancePointToPolygon(const Vector2D& point, const Polygon2D& polygon, Vector2D* closestPoint = nullptr);

real_t distancePointToCircle(const Vector2D& point, const Circle2D& circle, Vector2D* closestPoint = nullptr);


// Intersection calculation algorithms

//    --- Rays ---

bool intersectRayWithBBox(const Ray2D& ray, const BBox2D& bbox, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersectRayWithCircle(const Ray2D& ray, const Circle2D& circle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersectRayWithTriangle(const Ray2D& ray, const Triangle2D& triangle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersectRayWithRectangle(const Ray2D& ray, const Rectangle2D& rectangle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersectRayWithPolygon(const Ray2D& ray, const Polygon2D& polygon, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

bool intersectRayWithShape(const Ray2D& ray, const IBaseShape2D& shape, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D* hitInfo = nullptr);

//    --- Segments ---

bool intersectSegmentWithSegment(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegment(const Line2D& line1, const Line2D& line2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegmentStrict(const Line2D& line1, const Line2D& line2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithBBox(const Line2D& line, const BBox2D& bbox, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithCircle(const Line2D& line, const Circle2D& circle, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithTriangle(const Line2D& line, const Triangle2D& triangle, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithRectangle(const Line2D& line, const Rectangle2D& rectangle, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithPolygon(const Line2D& line, const Polygon2D& polygon, HitInfo2D* hitInfo = nullptr);

// --- Bounding Box ---

bool intersectBBoxWithBBox(const BBox2D& b1, const BBox2D& b2);

bool intersectBBoxWithTriangle(const BBox2D& bbox, const Triangle2D& triangle);

bool intersectBBoxWithRectangle(const BBox2D& bbox, const Rectangle2D& rectangle);

bool intersectBBoxWithPolygon(const BBox2D& bbox, const Polygon2D& polygon);

bool intersectBBoxWithCircle(const BBox2D& bbox, const Circle2D& circle);

// --- Triangles ---

bool intersectTriangleWithTriangle(const Triangle2D& triangle1, const Triangle2D& triangle2);

bool intersectTriangleWithRectangle(const Triangle2D& triangle, const Rectangle2D& rectangle);

bool intersectTriangleWithPolygon(const Triangle2D& triangle, const ConvexPolygon2D& polygon);

bool intersectTriangleWithCircle(const Triangle2D& triangle, const Circle2D& circle);

// --- Rectangles ---

bool intersectRectangleWithRectangle(const Rectangle2D& rectangle1, const Rectangle2D& rectangle2);

bool intersectRectangleWithPolygon(const Rectangle2D& rectangle, const ConvexPolygon2D& polygon);

bool intersectRectangleWithCircle(const Rectangle2D& rectangle, const Circle2D& circle);

// --- Polygons --- 

bool intersectConvexPolygonWithConvexPolygon(const std::vector<Vector2D>& p1, const std::vector<Vector2D>& p2);

bool intersectPolygonWithPolygon(const ConvexPolygon2D& polygon1, const ConvexPolygon2D& polygon2);

bool intersectPolygonWithCircle(const ConvexPolygon2D& polygon, const Circle2D& circle);

// --- Circles ---

bool intersectCircleWithCircle(const Circle2D& circle1, const Circle2D& circle2);

// --- Generic ---

bool intersect(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2);

} // namespace Math

} // namespace Arns