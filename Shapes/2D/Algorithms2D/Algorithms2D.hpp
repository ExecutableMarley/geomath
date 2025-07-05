/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector2D.hpp"
//#include "../Line2D.hpp"

namespace Utility
{

namespace Math
{

// Forward declarations
class Ray2D;
class Line2D;
class BBox2D;
class IShape2D;
class Circle;
class Triangle; // rename?
class Rectangle;
class Polygon;

struct HitInfo2D
{
    float t;
    Vector2D intersectionPoint;
    //Vector2D normal;
    //IShape2D* shape;
};

//

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd);

bool isPointOnSegment(const Vector2D& point, const Line2D& line);

// Distance calculation algorithms

float distancePointToLine(const Vector2D& point, const Vector2D& lineStart, const Vector2D& lineEnd, Vector2D* closestPoint = nullptr);

float distancePointToLine(const Vector2D& point, const Line2D& line, Vector2D* closestPoint = nullptr);

float distanceLineToLine(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

float distanceLineToLine(const Line2D& line1, const Line2D& line2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

// Intersection calculation algorithms

bool intersectRayWithBBox(const Ray2D& ray, const BBox2D& bbox, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectRayWithCircle(const Ray2D& ray, const Circle& circle, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectRayWithTriangle(const Ray2D& ray, const Triangle& triangle, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectRayWithRectangle(const Ray2D& ray, const Rectangle& rectangle, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectRayWithPolygon(const Ray2D& ray, const Polygon& polygon, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectRayWithShape(const Ray2D& ray, const IShape2D& shape, float t_min, float t_max, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithSegmentStrict(const Line2D& line1, const Line2D& line2, HitInfo2D* hitInfo = nullptr);

bool intersectSegmentWithPolygon(const Line2D& line, const Polygon& polygon, HitInfo2D* hitInfo = nullptr);
} // namespace Math

} // namespace Utility