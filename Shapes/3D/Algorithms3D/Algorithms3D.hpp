/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector3D.hpp"

namespace Arns
{

namespace Math
{

// Forward declarations

class Ray3D;
class Line3D;
class BBox3D;
class ISHape3D;
class Plane;
class Triangle3D;
class Sphere;
class Cylinder;
class Capsule;

struct HitInfo3D
{
    real_t t;
    Vector3D intersectionPoint;
    //Vector3D normal;
    //IShape3D* shape;
};

// Distance calculation algorithms

real_t distancePointToLine(const Vector3D& point, const Vector3D& lineStart, const Vector3D& lineEnd, Vector3D* closestPoint = nullptr);

real_t distancePointToLine(const Vector3D& point, const Line3D& line, Vector3D* closestPoint = nullptr);

real_t distanceLineToLine(const Vector3D& s1, const Vector3D& s2, const Vector3D& k1, const Vector3D& k2, Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);

real_t distanceLineToLine(const Line3D& line1, const Line3D& line2, Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);

// Intersection calculation algorithms

bool intersectRayWithBBox(const Ray3D& ray, const BBox3D& bbox, real_t t_min, real_t t_max, HitInfo3D* hitInfo = NULL);

bool intersectRayWithSphere(const Ray3D& ray, const Sphere& sphere, real_t t_min, real_t t_max, HitInfo3D* hitInfo = NULL);

bool intersectRayWithPlane(const Ray3D& ray, const Plane& plane, real_t t_min, real_t t_max, HitInfo3D* hitInfo = NULL);

bool intersectRayWithTriangle(const Ray3D& ray, const Triangle3D& triangle, real_t t_min, real_t t_max, HitInfo3D* hitInfo = NULL);

bool intersectRayWithCylinder(const Ray3D& ray, const Cylinder& cylinder, real_t t_min, real_t t_max, HitInfo3D* hitInfo = NULL);

bool intersectRayWithCapsule(const Ray3D& ray, const Capsule& capsule, real_t t_min, real_t t_max);
//Todo: Doesn't handle null
bool intersectRayWithCapsule(const Ray3D& ray, const Capsule& capsule, real_t t_min, real_t t_max, HitInfo3D* hitInfo = NULL);

} // namespace Math

} // namespace Arns