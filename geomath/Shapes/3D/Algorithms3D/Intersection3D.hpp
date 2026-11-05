/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <concepts>
#include <limits>

#include "Geometry/Vector3D.hpp"

namespace Arns::geomath
{

class Ray3D;
class BBox3D;
class Plane3D;
class Triangle3D;
class Sphere3D;
class Cylinder3D;
class Capsule3D;

struct HitInfo3D
{
    real_t t;
    Vector3D intersectionPoint;
    //Vector3D normal;
    //IShape3D* shape;
};

namespace detail
{

bool rayBoundsInterval(const Ray3D& ray, const BBox3D& bounds, real_t t_min, real_t t_max, real_t& enter, real_t& exit);
bool intersectRayWithBBox(const Ray3D& ray, const BBox3D& bbox, real_t t_min, real_t t_max, HitInfo3D* hitInfo = nullptr);
bool intersectRayWithSphere(const Ray3D& ray, const Sphere3D& sphere, real_t t_min, real_t t_max, HitInfo3D* hitInfo = nullptr);
bool intersectRayWithPlane(const Ray3D& ray, const Plane3D& plane, real_t t_min, real_t t_max, HitInfo3D* hitInfo = nullptr);
bool intersectRayWithTriangle(const Ray3D& ray, const Triangle3D& triangle, real_t t_min, real_t t_max, HitInfo3D* hitInfo = nullptr);
bool intersectRayWithCylinder(const Ray3D& ray, const Cylinder3D& cylinder, real_t t_min, real_t t_max, HitInfo3D* hitInfo = nullptr);
bool intersectRayWithCapsule(const Ray3D& ray, const Capsule3D& capsule, real_t t_min, real_t t_max, HitInfo3D* hitInfo = nullptr);

inline bool intersectRayImpl(const Ray3D& ray, const BBox3D& bbox, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    return intersectRayWithBBox(ray, bbox, tMin, tMax, hitInfo);
}

inline bool intersectRayImpl(const Ray3D& ray, const Sphere3D& sphere, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    return intersectRayWithSphere(ray, sphere, tMin, tMax, hitInfo);
}

inline bool intersectRayImpl(const Ray3D& ray, const Plane3D& plane, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    return intersectRayWithPlane(ray, plane, tMin, tMax, hitInfo);
}

inline bool intersectRayImpl(const Ray3D& ray, const Triangle3D& triangle, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    return intersectRayWithTriangle(ray, triangle, tMin, tMax, hitInfo);
}

inline bool intersectRayImpl(const Ray3D& ray, const Cylinder3D& cylinder, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    return intersectRayWithCylinder(ray, cylinder, tMin, tMax, hitInfo);
}

inline bool intersectRayImpl(const Ray3D& ray, const Capsule3D& capsule, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    return intersectRayWithCapsule(ray, capsule, tMin, tMax, hitInfo);
}

template <typename Object>
concept HasRayIntersectionImpl3D = requires(const Ray3D& ray, const Object& object, real_t tMin, real_t tMax, HitInfo3D* hitInfo)
{
    { intersectRayImpl(ray, object, tMin, tMax, hitInfo) } -> std::convertible_to<bool>;
};

}

//=======================================
// Public Intersection dispatch functions
//=======================================

template <typename Object> requires detail::HasRayIntersectionImpl3D<Object>
bool raycast(const Ray3D& ray, const Object& object,
             real_t tMin = real_t{0}, real_t tMax = std::numeric_limits<real_t>::max(),
             HitInfo3D* hitInfo = nullptr)
{
    return detail::intersectRayImpl(ray, object, tMin, tMax, hitInfo);
}

} // namespace Arns::geomath