/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Intersection3D.hpp"
#include "../Ray3D.hpp"
#include "../BBox3D.hpp"
#include "../Plane.hpp"
#include "../Triangle3D.hpp"
#include "../Sphere.hpp"
#include "../Cylinder.hpp"
#include "../Capsule.hpp"

#include <algorithm>

namespace Arns::geomath
{

bool detail::rayBoundsInterval(const Ray3D& ray, const BBox3D& bounds,
    real_t t_min, real_t t_max, real_t& enter, real_t& exit)
{
    if (t_min > t_max)
        return false;

    enter = t_min;
    exit = t_max;
    for (size_t axis = 0; axis < Vector3D::dimensions(); ++axis)
    {
        const real_t origin = ray.m_origin[axis];
        const real_t direction = ray.m_direction[axis];
        if (direction == real_t{0})
        {
            if (origin < bounds.m_min[axis] || origin > bounds.m_max[axis])
                return false;
            continue;
        }

        real_t first = (bounds.m_min[axis] - origin) / direction;
        real_t second = (bounds.m_max[axis] - origin) / direction;
        if (first > second)
            std::swap(first, second);

        enter = std::max(enter, first);
        exit = std::min(exit, second);
        if (enter > exit)
            return false;
    }
    return true;
}

bool detail::intersectRayWithBBox(const Ray3D& ray, const BBox3D& bbox, real_t t_min, real_t t_max, HitInfo3D* hitInfo)
{
    real_t enter;
    real_t exit;
    if (!detail::rayBoundsInterval(ray, bbox, t_min, t_max, enter, exit) || enter >= exit)
        return false;

    if (hitInfo)
    {
        hitInfo->t = enter;
        hitInfo->intersectionPoint = ray.pointAt(enter);
    }

    return true;
}

bool detail::intersectRayWithSphere(const Ray3D& ray, const Sphere3D& sphere, real_t t_min, real_t t_max, HitInfo3D* hitInfo)
{
    Vector3D oc = ray.m_origin - sphere.m_center;
    real_t a = ray.m_direction.lengthSquared();
    real_t b = oc.dot(ray.m_direction);
    real_t c = oc.lengthSquared() - sphere.m_radius * sphere.m_radius;
    real_t discriminant = b * b - a * c;

    if (discriminant > real_t{0})
    {
        real_t temp = (-b - sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min)
        {
            if (hitInfo)
            {
                hitInfo->t = temp;
                hitInfo->intersectionPoint = ray.pointAt(temp);
            }
            return true;
        }

        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min)
        {
            if (hitInfo)
            {
                hitInfo->t = temp;
                hitInfo->intersectionPoint = ray.pointAt(temp);
            }
            return true;
        }
    }
    return false;
}

bool detail::intersectRayWithPlane(const Ray3D& ray, const Plane3D& plane, real_t t_min, real_t t_max, HitInfo3D* hitInfo)
{
    real_t denominator = plane.m_normal.dot(ray.m_direction);

    if (approximatelyZero(denominator))
        return false;

    real_t t = (plane.m_normal.dot(plane.m_planePoint) - plane.m_normal.dot(ray.m_origin)) / denominator;

    if (t > t_min && t < t_max)
    {
        if (hitInfo)
        {
            hitInfo->t = t;
            hitInfo->intersectionPoint = ray.pointAt(t);
        }
        return true;
    }
    return false;
}

bool detail::intersectRayWithTriangle(const Ray3D& ray, const Triangle3D& triangle, real_t t_min, real_t t_max, HitInfo3D* hitInfo)
{
    Vector3D edge1 = triangle.m_b - triangle.m_a;
    Vector3D edge2 = triangle.m_c - triangle.m_a;
    Vector3D h = ray.m_direction.cross(edge2);
    real_t a = edge1.dot(h);

    if (approximatelyZero(a))
        return false;

    real_t f = real_t{1} / a;
    Vector3D s = ray.m_origin - triangle.m_a;
    real_t u = f * s.dot(h);

    if (u < real_t{0} || u > real_t{1})
        return false;

    Vector3D q = s.cross(edge1);
    real_t v = f * ray.m_direction.dot(q);

    if (v < real_t{0} || u + v > real_t{1})
        return false;

    real_t t = f * edge2.dot(q);

    if (t > t_min && t < t_max)
    {
        if (hitInfo)
        {
            hitInfo->t = t;
            hitInfo->intersectionPoint = ray.pointAt(t);
        }
        return true;
    }
    return false;
}

bool detail::intersectRayWithCylinder(const Ray3D& ray, const Cylinder3D& cylinder, real_t t_min, real_t t_max, HitInfo3D* hitInfo)
{
    Vector3D d = cylinder.m_endPoint - cylinder.m_startPoint;
    Vector3D m = ray.m_origin - cylinder.m_startPoint;
    Vector3D n = ray.m_direction;

    real_t dd = d.dot(d);
    real_t nd = n.dot(d);
    real_t md = m.dot(d);

    real_t a = dd - nd * nd;
    real_t b = dd * m.dot(n) - md * nd;
    real_t c = dd * m.dot(m) - md * md - cylinder.m_radius * cylinder.m_radius * dd;

    real_t discriminant = b * b - a * c;

    if (discriminant < real_t{0})
        return false;

    real_t sqrtDiscriminant = sqrt(discriminant);
    real_t t0 = (-b - sqrtDiscriminant) / a;
    real_t t1 = (-b + sqrtDiscriminant) / a;

    if (t0 > t1)
        std::swap(t0, t1);

    if (t1 < t_min || t0 > t_max)
        return false;

    real_t t = t0;
    if (t < t_min)
        t = t1;
    if (t < t_min || t > t_max)
        return false;

    Vector3D hitPoint = ray.pointAt(t);
    Vector3D projOnCylinderAxis = cylinder.m_startPoint + d * ((hitPoint - cylinder.m_startPoint).dot(d) / dd);

    if ((projOnCylinderAxis - cylinder.m_startPoint).dot(d) < real_t{0} ||
        (projOnCylinderAxis - cylinder.m_endPoint).dot(d) > real_t{0})
        return false;

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = hitPoint;
    }

    return true;
}

bool intersectRayWithCapsuleOptimized(const Ray3D& ray, const Capsule3D& capsule, real_t t_min, real_t t_max)
{
#if 0
    Vector3D startPoint = ray.pointAt(t_min);
    Vector3D endPoint = ray.pointAt(t_max);
    real_t distance = detail::distanceSegmentToSegment(startPoint, endPoint, capsule.m_startPoint, capsule.m_endPoint);
    return distance <= capsule.m_radius;
#else
    Sphere3D sphere1(capsule.m_startPoint, capsule.m_radius);
    Sphere3D sphere2(capsule.m_endPoint, capsule.m_radius);
    Cylinder3D cylinder(capsule.m_startPoint, capsule.m_endPoint, capsule.m_radius);

    if (detail::intersectRayWithSphere(ray, sphere1, t_min, t_max, nullptr))
        return true;
    if (detail::intersectRayWithSphere(ray, sphere2, t_min, t_max, nullptr))
        return true;
    if (detail::intersectRayWithCylinder(ray, cylinder, t_min, t_max, nullptr))
        return true;

    return false;
#endif
}

bool detail::intersectRayWithCapsule(const Ray3D& ray, const Capsule3D& capsule, real_t t_min, real_t t_max, HitInfo3D* hitInfo)
{
    if (hitInfo == nullptr) 
        return intersectRayWithCapsuleOptimized(ray, capsule, t_min, t_max);

    Sphere3D sphere1(capsule.m_startPoint, capsule.m_radius);
    Sphere3D sphere2(capsule.m_endPoint, capsule.m_radius);
    Cylinder3D cylinder(capsule.m_startPoint, capsule.m_endPoint, capsule.m_radius);
    HitInfo3D localHitInfo;
    HitInfo3D* output = hitInfo ? hitInfo : &localHitInfo;

    bool hit = false;
    if (intersectRayWithSphere(ray, sphere1, t_min, t_max, output))
    {
        t_max = output->t;
        hit = true;
    }
    if (intersectRayWithSphere(ray, sphere2, t_min, t_max, output))
    {
        t_max = output->t;
        hit = true;
    }
    if (intersectRayWithCylinder(ray, cylinder, t_min, t_max, output))
        hit = true;

    return hit;
}

} // namespace Arns::geomath