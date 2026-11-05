/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"

namespace Arns
{

namespace geomath
{

class Plane3D
{
public:
    Vector3D m_planePoint;
    Vector3D m_normal;

    Plane3D() : m_planePoint(), m_normal() {}

    Plane3D(const Vector3D &planePoint, const Vector3D &normal) : m_planePoint(planePoint), m_normal(normal) {}

    Plane3D(const Vector3D &point1, const Vector3D &point2, const Vector3D &point3)
    {
        Vector3D v1 = point2 - point1;
        Vector3D v2 = point3 - point1;
        m_normal = v1.cross(v2).normalize();
        m_planePoint = point1;
    }

    real_t distance(const Vector3D &point) const
    {
        return m_normal.dot(point - m_planePoint);
    }

    Vector3D project(const Vector3D &point) const
    {
        return point - m_normal * distance(point);
    }

    bool isParallel(const Plane3D &plane) const
    {
        return m_normal.isParallel(plane.m_normal);
    }

    bool isOrthogonal(const Plane3D &plane) const
    {
        return m_normal.isOrthogonal(plane.m_normal);
    }

    bool isOnPlane(const Vector3D &point) const
    {
        return approximatelyZero(distance(point));
    }

    bool contains(const Vector3D &point) const
    {
        return isOnPlane(point);
    }

    bool intersects(const Plane3D &plane) const
    {
        return !isParallel(plane);
    }

    static Plane3D fromAxisNormal(const Vector3D &point, int axis)
    {
        Vector3D normal(0,0,0);
        normal[axis] = 1.f;
        return Plane3D(point, normal);
    }
};



} // namespace Math

} // namespace Arns