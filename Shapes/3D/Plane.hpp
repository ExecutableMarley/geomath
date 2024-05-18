/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector3D.hpp"

namespace Utility
{

namespace Math
{

class Plane
{
public:
    Vector3D m_planePoint;
    Vector3D m_normal;

    Plane() : m_planePoint(), m_normal() {}

    Plane(const Vector3D &planePoint, const Vector3D &normal) : m_planePoint(planePoint), m_normal(normal) {}

    Plane(const Vector3D &point1, const Vector3D &point2, const Vector3D &point3)
    {
        Vector3D v1 = point2 - point1;
        Vector3D v2 = point3 - point1;
        m_normal = v1.cross(v2).normalize();
        m_planePoint = point1;
    }

    float distance(const Vector3D &point) const
    {
        return m_normal.dot(point - m_planePoint);
    }

    Vector3D project(const Vector3D &point) const
    {
        return point - m_normal * distance(point);
    }

    bool isParallel(const Plane &plane) const
    {
        return m_normal.isParallel(plane.m_normal);
    }
};



} // namespace Math

} // namespace Utility