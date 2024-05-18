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
    Vector3D planePoint;
    Vector3D normal;

    Plane() : planePoint(), normal() {}

    Plane(const Vector3D &planePoint, const Vector3D &normal) : planePoint(planePoint), normal(normal) {}

    Plane(const Vector3D &point1, const Vector3D &point2, const Vector3D &point3)
    {
        Vector3D v1 = point2 - point1;
        Vector3D v2 = point3 - point1;
        normal = v1.cross(v2).normalize();
        planePoint = point1;
    }

    float distance(const Vector3D &point) const
    {
        return normal.dot(point - planePoint);
    }

    Vector3D project(const Vector3D &point) const
    {
        return point - normal * distance(point);
    }

    bool isParallel(const Plane &plane) const
    {
        return normal.isParallel(plane.normal);
    }
};



} // namespace Math

} // namespace Utility