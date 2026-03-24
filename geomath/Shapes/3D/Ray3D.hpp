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

class Ray3D
{
public:
    Vector3D m_origin;
    Vector3D m_direction;

    Ray3D() : m_origin(), m_direction() {}

    Ray3D(const Vector3D &origin, const Vector3D &direction) : m_origin(origin), m_direction(direction) {}

    Vector3D origin() const
    {
        return m_origin;
    }

    Vector3D direction() const
    {
        return m_direction;
    }

    Vector3D pointAt(real_t t) const
    {
        return m_origin + m_direction * t;
    }

    real_t closestParameter(const Vector3D &point) const
    {
        return (point - m_origin).dot(m_direction) / m_direction.lengthSquared();
    }
};


} // namespace Math

} // namespace Arns