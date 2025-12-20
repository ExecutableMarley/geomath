/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"

namespace Arns
{

namespace Math
{

class Ray2D
{
public:
    Vector2D m_origin;
    Vector2D m_direction;

    Ray2D() : m_origin(), m_direction() {}

    Ray2D(const Vector2D &origin, const Vector2D &direction) : m_origin(origin), m_direction(direction) {}

    Vector2D origin() const
    {
        return m_origin;
    }

    Vector2D direction() const
    {
        return m_direction;
    }

    Vector2D pointAt(real_t t) const
    {
        return m_origin + m_direction * t;
    }

    real_t closestParameter(const Vector2D &point) const
    {
        return (point - m_origin).dot(m_direction) / m_direction.lengthSquared();
    }
};

} // namespace Math

} // namespace Arns