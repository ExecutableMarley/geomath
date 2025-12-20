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

class Line2D
{
public:
    Vector2D m_start;
    Vector2D m_end;

    Line2D() : m_start(), m_end() {}

    Line2D(Vector2D startPoint, Vector2D endPoint) : m_start(startPoint), m_end(endPoint) {}

    Line2D(Vector2D startPoint, Vector2D direction, real_t length) : m_start(startPoint), m_end(startPoint + direction.createNormalized() * length) {}

    real_t length() const
    {
        return (m_start - m_end).length();
    }

    Vector2D direction() const
    {
        return (m_end - m_start).normalize();
    }

    Vector2D deltaVector() const
    {
        return m_end - m_start;
    }

    Vector2D pointAt(real_t t) const
    {
        return m_start + (m_end - m_start) * t;
    }

    Vector2D normal() const
    {
        return Vector2D(m_end.y - m_start.y, m_start.x - m_end.x).normalize();
    }

    real_t closestParameter(const Vector2D &point) const
    {
        Vector2D deltaVector = m_end - m_start;
        const real_t lengthSquared = deltaVector.lengthSquared();
        return lengthSquared == 0.f ? 0.f : (point - m_start).dot(deltaVector) / lengthSquared;
    }
};

} // namespace Math

} // namespace Arns