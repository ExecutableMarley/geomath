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

class Line3D
{
public:
    Vector3D m_start;
    Vector3D m_end;

    Line3D() : m_start(), m_end() {}

    Line3D(const Vector3D &start, const Vector3D &end) : m_start(start), m_end(end) {}

    Line3D(const Vector3D &start, const Vector3D &direction, real_t length) : m_start(start), m_end(start + direction.createNormalized() * length) {}

    real_t length() const
    {
        return (m_end - m_start).length();
    }

    Vector3D direction() const
    {
        return (m_end - m_start).normalize();
    }

    Vector3D deltaVector() const
    {
        return m_end - m_start;
    }

    Vector3D pointAt(real_t t) const
    {
        return m_start + (m_end - m_start) * t;
    }

    real_t closestParameter(const Vector3D &point) const
    {
        Vector3D deltaVector = m_end - m_start;
        return (point - m_start).dot(deltaVector) / (deltaVector).lengthSquared();
    }
};

} // namespace Math

} // namespace Arns