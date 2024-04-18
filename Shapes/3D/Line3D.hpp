#pragma once

#include <math.h>

#include "Vector3D.hpp"

namespace Utility
{

namespace Math
{

struct Line3D
{
    Vector3D m_start;
    Vector3D m_end;

    Line3D() : m_start(), m_end() {}

    Line3D(const Vector3D &start, const Vector3D &end) : m_start(start), m_end(end) {}

    Line3D(const Vector3D &start, const Vector3D &direction, float length) : m_start(start), m_end(start + direction.createNormalized() * length) {}

    float length() const
    {
        return (m_end - m_start).length();
    }

    Vector3D direction() const
    {
        return (m_end - m_start).normalize();
    }
};

} // namespace Math

} // namespace Utility