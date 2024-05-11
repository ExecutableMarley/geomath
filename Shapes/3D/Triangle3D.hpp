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

class Triangle3D
{
public:
    Vector3D m_a;
    Vector3D m_b;
    Vector3D m_c;

    Triangle3D() : m_a(), m_b(), m_c() {}

    Triangle3D(const Vector3D &a, const Vector3D &b, const Vector3D &c) : m_a(a), m_b(b), m_c(c) {}

    float volume() const
    {
        return 0.f;
    }

    float surfaceArea() const
    {
        return 0.5f * (m_b - m_a).cross(m_c - m_a).length();
    }

    Vector3D centroid() const
    {
        return (m_a + m_b + m_c) / 3.0f;
    }

    Triangle3D& translate(const Vector3D &translation)
    {
        m_a += translation;
        m_b += translation;
        m_c += translation;
        return *this;
    }

    Vector3D normal() const
    {
        return (m_b - m_a).cross(m_c - m_a).normalize();
    }

    //Todo: Problematic
    bool contains(const Vector3D &point) const
    {
        Vector3D n = normal();
        Vector3D ab = m_b - m_a;
        Vector3D bc = m_c - m_b;
        Vector3D ca = m_a - m_c;

        Vector3D ap = point - m_a;
        Vector3D bp = point - m_b;
        Vector3D cp = point - m_c;

        if (n.dot(ab.cross(ap)) > 0.0f && n.dot(bc.cross(bp)) > 0.0f && n.dot(ca.cross(cp)) > 0.0f)
        {
            return true;
        }
        return false;
    }
};

} // namespace Math

} // namespace Utility