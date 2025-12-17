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

    Triangle3D& scale(float scaleFactor)
    {
        Vector3D centroid = this->centroid();
        m_a = centroid + scaleFactor * (m_a - centroid);
        m_b = centroid + scaleFactor * (m_b - centroid);
        m_c = centroid + scaleFactor * (m_c - centroid);
        return *this;
    }

    Vector3D normal() const
    {
        return (m_b - m_a).cross(m_c - m_a).normalize();
    }

    bool contains(const Vector3D &point) const
    {
        Vector3D ab = m_b - m_a;
        Vector3D bc = m_c - m_b;
        Vector3D ca = m_a - m_c;

        Vector3D ap = point - m_a;
        Vector3D bp = point - m_b;
        Vector3D cp = point - m_c;

        Vector3D crossAB_AP = ab.cross(ap);
        Vector3D crossBC_BP = bc.cross(bp);
        Vector3D crossCA_CP = ca.cross(cp);

        float cross1 = crossAB_AP.dot(crossBC_BP);
        float cross2 = crossBC_BP.dot(crossCA_CP);
        float cross3 = crossCA_CP.dot(crossAB_AP);

        if ((cross1 >= 0 && cross2 >= 0 && cross3 >= 0) 
            || (cross1 <= 0 && cross2 <= 0 && cross3 <= 0))
        {
            return true;
        }
        return false;
    }
};

} // namespace Math

} // namespace Arns