/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "BBox3D.hpp"
#include "IShape3D.hpp"

namespace Arns
{

namespace Math
{

class Cylinder : public IFiniteShape3D
{
public:
    Vector3D m_startPoint;
    Vector3D m_endPoint;
    real_t m_radius;

    Cylinder() : m_startPoint(), m_endPoint(), m_radius(0) {}

    Cylinder(const Vector3D &startPoint, const Vector3D &endPoint, real_t radius) : m_startPoint(startPoint), m_endPoint(endPoint), m_radius(radius) {}

    ShapeType3D type() const override
    {
        return ShapeType3D::CYLINDER;
    }

    real_t volume() const override
    {
        return PI * m_radius * m_radius * (m_endPoint - m_startPoint).length();
    }

    real_t surfaceArea() const override
    {
        return 2.0f * PI * m_radius * (m_endPoint - m_startPoint).length() + 2.0f * PI * m_radius * m_radius;
    }

    Vector3D centroid() const override
    {
        return (m_startPoint + m_endPoint) * 0.5f;
    }

    bool contains(const Vector3D &point) const override
    {
        // Todo: Calculating the closest parameter t should be done in a separate function
        const Vector3D delta = m_endPoint - m_startPoint;
        real_t t = (point - m_startPoint).dot(delta) / (delta).lengthSquared();

        Vector3D closestPoint = m_startPoint + delta * t;
        if (t >= 0.0f && t <= 1.0f && (point - closestPoint).lengthSquared() <= m_radius * m_radius)
        {
            return true;
        }
        return false;
    }

    Cylinder& translate(const Vector3D &translation) override
    {
        m_startPoint += translation;
        m_endPoint += translation;
        return *this;
    }

    Cylinder& scale(real_t scaleFactor)
    {
        Vector3D centroid = this->centroid();
        m_startPoint = centroid + scaleFactor * (m_startPoint - centroid);
        m_endPoint   = centroid + scaleFactor * (m_endPoint - centroid);
        m_radius    *= scaleFactor;
        return *this;
    }

    BBox3D boundingBox() const override
    {
        const Vector3D min = Vector3D::min(m_startPoint, m_endPoint) - Vector3D(m_radius, m_radius, m_radius);
        const Vector3D max = Vector3D::max(m_startPoint, m_endPoint) + Vector3D(m_radius, m_radius, m_radius);
        return BBox3D(min, max);
    }
};


} // namespace Math

} // namespace Arns