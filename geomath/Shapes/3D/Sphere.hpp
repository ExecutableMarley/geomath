/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "BBox3D.hpp"
#include "IShape3D.hpp"

namespace Arns
{

namespace geomath
{

class Sphere3D : public IFiniteShape3D
{
public:
    Vector3D m_center;
    real_t m_radius;

    Sphere3D() : m_center(), m_radius(0) {}

    Sphere3D(const Vector3D &center, real_t radius) : m_center(center), m_radius(radius) {}

    ShapeType3D type() const override
    {
        return ShapeType3D::SPHERE;
    }

    real_t volume() const override
    {
        return 4.0f / 3.0f * PI * m_radius * m_radius * m_radius;
    }

    real_t surfaceArea() const override
    {
        return 4.0f * PI * m_radius * m_radius;
    }

    Vector3D centroid() const override
    {
        return m_center;
    }

    bool contains(const Vector3D &point) const override
    {
        return (point - m_center).lengthSquared() <= m_radius * m_radius;
    }

    bool contains(const Sphere3D &sphere) const
    {
        return (m_center - sphere.m_center).lengthSquared() + sphere.m_radius * sphere.m_radius <= m_radius * m_radius;
    }

    Sphere3D& translate(const Vector3D &translation) override
    {
        m_center += translation;
        return *this;
    }

    Sphere3D& scale(real_t scaleFactor)
    {
        m_radius *= scaleFactor;
        return *this;
    }
    
    bool intersects(const Sphere3D &sphere) const
    {
        return (m_center - sphere.m_center).lengthSquared() <= (m_radius + sphere.m_radius) * (m_radius + sphere.m_radius);
    }

    BBox3D boundingBox() const override
    {
        return BBox3D(m_center - Vector3D(m_radius, m_radius, m_radius), m_center + Vector3D(m_radius, m_radius, m_radius));
    }
};

} // namespace Math

} // namespace Arns