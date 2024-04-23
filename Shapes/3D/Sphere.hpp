#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Vector3D.hpp"
#include "BBox3D.hpp"

namespace Utility
{

namespace Math
{

class Sphere
{
public:
    Vector3D m_center;
    float m_radius;

    Sphere() : m_center(), m_radius(0) {}

    Sphere(const Vector3D &center, float radius) : m_center(center), m_radius(radius) {}

    float volume() const
    {
        return 4.0f / 3.0f * PI * m_radius * m_radius * m_radius;
    }

    bool contains(const Vector3D &point) const
    {
        return (point - m_center).lengthSquared() <= m_radius * m_radius;
    }

    bool contains(const Sphere &sphere) const
    {
        return (m_center - sphere.m_center).lengthSquared() + sphere.m_radius * sphere.m_radius <= m_radius * m_radius;
    }
    
    bool intersects(const Sphere &sphere) const
    {
        return (m_center - sphere.m_center).lengthSquared() <= (m_radius + sphere.m_radius) * (m_radius + sphere.m_radius);
    }

    BBox3D boundingBox() const
    {
        return BBox3D(m_center - Vector3D(m_radius, m_radius, m_radius), m_center + Vector3D(m_radius, m_radius, m_radius));
    }
};

} // namespace Math

} // namespace Utility