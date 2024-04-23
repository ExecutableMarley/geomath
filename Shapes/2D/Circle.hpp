#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Vector2D.hpp"
#include "BBox2D.hpp"

namespace Utility
{

namespace Math
{

struct Circle
{
    Vector2D m_center;
    float m_radius;

    Circle() : m_center(), m_radius(0) {}

    Circle(const Vector2D &center, float radius) : m_center(center), m_radius(radius) {}

    float area() const
    {
        return PI * m_radius * m_radius;
    }

    float circumference() const
    {
        return 2.0f * PI * m_radius;
    }

    float perimeter() const
    {
        return circumference();
    }

    Vector2D centroid() const
    {
        return m_center;
    }

    Circle& translate(const Vector2D &translation)
    {
        m_center += translation;
        return *this;
    }

    Circle& rotate(float angle)
    {
        return *this;
    }

    bool contains(const Vector2D &point) const
    {
        return (point - m_center).lengthSquared() <= m_radius * m_radius;
    }

    bool contains(const Circle &circle) const
    {
        return (m_center - circle.m_center).lengthSquared() + circle.m_radius * circle.m_radius <= m_radius * m_radius;
    }

    bool intersects(const Circle &circle) const
    {
        return (m_center - circle.m_center).lengthSquared() <= (m_radius + circle.m_radius) * (m_radius + circle.m_radius);
    }

    BBox2D boundingBox() const
    {
        return BBox2D(m_center - Vector2D(m_radius, m_radius), m_center + Vector2D(m_radius, m_radius));
    }
};

} // namespace Math

} // namespace Utility