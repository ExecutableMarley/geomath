/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <cmath>
#include <math.h>
#include <vector>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Utility
{

namespace Math
{

class Circle : public IFiniteShape2D
{
public:
    Vector2D m_center;
    float m_radius;

    Circle() : m_center(), m_radius(0) {}

    Circle(const Vector2D &center, float radius) : m_center(center), m_radius(radius) {}

    ShapeType2D type() const override
    {
        return SHAPE2D_CIRCLE;
    }

    std::vector<Vector2D> getVertices(int resolution = 32) const
    {
        std::vector<Vector2D> vertices;
        vertices.reserve(resolution);

        const float step = 2.0f * static_cast<float>(PI) / resolution;
        for (int i = 0; i < resolution; ++i)
        {
            float angle = i * step;
            float x = m_center.x + m_radius * std::cos(angle);
            float y = m_center.y + m_radius * std::sin(angle);
            vertices.emplace_back(x, y);
        }

        return vertices;
    }

    float area() const override
    {
        return PI * m_radius * m_radius;
    }

    float circumference() const
    {
        return 2.0f * PI * m_radius;
    }

    float perimeter() const override
    {
        return circumference();
    }

    Vector2D centroid() const override
    {
        return m_center;
    }

    Circle& translate(const Vector2D &translation) override
    {
        m_center += translation;
        return *this;
    }

    Circle& rotate(float angle)
    {
        return *this;
    }

    Circle& rotate(float angle, const Vector2D& point)
    {
        this->m_center.rotateAround(angle, point);
        return *this;
    }

    Circle& scale(float factor)
    {
        m_radius *= factor;
        return *this;
    }

    bool contains(const Vector2D &point) const override
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

    BBox2D boundingBox() const override
    {
        return BBox2D(m_center - Vector2D(m_radius, m_radius), m_center + Vector2D(m_radius, m_radius));
    }
};

} // namespace Math

} // namespace Utility