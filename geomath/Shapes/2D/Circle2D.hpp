/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <cmath>
#include <math.h>
#include <vector>
#include <memory>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Arns
{

namespace Math
{

class Circle2D : public IFiniteShape2D
{
public:
    Vector2D m_center;
    real_t m_radius;

    Circle2D() : m_center(), m_radius(0) {}

    Circle2D(const Vector2D &center, real_t radius) : m_center(center), m_radius(radius) {}

    ShapeType2D type() const override
    {
        return SHAPE2D_CIRCLE;
    }

    static constexpr ShapeType2D shapeType = SHAPE2D_CIRCLE;

    std::vector<Vector2D> getVertices(int resolution = 32) const
    {
        std::vector<Vector2D> vertices;
        vertices.reserve(resolution);

        const real_t step = 2.0f * static_cast<real_t>(PI) / resolution;
        for (int i = 0; i < resolution; ++i)
        {
            real_t angle = i * step;
            real_t x = m_center.x + m_radius * std::cos(angle);
            real_t y = m_center.y + m_radius * std::sin(angle);
            vertices.emplace_back(x, y);
        }

        return vertices;
    }

    real_t area() const override
    {
        return PI * m_radius * m_radius;
    }

    real_t circumference() const
    {
        return 2.0f * PI * m_radius;
    }

    real_t perimeter() const override
    {
        return circumference();
    }

    Vector2D centroid() const override
    {
        return m_center;
    }

    Circle2D& translate(const Vector2D &translation) override
    {
        m_center += translation;
        return *this;
    }

    Circle2D& rotate(real_t angle)
    {
        return *this;
    }

    Circle2D& rotate(real_t angle, const Vector2D& point)
    {
        this->m_center.rotateAround(angle, point);
        return *this;
    }

    Circle2D& scale(real_t factor)
    {
        m_radius *= factor;
        return *this;
    }

    Circle2D copy() const
    {
        return *this;
    }

    std::unique_ptr<IFiniteShape2D> clone() const
    {
        return std::make_unique<Circle2D>(*this);        
    }

    bool contains(const Vector2D &point) const override
    {
        return (point - m_center).lengthSquared() <= m_radius * m_radius;
    }

    bool contains(const Circle2D &circle) const
    {
        return (m_center - circle.m_center).lengthSquared() + circle.m_radius * circle.m_radius <= m_radius * m_radius;
    }

    bool intersects(const Circle2D &circle) const
    {
        return (m_center - circle.m_center).lengthSquared() <= (m_radius + circle.m_radius) * (m_radius + circle.m_radius);
    }

    BBox2D boundingBox() const override
    {
        return BBox2D(m_center - Vector2D(m_radius, m_radius), m_center + Vector2D(m_radius, m_radius));
    }
};

} // namespace Math

} // namespace Arns