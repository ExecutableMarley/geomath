/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <stdexcept>
#include <vector>

#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Utility
{

namespace Math
{

class Rectangle : public IFiniteShape2D
{ 
public:
    Vector2D m_a;
    Vector2D m_b;
    Vector2D m_c;
    Vector2D m_d;

    Rectangle() : m_a(), m_b(), m_c(), m_d() {}

    Rectangle(const Vector2D &a, const Vector2D &b, const Vector2D &c, const Vector2D &d) : m_a(a), m_b(b), m_c(c), m_d(d) {}

    Rectangle(const Vector2D &pos, float width, float height) : m_a(pos), m_b(pos + Vector2D(width, 0)), m_c(pos + Vector2D(width, height)), m_d(pos + Vector2D(0, height)) {}

    ShapeType2D type() const override
    {
        return SHAPE2D_RECTANGLE;
    }

    size_t vertexCount() const
    {
        return 4;
    }

    const Vector2D& vertexAt(size_t index) const
    {
        if (index > 3)
            throw std::out_of_range("Index out of range");

        return (&m_a)[index];
    }

    Vector2D& vertexAt(size_t index)
    {
        if (index > 3)
            throw std::out_of_range("Index out of range");
        return (&m_a)[index];
    }

    std::vector<Vector2D> getVertices() const
    {
        return {m_a, m_b, m_c, m_d};
    }

    float width() const
    {
        return (m_b - m_a).length();
    }

    float height() const
    {
        return (m_d - m_a).length();
    }

    float area() const override
    {
        return width() * height();
    }

    float perimeter() const override
    {
        return 2.0f * (width() + height());
    }

    Vector2D centroid() const override
    {
        return (m_a + m_b + m_c + m_d) / 4.0f;
    }

    Rectangle& translate(const Vector2D &translation) override
    {
        m_a += translation;
        m_b += translation;
        m_c += translation;
        m_d += translation;
        return *this;
    }

    Rectangle& rotate(float angle)
    {
        const Vector2D centroid = this->centroid();
        m_a.rotateAround(angle, centroid);
        m_b.rotateAround(angle, centroid);
        m_c.rotateAround(angle, centroid);
        m_d.rotateAround(angle, centroid);
        return *this;
    }

    Rectangle& rotate(float angle, const Vector2D& point)
    {
        m_a.rotateAround(angle, point);
        m_b.rotateAround(angle, point);
        m_c.rotateAround(angle, point);
        m_d.rotateAround(angle, point);
        return *this;
    }

    Rectangle& scale(float factor)
    {
        const Vector2D centroid = this->centroid();
        m_a = centroid + (m_a - centroid) * factor;
        m_b = centroid + (m_b - centroid) * factor;
        m_c = centroid + (m_c - centroid) * factor;
        m_d = centroid + (m_d - centroid) * factor;
        return *this;
    }

    bool contains(const Vector2D &point) const override
    {
        const Vector2D ab = m_b - m_a;
        const Vector2D ap = point - m_a;
        const Vector2D bc = m_c - m_b;
        const Vector2D bp = point - m_b;
        const Vector2D cd = m_d - m_c;
        const Vector2D cp = point - m_c;
        const Vector2D da = m_a - m_d;
        const Vector2D dp = point - m_d;

        return ab.cross(ap) >= 0 && bc.cross(bp) >= 0 && cd.cross(cp) >= 0 && da.cross(dp) >= 0;
    }

    bool contains(const Rectangle &rectangle) const
    {
        return contains(rectangle.m_a) && contains(rectangle.m_b) && contains(rectangle.m_c) && contains(rectangle.m_d);
    }

    BBox2D boundingBox() const override
    {
        return BBox2D(Vector2D::min(m_a, m_b, m_c, m_d), Vector2D::max(m_a, m_b, m_c, m_d));
    }

    Vector2D& operator[](size_t index)
    {
        return vertexAt(index);
    }

    const Vector2D& operator[](size_t index) const
    {
        return vertexAt(index);
    }
};

} // namespace Math

} // namespace Utility