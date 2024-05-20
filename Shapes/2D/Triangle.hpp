/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Utility
{

namespace Math
{

class Triangle : public IShape2D
{
public:
    Vector2D m_a;
    Vector2D m_b;
    Vector2D m_c;

    Triangle() : m_a(), m_b(), m_c() {}

    Triangle(const Vector2D &a, const Vector2D &b, const Vector2D &c) : m_a(a), m_b(b), m_c(c) {}

    ShapeType2D type() const override
    {
        return SHAPE2D_TRIANGLE;
    }

    float area() const override
    {
        return 0.5f * fabs((m_a.x - m_c.x) * (m_b.y - m_c.y) - (m_b.x - m_c.x) * (m_a.y - m_c.y));
    }

    float perimeter() const override
    {
        return (m_a - m_b).length() + (m_b - m_c).length() + (m_c - m_a).length();
    }

    Vector2D centroid() const override
    {
        return (m_a + m_b + m_c) / 3.0f;
    }

    Triangle& translate(const Vector2D &translation) override
    {
        m_a += translation;
        m_b += translation;
        m_c += translation;
        return *this;
    }

    Triangle& rotate(float angle)
    {
        const Vector2D centroid = this->centroid();
        m_a.rotateAround(angle, centroid);
        m_b.rotateAround(angle, centroid);
        m_c.rotateAround(angle, centroid);
        return *this;
    }

    bool contains(const Vector2D &point) const override
    {
        const float areaABC = area();
        const float areaPBC = 0.5f * fabs((m_b.x - point.x) * (m_c.y - point.y) - (m_c.x - point.x) * (m_b.y - point.y));
        const float areaPCA = 0.5f * fabs((m_c.x - point.x) * (m_a.y - point.y) - (m_a.x - point.x) * (m_c.y - point.y));
        const float areaPAB = 0.5f * fabs((m_a.x - point.x) * (m_b.y - point.y) - (m_b.x - point.x) * (m_a.y - point.y));
        return approximatelyEqual(areaABC, areaPBC + areaPCA + areaPAB);
    }

    bool contains(const Triangle &triangle) const
    {
        return contains(triangle.m_a) && contains(triangle.m_b) && contains(triangle.m_c);
    }

    BBox2D boundingBox() const override
    {
        return BBox2D(Vector2D(fminf(fminf(m_a.x, m_b.x), m_c.x), fminf(fminf(m_a.y, m_b.y), m_c.y)),
                      Vector2D(fmaxf(fmaxf(m_a.x, m_b.x), m_c.x), fmaxf(fmaxf(m_a.y, m_b.y), m_c.y)));
    }
};

} // namespace Math

} // namespace Utility