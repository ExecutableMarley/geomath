/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector2D.hpp"

namespace Utility
{

namespace Math
{

class BBox2D
{
public:
    Vector2D m_min;
    Vector2D m_max;

    BBox2D() : m_min(), m_max() {}
    
    BBox2D(const Vector2D &min, const Vector2D &max) : m_min(min), m_max(max) {}

    BBox2D(const Vector2D &point, float width = 0.f, float height = 0.f) : m_min(point), m_max(point.x + width, point.y + height) {}

    float width() const
    {
        return m_max.x - m_min.x;
    }

    float height() const
    {
        return m_max.y - m_min.y;
    }

    float area() const
    {
        return width() * height();
    }

    float perimeter() const
    {
        return 2.0f * (width() + height());
    }

    Vector2D centroid() const
    {
        return (m_min + m_max) * 0.5f;
    }

    BBox2D& translate(const Vector2D &translation)
    {
        m_min += translation;
        m_max += translation;
        return *this;
    }

    BBox2D& scale(float factor)
    {
        Vector2D center = centroid();
        Vector2D halfSize = (m_max - m_min) * 0.5f;
        m_min = center - halfSize * factor;
        m_max = center + halfSize * factor;
        return *this;
    }

    BBox2D& encapsulate(const Vector2D &point)
    {
        m_min = Vector2D::min(m_min, point);
        m_max = Vector2D::max(m_max, point);
        return *this;
    }

    bool contains(const Vector2D &point) const
    {
        return point.x >= m_min.x && point.x <= m_max.x && point.y >= m_min.y && point.y <= m_max.y;
    }

    bool contains(const BBox2D &rectangle) const
    {
        return contains(rectangle.m_min) && contains(rectangle.m_max);
    }

    bool intersects(const BBox2D &rectangle) const
    {
        return m_min.x <= rectangle.m_max.x && m_max.x >= rectangle.m_min.x && m_min.y <= rectangle.m_max.y && m_max.y >= rectangle.m_min.y;
    }

    // Static functions

    static BBox2D merge(const BBox2D &box1, const BBox2D &box2)
    {
        return BBox2D(Vector2D::min(box1.m_min, box2.m_min), Vector2D::max(box1.m_max, box2.m_max));
    }
};

} // namespace Math

} // namespace Utility