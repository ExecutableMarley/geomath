/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <algorithm>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"

namespace Arns
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

    BBox2D(const Vector2D &point, real_t width = 0.f, real_t height = 0.f) : m_min(point), m_max(point.x + width, point.y + height) {}

    // --- Derived geometry ---

    real_t width() const
    {
        return m_max.x - m_min.x;
    }

    real_t height() const
    {
        return m_max.y - m_min.y;
    }

    real_t area() const
    {
        return width() * height();
    }

    real_t perimeter() const
    {
        return 2.0f * (width() + height());
    }

    Vector2D centroid() const
    {
        return (m_min + m_max) * 0.5f;
    }

    Vector2D bottomLeft() const
    {
        return { m_min.x, m_min.y };
    }

    Vector2D bottomRight() const
    {
        return { m_max.x, m_min.y };
    }

    Vector2D topRight() const
    {
        return { m_max.x, m_max.y };
    }

    Vector2D topLeft() const
    {
        return { m_min.x, m_max.y };
    }

    // --- Transform / modification ---

    BBox2D& translate(const Vector2D &translation)
    {
        m_min += translation;
        m_max += translation;
        return *this;
    }

    BBox2D& scale(real_t factor)
    {
        Vector2D center = centroid();
        Vector2D halfSize = (m_max - m_min) * 0.5f;
        m_min = center - halfSize * factor;
        m_max = center + halfSize * factor;
        return *this;
    }

    BBox2D& pad(real_t amount)
    {
        m_min = m_min - Vector2D(amount, amount);
        m_max = m_max + Vector2D(amount, amount);
        return *this;
    }

    BBox2D& encapsulate(const Vector2D &point)
    {
        m_min = Vector2D::min(m_min, point);
        m_max = Vector2D::max(m_max, point);
        return *this;
    }

    BBox2D& encapsulate(const BBox2D &bbox)
    {
        m_min = Vector2D::min(m_min, bbox.m_min);
        m_max = Vector2D::max(m_max, bbox.m_max);
        return *this;
    }

    // --- ---

    Vector2D closestPoint(const Vector2D& p) const
    {
        return Vector2D(
            std::max(m_min.x, std::min(p.x, m_max.x)),
            std::max(m_min.y, std::min(p.y, m_max.y))
        );
    }

    real_t minDistanceSquared(const Vector2D& point) const
    {
        real_t dx = std::max(0.0f, std::max(m_min.x - point.x, point.x - m_max.x));
        real_t dy = std::max(0.0f, std::max(m_min.y - point.y, point.y - m_max.y));
        return dx * dx + dy * dy;
    }

    real_t maxDistanceSquared(const Vector2D& point) const
    {
        real_t dx1 = point.x - m_min.x;
        real_t dx2 = point.x - m_max.x;
        real_t dy1 = point.y - m_min.y;
        real_t dy2 = point.y - m_max.y;
        return std::max(dx1*dx1, dx2*dx2) + std::max(dy1*dy1, dy2*dy2);
    }

    // --- ---

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

    // --- Static factory Functions ---

    template <typename Iter>
    static BBox2D fromRange(Iter begin, Iter end)
    {
        if (begin == end) return {};

        Vector2D min = *begin;
        Vector2D max = *begin;
        for (auto it = std::next(begin); it != end; it++)
        {
            min = Vector2D::min(min, *it);
            max = Vector2D::max(max, *it);
        }
        return BBox2D(min, max);
    }

    static BBox2D fromTwoPoints(const Vector2D& a, const Vector2D& b)
    {
        return BBox2D(Vector2D::min(a, b), Vector2D::max(a, b));
    }

    static BBox2D merge(const BBox2D &box1, const BBox2D &box2)
    {
        return BBox2D(Vector2D::min(box1.m_min, box2.m_min), Vector2D::max(box1.m_max, box2.m_max));
    }
};

} // namespace Math

} // namespace Arns