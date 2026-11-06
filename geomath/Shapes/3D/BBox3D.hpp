/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include <algorithm>

namespace Arns
{

namespace geomath
{

class BBox3D
{
public:
    Vector3D m_min;
    Vector3D m_max;

    // --- Constructors ---

    BBox3D() : m_min(), m_max() {}

    BBox3D(const Vector3D &min, const Vector3D &max) : m_min(min), m_max(max) {}

    BBox3D(const Vector3D& point) : m_min(point), m_max(point) {}

    // --- Derived geometry ---

    Vector3D halfSize() const
    {
        return (m_max - m_min) * 0.5f;
    }

    Vector3D size() const
    {
        return m_max - m_min;
    }

    real_t volume() const
    {
        const Vector3D s = size();
        return s.x * s.y * s.z;
    }

    real_t surfaceArea() const
    {
        const Vector3D s = size();
        return 2.0f * (s.x * s.y + s.x * s.z + s.y * s.z);
    }

    Vector3D centroid() const
    {
        return (m_min + m_max) * 0.5f;
    }

    // --- State Queries ---

    // --- Transform / modification ---

    BBox3D& translate(const Vector3D &translation)
    {
        m_min += translation;
        m_max += translation;
        return *this;
    }

    BBox3D& scale(real_t scaleFactor)
    {
        Vector3D center = centroid();
        Vector3D vHalfSize = this->halfSize();
        m_min = center - vHalfSize * scaleFactor;
        m_max = center + vHalfSize * scaleFactor;
        return *this;
    }

    BBox3D& encapsulate(const Vector3D &point)
    {
        m_min = Vector3D::min(m_min, point);
        m_max = Vector3D::max(m_max, point);
        return *this;
    }

    BBox3D& encapsulate(const BBox3D &box)
    {
        m_min = Vector3D::min(m_min, box.m_min);
        m_max = Vector3D::max(m_max, box.m_max);
        return *this;
    }

    /*
    Todo: Test this
    Vector3D normalAt(const Vector3D &point) const
    {
        const Vector3D half = halfSize();
        const Vector3D d = (point - m_min) / half;
        const Vector3D absD = Vector3D(fabs(d.x), fabs(d.y), fabs(d.z));
        const Vector3D sign = Vector3D(d.x < 0.0f ? -1.0f : 1.0f, d.y < 0.0f ? -1.0f : 1.0f, d.z < 0.0f ? -1.0f : 1.0f);
        const Vector3D normal = Vector3D(absD.x == 1.0f ? sign.x : 0.0f, absD.y == 1.0f ? sign.y : 0.0f, absD.z == 1.0f ? sign.z : 0.0f);
        return normal;
    }*/

    // --- Lifecycle / Factory Methods ---

    // --- ---

    real_t minDistanceSquared(const Vector3D& point) const
    {
        const real_t dx = std::max(real_t{0}, std::max(m_min.x - point.x, point.x - m_max.x));
        const real_t dy = std::max(real_t{0}, std::max(m_min.y - point.y, point.y - m_max.y));
        const real_t dz = std::max(real_t{0}, std::max(m_min.z - point.z, point.z - m_max.z));
        return dx * dx + dy * dy + dz * dz;
    }

    real_t minDistance(const Vector3D& point) const
    {
        return sqrt(minDistanceSquared(point));
    }

    // --- ---

    bool contains(const Vector3D &point) const
    {
        return point.x >= m_min.x && point.x <= m_max.x && point.y >= m_min.y && point.y <= m_max.y && point.z >= m_min.z && point.z <= m_max.z;
    }

    bool contains(const BBox3D &box) const
    {
        return contains(box.m_min) && contains(box.m_max);
    }

    bool intersects(const BBox3D &box) const
    {
        return m_min.x <= box.m_max.x && m_max.x >= box.m_min.x && m_min.y <= box.m_max.y && m_max.y >= box.m_min.y && m_min.z <= box.m_max.z && m_max.z >= box.m_min.z;
    }

    // --- Comparison Operators ---

    // --- Stream Output ---

    // --- Static factory Functions ---

    static BBox3D fromMinSize(const Vector3D &min, const Vector3D &size)
    {
        return BBox3D(min, min + size);
    }

    static BBox3D fromPoint(const Vector3D &point)
    {
        return BBox3D(point, point);
    }

    static BBox3D fromTwoPoints(const Vector3D& a, const Vector3D& b)
    {
        return BBox3D(Vector3D::min(a, b), Vector3D::max(a, b));
    }

    template <typename... Points> requires (std::same_as<std::remove_cvref_t<Points>, Vector3D> && ...)
    static BBox3D fromPoints(const Vector3D& first, const Points&... points)
    {
        BBox3D result(first);
        (result.encapsulate(points), ...);
        return result;
    }

    template <typename Iter>
    static BBox3D fromRange(Iter begin, Iter end)
    {
        if (begin == end) return {};

        Vector3D min = *begin;
        Vector3D max = *begin;
        for (auto it = std::next(begin); it != end; it++)
        {
            min = Vector3D::min(min, *it);
            max = Vector3D::max(max, *it);
        }
        return BBox3D(min, max);
    }

    static BBox3D merge(const BBox3D &box1, const BBox3D &box2)
    {
        return BBox3D(Vector3D::min(box1.m_min, box2.m_min), Vector3D::max(box1.m_max, box2.m_max));
    }

    static BBox3D FromPointAndRadius(const Vector3D& center, real_t radius)
    {
        return BBox3D(
            { center.x - radius, center.y - radius, center.z - radius },
            { center.x + radius, center.y + radius, center.z + radius }
        );
    }

    static BBox3D sweptBounds(const BBox3D &bounds, const Vector3D &displacement)
    {
        return {
            Vector3D::min(bounds.m_min, bounds.m_min + displacement),
            Vector3D::max(bounds.m_max, bounds.m_max + displacement)};
    }
};

} // namespace Math

} // namespace Arns
