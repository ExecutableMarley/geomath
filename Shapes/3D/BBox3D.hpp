/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector3D.hpp"

namespace Arns
{

namespace Math
{

class BBox3D
{
public:
    Vector3D m_min;
    Vector3D m_max;

    BBox3D() : m_min(), m_max() {}

    BBox3D(const Vector3D &min, const Vector3D &max) : m_min(min), m_max(max) {}

    BBox3D(const Vector3D &center, float halfSize) : m_min(center - Vector3D(halfSize, halfSize, halfSize)), m_max(center + Vector3D(halfSize, halfSize, halfSize)) {}

    Vector3D halfSize() const
    {
        return (m_max - m_min) * 0.5f;
    }

    Vector3D size() const
    {
        return m_max - m_min;
    }

    float volume() const
    {
        const Vector3D s = size();
        return s.x * s.y * s.z;
    }

    float surfaceArea() const
    {
        const Vector3D s = size();
        return 2.0f * (s.x * s.y + s.x * s.z + s.y * s.z);
    }

    Vector3D centroid() const
    {
        return (m_min + m_max) * 0.5f;
    }

    BBox3D& translate(const Vector3D &translation)
    {
        m_min += translation;
        m_max += translation;
        return *this;
    }

    BBox3D& scale(float scaleFactor)
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

    // Static functions

    static BBox3D merge(const BBox3D &box1, const BBox3D &box2)
    {
        return BBox3D(Vector3D::min(box1.m_min, box2.m_min), Vector3D::max(box1.m_max, box2.m_max));
    }
};

} // namespace Math

} // namespace Arns
