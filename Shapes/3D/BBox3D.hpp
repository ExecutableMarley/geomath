#pragma once

#include <math.h>

#include "Vector3D.hpp"

namespace Utility
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

    //BBox3D(const Vector3D &center, const Vector3D &halfSize) : m_min(center - halfSize), m_max(center + halfSize) {}

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
};

} // namespace Math

}
