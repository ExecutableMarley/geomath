/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "BBox3D.hpp"
#include "Algorithms3D/Algorithms3D.hpp"

namespace Utility
{

namespace Math
{

struct Capsule
{
    Vector3D m_startPoint;
    Vector3D m_endPoint;
    float m_radius;

    Capsule() : m_startPoint(), m_endPoint(), m_radius(0) {}
    
    Capsule(const Vector3D &startPoint, const Vector3D &endPoint, float radius) : m_startPoint(startPoint), m_endPoint(endPoint), m_radius(radius) {}

    float volume() const
    {
        return PI * m_radius * m_radius * (4.0f / 3.0f * m_radius + (m_endPoint - m_startPoint).length());
    }

    float surfaceArea() const
    {
        return 2.0f * PI * m_radius * (2.0f * m_radius + (m_endPoint - m_startPoint).length());
    }

    Vector3D centroid() const
    {
        return (m_startPoint + m_endPoint) * 0.5f;
    }

    bool contains(const Vector3D &point) const
    {
        return distancePointToLine(m_startPoint, m_endPoint, point) < m_radius;
    }

    BBox3D boundingBox() const
    {
        const Vector3D min = Vector3D::min(m_startPoint, m_endPoint) - Vector3D(m_radius, m_radius, m_radius);
        const Vector3D max = Vector3D::max(m_startPoint, m_endPoint) + Vector3D(m_radius, m_radius, m_radius);
        return BBox3D(min, max);
    }
};


} // namespace Math

} // namespace Utility