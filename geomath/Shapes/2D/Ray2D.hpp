/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"

namespace Arns
{

namespace Math
{

class Ray2D
{
public:

    // --- Constructors ---

    Vector2D m_origin;
    Vector2D m_direction;

    Ray2D() : m_origin(), m_direction() {}

    Ray2D(const Vector2D &origin, const Vector2D &direction) : m_origin(origin), m_direction(direction) {}

    // --- Geometric Properties ---

    Vector2D origin() const
    {
        return m_origin;
    }

    Vector2D direction() const
    {
        return m_direction;
    }

    Vector2D pointAt(real_t t) const
    {
        return m_origin + m_direction * t;
    }

    real_t closestParameter(const Vector2D &point) const
    {
        return (point - m_origin).dot(m_direction) / m_direction.lengthSquared();
    }

    // --- Transform / modification ---

    Ray2D& translate(const Vector2D &translation)
    {
        m_origin += translation;
        return *this;
    }

    Ray2D& rotate(real_t angle, const Vector2D& point)
    {
        m_origin.rotateAround(angle, point);
        m_direction.rotateAround(angle, point);
        return *this;
    }

    Ray2D& rotate(real_t angle)
    {
        const Vector2D c = pointAt(real_t{0.5});
        return rotate(angle, c);
    }

    Ray2D& scale(real_t factor, const Vector2D& point)
    {
        m_origin = point + (m_origin - point) * factor;
        return *this;
    }

    Ray2D& scale(real_t factor)
    {
        const Vector2D c = pointAt(real_t{0.5});
        return scale(factor, c);
    }

    Ray2D& transform(const Matrix3x3& matrix)
    {
        m_origin = matrix.transformPoint(m_origin);
        m_direction = matrix.transformVector(m_direction);
        return *this;
    }

    // --- Lifecycle / Factory Methods ---

    Ray2D copy() const
    {
        return *this;
    }

    // --- Comparison Operators ---

    bool operator==(const Ray2D& other) const
    {
        return this->m_origin == other.m_origin && this->m_direction == other.m_direction;
    }

    bool operator!=(const Ray2D& other) const
    {
        return !(*this == other);
    }
};

} // namespace Math

} // namespace Arns