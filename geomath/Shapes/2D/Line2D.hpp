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

class Line2D
{
public:
    Vector2D m_start;
    Vector2D m_end;

    // --- Constructors ---

    Line2D() : m_start(), m_end() {}

    Line2D(Vector2D startPoint, Vector2D endPoint) : m_start(startPoint), m_end(endPoint) {}

    Line2D(Vector2D startPoint, Vector2D direction, real_t length) : m_start(startPoint), m_end(startPoint + direction.createNormalized() * length) {}

    // --- Geometric Properties ---

    real_t length() const
    {
        return (m_start - m_end).length();
    }

    Vector2D origin() const
    {
        return m_start;
    }

    Vector2D direction() const
    {
        return (m_end - m_start).normalize();
    }

    Vector2D deltaVector() const
    {
        return m_end - m_start;
    }

    Vector2D pointAt(real_t t) const
    {
        return m_start + (m_end - m_start) * t;
    }

    Vector2D normal() const
    {
        return Vector2D(m_end.y - m_start.y, m_start.x - m_end.x).normalize();
    }

    real_t closestParameter(const Vector2D &point) const
    {
        Vector2D deltaVector = m_end - m_start;
        const real_t lengthSquared = deltaVector.lengthSquared();
        return approximatelyZero(lengthSquared) ? real_t{0} : (point - m_start).dot(deltaVector) / lengthSquared;
    }

    // --- Transform / modification ---

    Line2D& translate(const Vector2D &translation)
    {
        m_start += translation;
        m_end += translation;
        return *this;
    }

    Line2D& rotate(real_t angle, const Vector2D& point)
    {
        m_start.rotateAround(angle, point);
        m_end.rotateAround(angle, point);
        return *this;
    }

    Line2D& rotate(real_t angle)
    {
        const Vector2D c = (m_start + m_end) * real_t{0.5};
        return rotate(angle, c);
    }

    Line2D& scale(real_t factor, const Vector2D& point)
    {
        m_start = point + (m_start - point) * factor;
        m_end = point + (m_end - point) * factor;
        return *this;
    }

    Line2D& scale(real_t factor)
    {
        const Vector2D c = (m_start + m_end) * real_t{0.5};
        return scale(factor, c);
    }

    Line2D& transform(const Matrix3x3& matrix)
    {
        m_start = matrix.transformPoint(m_start);
        m_end = matrix.transformPoint(m_end);
        return *this;
    }

    // --- Lifecycle / Factory Methods ---

    Line2D copy() const
    {
        return *this;
    }

    // --- Operators ---

    bool operator==(const Line2D& other) const
    {
        return this->m_start == other.m_start && this->m_end == other.m_end;
    }

    bool operator!=(const Line2D& other) const
    {
        return !(*this == other);
    }
};

} // namespace Math

} // namespace Arns