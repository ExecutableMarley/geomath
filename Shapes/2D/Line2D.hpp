#pragma once

#include <math.h>

#include "Vector2D.hpp"

namespace Utility
{

namespace Math
{

class Line2D
{
public:
    Vector2D m_start;
    Vector2D m_end;

    Line2D() : m_start(), m_end() {}

    Line2D(Vector2D startPoint, Vector2D endPoint) : m_start(startPoint), m_end(endPoint) {}

    Line2D(Vector2D startPoint, Vector2D direction, float length) : m_start(startPoint), m_end(startPoint + direction.createNormalized() * length) {}

    float length() const
    {
        return (m_start - m_end).length();
    }

    Vector2D direction() const
    {
        return (m_end - m_start).normalize();
    }

    Vector2D deltaVector() const
    {
        return m_end - m_start;
    }
};

} // namespace Math

} // namespace Utility