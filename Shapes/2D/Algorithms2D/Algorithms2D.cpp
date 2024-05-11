/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms2D.hpp"

namespace Utility
{

namespace Math
{

float distancePointToLine(const Vector2D& point, const Vector2D& lineStart, const Vector2D& lineEnd, Vector2D* closestPoint = nullptr)
{
    return distancePointToLine(point, Line2D(lineStart, lineEnd), closestPoint);

    const Vector2D r = lineEnd - lineStart;
    const Vector2D s = point - lineStart;
    const float t = s.dot(r) / r.dot(r);
    const Vector2D projection = lineStart + r * t;
    if (closestPoint)
        *closestPoint = projection;
    return (point - projection).length();
}

float distancePointToLine(const Vector2D& point, const Line2D& line, Vector2D* closestPoint = nullptr)
{
    return distancePointToLine(point, line.m_start, line.m_end, closestPoint);
}

float distanceLineToLine(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr)
{
    const Vector2D delta1 = s2 - s1;
    const Vector2D delta2 = k2 - k1;
    float lengthSquared1 = delta1.lengthSquared();
    float lengthSquared2 = delta2.lengthSquared();

    if (lengthSquared1 == 0 || lengthSquared2 == 0)
    {
        if (closestPoint1)
            *closestPoint1 = s1;
        if (closestPoint2)
            *closestPoint2 = k1;
        return (s1 - k1).length();
    }

    const Vector2D startDelta = s1 - k1;
    
    const float a = delta1.dot(startDelta);
    const float b = delta2.dot(startDelta);
    const float c = delta1.dot(delta2);

    const float d = lengthSquared1 * lengthSquared2 - c * c;

    float t = d != 0 ? (b * c - a * lengthSquared2) / d : 0.0f; 
    float u = d != 0 ? (lengthSquared1 * b - c * a) / d : 0.0f;

    t = clamp(t, 0.0f, 1.0f);
    u = clamp(u, 0.0f, 1.0f);

    const Vector2D closestPointOnLine1 = s1 + delta1 * t;
    const Vector2D closestPointOnLine2 = k1 + delta2 * u;

    if (closestPoint1)
        *closestPoint1 = closestPointOnLine1;
    if (closestPoint2)
        *closestPoint2 = closestPointOnLine2;
    
    return (closestPointOnLine1 - closestPointOnLine2).length();
}

float distanceLineToLine(const Line2D& line1, const Line2D& line2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr)
{
    return distanceLineToLine(line1.m_start, line1.m_end, line2.m_start, line2.m_end, closestPoint1, closestPoint2);
}

} // namespace Math

} // namespace Utility