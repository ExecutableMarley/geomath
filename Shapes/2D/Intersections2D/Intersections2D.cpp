#include "Intersections2D.hpp"
#include "CommonMath.hpp"

namespace Utility
{

namespace Math
{

bool intersects(const Line2D& line1, const Line2D& line2, Vector2D* intersection)
{
    const Vector2D line1Delta = line1.deltaVector();
    const Vector2D line2Delta = line2.deltaVector();
    float denominator = line1Delta.cross(line2Delta);

    if (approximatelyZero(denominator))
        return false;

    const Vector2D lineToLine = line1.m_start - line2.m_start;

    float factor1 = lineToLine.cross(line2Delta) / denominator;
    float factor2 = lineToLine.cross(line1Delta) / denominator;

    if (factor1 >= 0.0f && factor1 <= 1.0f && factor2 >= 0.0f && factor2 <= 1.0f)
    {
        if (intersection)
            *intersection = line1.m_start + factor1 * line1Delta;
        return true;
    }
    return false;
}

bool intersects(const Line2D& line, const BBox2D& rectangle, Vector2D* intersection)
{
    Line2D top(rectangle.m_min, Vector2D(rectangle.m_max.x, rectangle.m_min.y));
    Line2D right(Vector2D(rectangle.m_max.x, rectangle.m_min.y), rectangle.m_max);
    Line2D bottom(rectangle.m_max, Vector2D(rectangle.m_min.x, rectangle.m_max.y));
    Line2D left(Vector2D(rectangle.m_min.x, rectangle.m_max.y), rectangle.m_min);

    if (intersects(line, top, intersection))
        return true;

    if (intersects(line, right, intersection))
        return true;

    if (intersects(line, bottom, intersection))
        return true;

    if (intersects(line, left, intersection))
        return true;

    return false;
}

bool intersects(const Line2D& line, const Circle& circle, Vector2D* intersection)
{
    const Vector2D lineDelta = line.deltaVector();
    const Vector2D f = line.m_start - circle.m_center;

    //float a = lineDelta.dot(lineDelta);
    const float a = lineDelta.lengthSquared();
    const float b = 2.0f * f.dot(lineDelta);
    //float c = f.dot(f) - circle.m_radius * circle.m_radius;
    const float c = f.lengthSquared() - circle.m_radius * circle.m_radius;

    float discriminant = b * b - 4.0f * a * c;

    if (discriminant < 0.0f)
        return false;

    discriminant = sqrt(discriminant);

    float t1 = (-b - discriminant) / (2.0f * a);
    float t2 = (-b + discriminant) / (2.0f * a);

    if (t1 >= 0.0f && t1 <= 1.0f)
    {
        if (intersection)
            *intersection = line.m_start + t1 * lineDelta;
        return true;
    }
    if (t2 >= 0.0f && t2 <= 1.0f)
    {
        if (intersection)
            *intersection = line.m_start + t2 * lineDelta;
        return true;
    }
    return false;
}

} // namespace Math

} // namespace Utility