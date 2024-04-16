#pragma once

#include <math.h>

namespace Utility
{

namespace Math
{

struct Vector2D
{
    float x;
    float y;

    Vector2D() : x(0), y(0) {}

    Vector2D(float x, float y) : x(x), y(y) {}


    float length() const
    {
        return sqrt(x * x + y * y);
    }

    float lengthSquared() const
    {
        return x * x + y * y;
    }

    Vector2D& normalize()
    {
        const float len = length();
        if (len != 0)
        {
            x /= len;
            y /= len;
        }
        return *this;
    }

    Vector2D createNormalized() const
    {
        const float len = length();
        if (len != 0)
        {
            return Vector2D(x / len, y / len);
        }
        return Vector2D();
    }

    Vector2D& resize(float newLength)
    {
        const float len = length();
        if (len != 0)
        {
            x *= newLength / len;
            y *= newLength / len;
        }
        return *this;
    }

    Vector2D createResized(float newLength) const
    {
        const float len = length();
        if (len != 0)
        {
            return Vector2D(x * newLength / len, y * newLength / len);
        }
        return Vector2D();
    }

    float dot(const Vector2D& other) const
    {
        return x * other.x + y * other.y;
    }

    float cross(const Vector2D& other) const
    {
        return x * other.y - y * other.x;
    }

    operator float*()
    {
        return &x;
    }

    operator const float*() const
    {
        return &x;
    }

    Vector2D operator+(const Vector2D& other) const
    {
        return Vector2D(x + other.x, y + other.y);
    }

    Vector2D operator-(const Vector2D& other) const
    {
        return Vector2D(x - other.x, y - other.y);
    }

    Vector2D operator*(float scalar) const
    {
        return Vector2D(x * scalar, y * scalar);
    }

    Vector2D operator/(float scalar) const
    {
        return Vector2D(x / scalar, y / scalar);
    }

    Vector2D& operator+=(const Vector2D& other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2D& operator-=(const Vector2D& other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2D& operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vector2D& operator/=(float scalar)
    {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    bool operator==(const Vector2D& other) const
    {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Vector2D& other) const
    {
        return x != other.x || y != other.y;
    }

    Vector2D operator-(const Vector2D& other) const
    {
        return Vector2D(x - other.x, y - other.y);
    }
};

} // namespace Math

} // namespace Utility