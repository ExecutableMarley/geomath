/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "../CommonMath.hpp"

namespace Arns
{

namespace Math
{

struct Vector2D
{
    real_t x;
    real_t y;

    Vector2D() : x(0), y(0) {}

    Vector2D(real_t x, real_t y) : x(x), y(y) {}


    real_t length() const
    {
        return sqrt(x * x + y * y);
    }

    real_t lengthSquared() const
    {
        return x * x + y * y;
    }

    bool isZero() const
    {
        return approximatelyZero(x) && approximatelyZero(y);
    }

    Vector2D& normalize()
    {
        const real_t len = length();
        if (len != 0)
        {
            x /= len;
            y /= len;
        }
        return *this;
    }

    Vector2D createNormalized() const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Vector2D(x / len, y / len);
        }
        return Vector2D();
    }

    Vector2D& resize(real_t newLength)
    {
        const real_t len = length();
        if (len != 0)
        {
            x *= newLength / len;
            y *= newLength / len;
        }
        return *this;
    }

    Vector2D createResized(real_t newLength) const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Vector2D(x * newLength / len, y * newLength / len);
        }
        return Vector2D();
    }

    Vector2D& clamp(const Vector2D& min, const Vector2D& max)
    {
        x = Arns::Math::clamp(x, min.x, max.x);
        y = Arns::Math::clamp(y, min.y, max.y);
        return *this;
    }

    Vector2D createPerpendicular() const
    {
        return Vector2D(-y, x);
    }

    Vector2D createUnitPerpendicular() const
    {
        return createPerpendicular().normalize();
    }

    real_t distance(const Vector2D& other) const
    {
        return (*this - other).length();
    }

    real_t distanceSquared(const Vector2D& other) const
    {
        return (*this - other).lengthSquared();
    }

    real_t dot(const Vector2D& other) const
    {
        return x * other.x + y * other.y;
    }

    real_t cross(const Vector2D& other) const
    {
        return x * other.y - y * other.x;
    }

    bool isNormalized() const
    {
        return approximatelyEqual(length(), 1);
    }

    bool isParallel(const Vector2D& other) const
    {
        return approximatelyZero(cross(other));
    }

    bool isOrthogonal(const Vector2D& other) const
    {
        return approximatelyZero(dot(other));
    }

    bool isPerpendicular(const Vector2D& other) const
    {
        return isOrthogonal(other);
    }

    Vector2D& rotate(real_t degree)
    {
        real_t cosAngle, sinAngle;
        sinCosDeg(degree, sinAngle, cosAngle);

        const real_t tempX = x;
        x = x * cosAngle - y * sinAngle;
        y = tempX * sinAngle + y * cosAngle;

        return *this;
    }

    Vector2D& rotateAround(real_t degree, const Vector2D& point)
    {
        return (*this -= point).rotate(degree) += point;
    }

    operator real_t*()
    {
        return &x;
    }

    operator const real_t*() const
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

    Vector2D operator*(real_t scalar) const
    {
        return Vector2D(x * scalar, y * scalar);
    }

    Vector2D operator/(real_t scalar) const
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

    Vector2D& operator*=(real_t scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vector2D& operator/=(real_t scalar)
    {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    bool operator==(const Vector2D& other) const
    {
        return approximatelyEqual(x, other.x) && approximatelyEqual(y, other.y);
    }

    bool operator!=(const Vector2D& other) const
    {
        return !approximatelyEqual(x, other.x) || !approximatelyEqual(y, other.y);
    }

    Vector2D operator-() const
    {
        return Vector2D(-x, -y);
    }

    // Static functions

    static Vector2D min(const Vector2D& a, const Vector2D& b)
    {
        return Vector2D(a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y);
    }

    template <typename... Args>
    static Vector2D min(const Vector2D& a, const Vector2D& b, Args... args)
    {
        return min(a, min(b, args...));
    }

    static Vector2D max(const Vector2D& a, const Vector2D& b)
    {
        return Vector2D(a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y);
    }

    template <typename... Args>
    static Vector2D max(const Vector2D& a, const Vector2D& b, Args... args)
    {
        return max(a, max(b, args...));
    }
};

inline Vector2D operator *(real_t scalar, const Vector2D& vector)
{
    return vector * scalar;
}

inline real_t orient2D(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool isCCW(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return approximatelyGreaterAbs(orient2D(a, b, c), 0);
}

inline bool isCW(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return approximatelyLessAbs(orient2D(a, b, c), 0);
}

inline bool isColinear(const Vector2D& a, const Vector2D& b, const Vector2D& c)
{
    return approximatelyZero(orient2D(a,b,c));
}

} // namespace Math

} // namespace Arns