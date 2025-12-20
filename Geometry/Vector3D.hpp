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

struct Vector3D
{
    real_t x;
    real_t y;
    real_t z;

    Vector3D() : x(0), y(0), z(0) {}

    Vector3D(real_t x, real_t y, real_t z) : x(x), y(y), z(z) {}


    real_t length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    real_t lengthSquared() const
    {
        return x * x + y * y + z * z;
    }

    bool isZero() const
    {
        return approximatelyZero(x) && approximatelyZero(y) && approximatelyZero(z);
    }

    Vector3D& normalize()
    {
        const real_t len = length();
        if (len != 0)
        {
            x /= len;
            y /= len;
            z /= len;
        }
        return *this;
    }

    Vector3D createNormalized() const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Vector3D(x / len, y / len, z / len);
        }
        return Vector3D();
    }

    Vector3D& resize(real_t newLength)
    {
        const real_t len = length();
        if (len != 0)
        {
            x *= newLength / len;
            y *= newLength / len;
            z *= newLength / len;
        }
        return *this;
    }

    Vector3D createResized(real_t newLength) const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Vector3D(x * newLength / len, y * newLength / len, z * newLength / len);
        }
        return Vector3D();
    }

    Vector3D& clamp(const Vector3D& min, const Vector3D& max)
    {
        x = Arns::Math::clamp(x, min.x, max.x);
        y = Arns::Math::clamp(y, min.y, max.y);
        z = Arns::Math::clamp(z, min.z, max.z);
        return *this;
    }

    real_t dot(const Vector3D& other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3D cross(const Vector3D& other) const
    {
        return Vector3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x);
    }

    bool isParallel(const Vector3D& other) const
    {
        return approximatelyZero(cross(other).length());
    }

    bool isOrthogonal(const Vector3D& other) const
    {
        return approximatelyZero(dot(other));
    }

    Vector3D& rotateAroundX(real_t angle)
    {
        const real_t cosAngle = cos(angle);
        const real_t sinAngle = sin(angle);
        const real_t newY = y * cosAngle - z * sinAngle;
        z = y * sinAngle + z * cosAngle;
        y = newY;
        return *this;
    }

    Vector3D& rotateAroundY(real_t angle)
    {
        const real_t cosAngle = cos(angle);
        const real_t sinAngle = sin(angle);
        const real_t newZ = z * cosAngle - x * sinAngle;
        x = z * sinAngle + x * cosAngle;
        z = newZ;
        return *this;
    }

    Vector3D& rotateAroundZ(real_t angle)
    {
        const real_t cosAngle = cos(angle);
        const real_t sinAngle = sin(angle);
        const real_t newX = x * cosAngle - y * sinAngle;
        y = x * sinAngle + y * cosAngle;
        x = newX;
        return *this;
    }

    Vector3D& rotate(real_t xAngle = 0.f, real_t yAngle = 0.f, real_t zAngle = 0.f)
    {
        if (xAngle != 0.f)
            rotateAroundX(xAngle);
        if (yAngle != 0.f)
            rotateAroundY(yAngle);
        if (zAngle != 0.f)
            rotateAroundZ(zAngle);
        return *this;
    }

    Vector3D& rotateAround(real_t xAngle, real_t yAngle, real_t zAngle, const Vector3D& point)
    {
        return (*this -= point).rotate(xAngle, yAngle, zAngle) += point;
    }

    operator real_t*()
    {
        return &x;
    }

    operator const real_t*() const
    {
        return &x;
    }

    Vector3D operator+(const Vector3D& other) const
    {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    Vector3D operator-(const Vector3D& other) const
    {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    Vector3D operator*(real_t scalar) const
    {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    Vector3D operator/(real_t scalar) const
    {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }

    Vector3D& operator+=(const Vector3D& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    Vector3D& operator-=(const Vector3D& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    Vector3D& operator*=(real_t scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Vector3D& operator/=(real_t scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    bool operator==(const Vector3D& other) const
    {
        return approximatelyEqual(x, other.x) && approximatelyEqual(y, other.y) && approximatelyEqual(z, other.z);
    }

    bool operator!=(const Vector3D& other) const
    {
        return !approximatelyEqual(x, other.x) || !approximatelyEqual(y, other.y) || !approximatelyEqual(z, other.z);
    }

    Vector3D operator-() const
    {
        return Vector3D(-x, -y, -z);
    }

    // Static functions

    static Vector3D min(const Vector3D& a, const Vector3D& b)
    {
        return Vector3D(
            a.x < b.x ? a.x : b.x,
            a.y < b.y ? a.y : b.y,
            a.z < b.z ? a.z : b.z);
    }

    template <typename... Args>
    static Vector3D min(const Vector3D& a, const Vector3D& b, Args... args)
    {
        return min(a, min(b, args...));
    }

    static Vector3D max(const Vector3D& a, const Vector3D& b)
    {
        return Vector3D(
            a.x > b.x ? a.x : b.x,
            a.y > b.y ? a.y : b.y,
            a.z > b.z ? a.z : b.z);
    }

    template <typename... Args>
    static Vector3D max(const Vector3D& a, const Vector3D& b, Args... args)
    {
        return max(a, max(b, args...));
    }
};

inline Vector3D operator *(real_t scalar, const Vector3D& vector)
{
    return vector * scalar;
}

} // namespace Math

} // namespace Arns