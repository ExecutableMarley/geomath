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
    float x;
    float y;
    float z;

    Vector3D() : x(0), y(0), z(0) {}

    Vector3D(float x, float y, float z) : x(x), y(y), z(z) {}


    float length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    float lengthSquared() const
    {
        return x * x + y * y + z * z;
    }

    bool isZero() const
    {
        return approximatelyZero(x) && approximatelyZero(y) && approximatelyZero(z);
    }

    Vector3D& normalize()
    {
        const float len = length();
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
        const float len = length();
        if (len != 0)
        {
            return Vector3D(x / len, y / len, z / len);
        }
        return Vector3D();
    }

    Vector3D& resize(float newLength)
    {
        const float len = length();
        if (len != 0)
        {
            x *= newLength / len;
            y *= newLength / len;
            z *= newLength / len;
        }
        return *this;
    }

    Vector3D createResized(float newLength) const
    {
        const float len = length();
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

    float dot(const Vector3D& other) const
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

    Vector3D& rotateAroundX(float angle)
    {
        const float cosAngle = cos(angle);
        const float sinAngle = sin(angle);
        const float newY = y * cosAngle - z * sinAngle;
        z = y * sinAngle + z * cosAngle;
        y = newY;
        return *this;
    }

    Vector3D& rotateAroundY(float angle)
    {
        const float cosAngle = cos(angle);
        const float sinAngle = sin(angle);
        const float newZ = z * cosAngle - x * sinAngle;
        x = z * sinAngle + x * cosAngle;
        z = newZ;
        return *this;
    }

    Vector3D& rotateAroundZ(float angle)
    {
        const float cosAngle = cos(angle);
        const float sinAngle = sin(angle);
        const float newX = x * cosAngle - y * sinAngle;
        y = x * sinAngle + y * cosAngle;
        x = newX;
        return *this;
    }

    Vector3D& rotate(float xAngle = 0.f, float yAngle = 0.f, float zAngle = 0.f)
    {
        if (xAngle != 0.f)
            rotateAroundX(xAngle);
        if (yAngle != 0.f)
            rotateAroundY(yAngle);
        if (zAngle != 0.f)
            rotateAroundZ(zAngle);
        return *this;
    }

    Vector3D& rotateAround(float xAngle, float yAngle, float zAngle, const Vector3D& point)
    {
        return (*this -= point).rotate(xAngle, yAngle, zAngle) += point;
    }

    operator float*()
    {
        return &x;
    }

    operator const float*() const
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

    Vector3D operator*(float scalar) const
    {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    Vector3D operator/(float scalar) const
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

    Vector3D& operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Vector3D& operator/=(float scalar)
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

inline Vector3D operator *(float scalar, const Vector3D& vector)
{
    return vector * scalar;
}

} // namespace Math

} // namespace Arns