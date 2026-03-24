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

struct Vector4D
{
    real_t x;
    real_t y;
    real_t z;
    real_t w;

    Vector4D() : x(0), y(0), z(0), w(0) {}

    Vector4D(real_t x, real_t y, real_t z, real_t w) : x(x), y(y), z(z), w(w) {}

    real_t length() const
    {
        return sqrt(x * x + y * y + z * z + w * w);
    }

    real_t lengthSquared() const
    {
        return x * x + y * y + z * z + w * w;
    }

    Vector4D& normalize()
    {
        const real_t len = length();
        if (len != 0)
        {
            x /= len;
            y /= len;
            z /= len;
            w /= len;
        }
        return *this;
    }

    Vector4D createNormalized() const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Vector4D(x / len, y / len, z / len, w / len);
        }
        return Vector4D();
    }

    Vector4D& resize(real_t newLength)
    {
        const real_t len = length();
        if (len != 0)
        {
            x *= newLength / len;
            y *= newLength / len;
            z *= newLength / len;
            w *= newLength / len;
        }
        return *this;
    }

    Vector4D createResized(real_t newLength) const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Vector4D(x * newLength / len, y * newLength / len, z * newLength / len, w * newLength / len);
        }
        return Vector4D();
    }

    operator real_t*()
    {
        return &x;
    }

    operator const real_t*() const
    {
        return &x;
    }

    Vector4D operator+(const Vector4D &other) const
    {
        return Vector4D(x + other.x, y + other.y, z + other.z, w + other.w);
    }

    Vector4D operator-(const Vector4D &other) const
    {
        return Vector4D(x - other.x, y - other.y, z - other.z, w - other.w);
    }

    Vector4D operator*(real_t scalar) const
    {
        return Vector4D(x * scalar, y * scalar, z * scalar, w * scalar);
    }

    Vector4D operator/(real_t scalar) const
    {
        return Vector4D(x / scalar, y / scalar, z / scalar, w / scalar);
    }

    Vector4D& operator+=(const Vector4D &other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        w += other.w;
        return *this;
    }

    Vector4D& operator-=(const Vector4D &other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        w -= other.w;
        return *this;
    }

    Vector4D& operator*=(real_t scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    Vector4D& operator/=(real_t scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        w /= scalar;
        return *this;
    }

    bool operator==(const Vector4D &other) const
    {
        return approximatelyEqual(x, other.x) && approximatelyEqual(y, other.y) && approximatelyEqual(z, other.z) && approximatelyEqual(w, other.w);
    }

    bool operator!=(const Vector4D &other) const
    {
        return !approximatelyEqual(x, other.x) || !approximatelyEqual(y, other.y) || !approximatelyEqual(z, other.z) || !approximatelyEqual(w, other.w);
    }

    Vector4D operator-() const
    {
        return Vector4D(-x, -y, -z, -w);
    }
};

} // namespace Math

} // namespace Arns