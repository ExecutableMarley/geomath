/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "../CommonMath.hpp"
#include <math.h>
#include <ostream>

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

    // --- Constructors ---

    constexpr Vector4D() : x(0), y(0), z(0), w(0) {}

    constexpr Vector4D(real_t x, real_t y, real_t z, real_t w) : x(x), y(y), z(z), w(w) {}

    // --- Geometric Properties ---

    real_t length() const
    {
        return sqrt(x * x + y * y + z * z + w * w);
    }

    real_t lengthSquared() const
    {
        return x * x + y * y + z * z + w * w;
    }

    // --- State Queries ---

    bool isZero() const
    {
        return approximatelyZero(x) && approximatelyZero(y) && approximatelyZero(z) && approximatelyZero(w);
    }

    bool isNormalized() const
    {
        return approximatelyEqual(length(), 1);
    }

    // --- Transform / Modification ---

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

    Vector4D& clamp(const Vector4D& min, const Vector4D& max)
    {
        x = Arns::Math::clamp(x, min.x, max.x);
        y = Arns::Math::clamp(y, min.y, max.y);
        z = Arns::Math::clamp(z, min.z, max.z);
        w = Arns::Math::clamp(w, min.w, max.w);
        return *this;
    }

    // --- Derived Vectors ---

    Vector4D copy() const
    {
        return *this;
    }

    // --- Conversion ---

    operator real_t*()
    {
        return &x;
    }

    operator const real_t*() const
    {
        return &x;
    }

    // --- Arithmetic Operators ---

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

    Vector4D operator-() const
    {
        return Vector4D(-x, -y, -z, -w);
    }

    // --- Comparison Operators ---

    bool operator==(const Vector4D &other) const
    {
        return approximatelyEqual(x, other.x) && approximatelyEqual(y, other.y) && approximatelyEqual(z, other.z) && approximatelyEqual(w, other.w);
    }

    bool operator!=(const Vector4D &other) const
    {
        return !approximatelyEqual(x, other.x) || !approximatelyEqual(y, other.y) || !approximatelyEqual(z, other.z) || !approximatelyEqual(w, other.w);
    }

    // --- Constants ---

    static constexpr Vector4D zero()  { return Vector4D(0, 0, 0, 0); }
    static constexpr Vector4D unitX() { return Vector4D(1, 0, 0, 0); }
    static constexpr Vector4D unitY() { return Vector4D(0, 1, 0, 0); }
    static constexpr Vector4D unitZ() { return Vector4D(0, 0, 1, 0); }
    static constexpr Vector4D unitW() { return Vector4D(0, 0, 0, 1); }

    // --- Utility Functions ---

    static Vector4D min(const Vector4D& a, const Vector4D& b)
    {
        return Vector4D(
            a.x < b.x ? a.x : b.x,
            a.y < b.y ? a.y : b.y,
            a.z < b.z ? a.z : b.z,
            a.w < b.w ? a.w : b.w);
    }

    template <typename... Args>
    static Vector4D min(const Vector4D& a, const Vector4D& b, Args... args)
    {
        return min(a, min(b, args...));
    }

    static Vector4D max(const Vector4D& a, const Vector4D& b)
    {
        return Vector4D(
            a.x > b.x ? a.x : b.x,
            a.y > b.y ? a.y : b.y,
            a.z > b.z ? a.z : b.z,
            a.w > b.w ? a.w : b.w);
    }

    template <typename... Args>
    static Vector4D max(const Vector4D& a, const Vector4D& b, Args... args)
    {
        return max(a, max(b, args...));
    }

    // --- Stream Output ---

    friend std::ostream& operator<<(std::ostream& os, const Vector4D& vec) {
        return os << "{" << vec.x << ", " << vec.y << ", " << vec.z << ", " << vec.w << "}";
    }
};

inline Vector4D operator *(real_t scalar, const Vector4D& vector)
{
    return vector * scalar;
}

} // namespace Math

} // namespace Arns