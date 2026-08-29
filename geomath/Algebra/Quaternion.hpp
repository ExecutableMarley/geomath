/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "Matrices/Matrix3x3.hpp"
#include "Matrices/Matrix4x4.hpp"

namespace Arns
{

namespace Math
{

struct Quaternion
{
    real_t x;
    real_t y;
    real_t z;
    real_t w;

    Quaternion() : x(0), y(0), z(0), w(1) {}

    Quaternion(real_t x, real_t y, real_t z, real_t w) : x(x), y(y), z(z), w(w) {}

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

    bool isNormalized() const
    {
        return approximatelyEqual(length(), 1);
    }

    // --- Transform / modification ---

    Quaternion& normalize()
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

    Quaternion& resize(real_t newLength)
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

    // --- Lifecycle / Factory Methods ---

    Quaternion copy() const
    {
        return *this;
    }

    Quaternion createNormalized() const
    {
        const real_t len = length();
        if (len != 0)
        {
            return Quaternion(x / len, y / len, z / len, w / len);
        }
        return Quaternion(0,0,0,0);
    }

    Quaternion conjugate() const
    {
        return Quaternion(-x, -y, -z, w);
    }

    Quaternion inverse() const
    {
        real_t len_sq = lengthSquared();
        if (len_sq != 0)
        {
            return conjugate() / len_sq;
        }
        return Quaternion(0,0,0,0);
    }

    //

    Vector3D rotate(const Vector3D &vec) const
    {
        Quaternion q = *this * Quaternion(vec.x, vec.y, vec.z, 0) * conjugate();
        return Vector3D(q.x, q.y, q.z);
    }

    // --- Operators ---

    Quaternion operator*(const Quaternion &other) const
    {
        return Quaternion(
            this->w * other.x + this->x * other.w + this->y * other.z - this->z * other.y, // x
            this->w * other.y - this->x * other.z + this->y * other.w + this->z * other.x, // y
            this->w * other.z + this->x * other.y - this->y * other.x + this->z * other.w, // z
            this->w * other.w - this->x * other.x - this->y * other.y - this->z * other.z  // w
        );
    }

    Quaternion& operator*=(const Quaternion &other)
    {
        *this = *this * other;
        return *this;
    }

    Quaternion operator*(real_t scalar) const
    {
        return Quaternion(x * scalar, y * scalar, z * scalar, w * scalar);
    }

    Quaternion operator/(real_t scalar) const
    {
        return Quaternion(x / scalar, y / scalar, z / scalar, w / scalar);
    }

    Quaternion& operator*=(real_t scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    Quaternion& operator/=(real_t scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        w /= scalar;
        return *this;
    }

    // --- Comparison Operators ---

    bool operator==(const Quaternion &other) const
    {
        return approximatelyEqual(x, other.x) && approximatelyEqual(y, other.y) && approximatelyEqual(z, other.z) && approximatelyEqual(w, other.w);
    }

    bool operator!=(const Quaternion &other) const
    {
        return !approximatelyEqual(x, other.x) || !approximatelyEqual(y, other.y) || !approximatelyEqual(z, other.z) || !approximatelyEqual(w, other.w);
    }

    // --- Conversation between Quaternion and Matrix ---

    Matrix3x3 toMatrix3x3() const
    {
        real_t xx = x * x;
        real_t yy = y * y;
        real_t zz = z * z;
        real_t xy = x * y;
        real_t xz = x * z;
        real_t yz = y * z;
        real_t wx = w * x;
        real_t wy = w * y;
        real_t wz = w * z;

        return Matrix3x3(
            real_t(1) - real_t(2) * (yy + zz), real_t(2) * (xy - wz), real_t(2) * (xz + wy),
            real_t(2) * (xy + wz), real_t(1) - real_t(2) * (xx + zz), real_t(2) * (yz - wx),
            real_t(2) * (xz - wy), real_t(2) * (yz + wx), real_t(1) - real_t(2) * (xx + yy)
        );
    }

    Matrix4x4 toMatrix4x4() const
    {
        real_t xx = x * x;
        real_t yy = y * y;
        real_t zz = z * z;
        real_t xy = x * y;
        real_t xz = x * z;
        real_t yz = y * z;
        real_t wx = w * x;
        real_t wy = w * y;
        real_t wz = w * z;

        return Matrix4x4(
            real_t(1) - real_t(2) * (yy + zz), real_t(2) * (xy - wz), real_t(2) * (xz + wy), 0,
            real_t(2) * (xy + wz), real_t(1) - real_t(2) * (xx + zz), real_t(2) * (yz - wx), 0,
            real_t(2) * (xz - wy), real_t(2) * (yz + wx), real_t(1) - real_t(2) * (xx + yy), 0,
            0, 0, 0, 1
        );
    }
};

} // namespace Math

} // namespace Arns