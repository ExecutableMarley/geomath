/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector3D.hpp"

namespace Utility
{

namespace Math
{

struct Quaternion
{
    float x;
    float y;
    float z;
    float w;

    Quaternion() : x(0), y(0), z(0), w(1) {}

    Quaternion(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

    float length() const
    {
        return sqrt(x * x + y * y + z * z + w * w);
    }

    float lengthSquared() const
    {
        return x * x + y * y + z * z + w * w;
    }

    Quaternion& normalize()
    {
        const float len = length();
        if (len != 0)
        {
            x /= len;
            y /= len;
            z /= len;
            w /= len;
        }
        return *this;
    }

    Quaternion createNormalized() const
    {
        const float len = length();
        if (len != 0)
        {
            return Quaternion(x / len, y / len, z / len, w / len);
        }
        return Quaternion(0,0,0,0);
    }

    Quaternion& resize(float newLength)
    {
        const float len = length();
        if (len != 0)
        {
            x *= newLength / len;
            y *= newLength / len;
            z *= newLength / len;
            w *= newLength / len;
        }
        return *this;
    }

    Quaternion conjugate() const
    {
        return Quaternion(-x, -y, -z, w);
    }

    Quaternion inverse() const
    {
        float len_sq = lengthSquared();
        if (len_sq != 0)
        {
            return conjugate() / len_sq;
        }
        return Quaternion(0,0,0,0);
    }

    Vector3D rotate(const Vector3D &vec) const
    {
        Quaternion q = *this * Quaternion(vec.x, vec.y, vec.z, 0) * conjugate();
        return Vector3D(q.x, q.y, q.z);
    }

    Quaternion operator*(const Quaternion &other) const
    {
        return Quaternion(
            this->w * other.w - this->x * other.x - this->y * other.y - this->z * other.z,
            this->w * other.x + this->x * other.w + this->y * other.z - this->z * other.y,
            this->w * other.y - this->x * other.z + this->y * other.w + this->z * other.x,
            this->w * other.z + this->x * other.y - this->y * other.x + this->z * other.w
        );
    }

    Quaternion& operator*=(const Quaternion &other)
    {
        *this = *this * other;
        return *this;
    }

    Quaternion operator*(float scalar) const
    {
        return Quaternion(x * scalar, y * scalar, z * scalar, w * scalar);
    }

    Quaternion operator/(float scalar) const
    {
        return Quaternion(x / scalar, y / scalar, z / scalar, w / scalar);
    }

    Quaternion& operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        w *= scalar;
        return *this;
    }

    Quaternion& operator/=(float scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        w /= scalar;
        return *this;
    }
};

} // namespace Math

} // namespace Utility