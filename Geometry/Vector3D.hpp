#pragma once

#include <math.h>

namespace Utility
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
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator!=(const Vector3D& other) const
    {
        return x != other.x || y != other.y || z != other.z;
    }

    Vector3D operator-() const
    {
        return Vector3D(-x, -y, -z);
    }
};

} // namespace Math

} // namespace Utility