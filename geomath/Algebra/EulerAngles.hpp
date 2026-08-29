/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"

namespace Arns
{

namespace Math
{

class EulerAngles
{
public:
    real_t m_pitch;
    real_t m_yaw;
    real_t m_roll;

    EulerAngles() : m_pitch(0), m_yaw(0), m_roll(0) {}

    EulerAngles(real_t pitch, real_t yaw, real_t roll) : m_pitch(pitch), m_yaw(yaw), m_roll(roll) {}

    // --- Geometric Properties ---

    real_t length() const
    {
        return sqrt(m_pitch * m_pitch + m_yaw * m_yaw + m_roll * m_roll);
    }

    real_t lengthSquared() const
    {
        return m_pitch * m_pitch + m_yaw * m_yaw + m_roll * m_roll;
    }

    real_t distance(const EulerAngles &other) const
    {
        return (*this - other).length();
    }

    // --- State Queries ---

    bool isNormalized() const
    {
        return m_pitch >= -180 && m_pitch <= 180 && m_yaw >= -180 && m_yaw <= 180 && m_roll >= -180 && m_roll <= 180;
    }

    // --- Transform / modification ---

    EulerAngles& normalize()
    {
        m_pitch = wrapValue(m_pitch, -180.f, 180.f);
        m_yaw   = wrapValue(m_yaw,   -180.f, 180.f);
        m_roll  = wrapValue(m_roll,  -180.f, 180.f);
        return *this;
    }

    EulerAngles& clampPitch(real_t min, real_t max)
    {
        m_pitch = fmin(max, fmax(min, m_pitch));
        return *this;
    }

    EulerAngles& clampYaw(real_t min, real_t max)
    {
        m_yaw = fmin(max, fmax(min, m_yaw));
        return *this;
    }

    EulerAngles& clampRoll(real_t min, real_t max)
    {
        m_roll = fmin(max, fmax(min, m_roll));
        return *this;
    }

    EulerAngles lerp(const EulerAngles &other, real_t t) const
    {
        return EulerAngles(m_pitch + (other.m_pitch - m_pitch) * t, m_yaw + (other.m_yaw - m_yaw) * t, m_roll + (other.m_roll - m_roll) * t).normalize();
    }

    // --- Operators ---

    EulerAngles operator+(const EulerAngles &other) const
    {
        return EulerAngles(m_pitch + other.m_pitch, m_yaw + other.m_yaw, m_roll + other.m_roll).normalize();
    }

    EulerAngles operator-(const EulerAngles &other) const
    {
        return EulerAngles(m_pitch - other.m_pitch, m_yaw - other.m_yaw, m_roll - other.m_roll).normalize();
    }

    EulerAngles operator*(real_t scalar) const
    {
        return EulerAngles(m_pitch * scalar, m_yaw * scalar, m_roll * scalar).normalize();
    }

    EulerAngles operator/(real_t scalar) const
    {
        return EulerAngles(m_pitch / scalar, m_yaw / scalar, m_roll / scalar).normalize();
    }

    EulerAngles& operator+=(const EulerAngles &other)
    {
        m_pitch += other.m_pitch;
        m_yaw += other.m_yaw;
        m_roll += other.m_roll;
        return this->normalize();
    }

    EulerAngles& operator-=(const EulerAngles &other)
    {
        m_pitch -= other.m_pitch;
        m_yaw -= other.m_yaw;
        m_roll -= other.m_roll;
        return this->normalize();
    }

    EulerAngles& operator*=(real_t scalar)
    {
        m_pitch *= scalar;
        m_yaw *= scalar;
        m_roll *= scalar;
        return this->normalize();
    }

    EulerAngles& operator/=(real_t scalar)
    {
        m_pitch /= scalar;
        m_yaw /= scalar;
        m_roll /= scalar;
        return this->normalize();
    }

    // --- Comparison Operators ---

    bool operator==(const EulerAngles &other) const
    {
        return approximatelyEqual(m_pitch, other.m_pitch) && approximatelyEqual(m_yaw, other.m_yaw) && approximatelyEqual(m_roll, other.m_roll);
    }

    bool operator!=(const EulerAngles &other) const
    {
        return !approximatelyEqual(m_pitch, other.m_pitch) || !approximatelyEqual(m_yaw, other.m_yaw) || !approximatelyEqual(m_roll, other.m_roll);
    }

    // Conversion between EulerAngles and Matrix

    Matrix3x3 toMatrix3x3() const
    {
        return Matrix3x3::createRotationXRads(degToRad(m_pitch)) * 
            Matrix3x3::createRotationYRads(degToRad(m_yaw)) * 
            Matrix3x3::createRotationZRads(degToRad(m_roll));
    }

    Matrix4x4 toMatrix4x4() const
    {
        return Matrix4x4::createRotationXRads(degToRad(m_pitch)) * 
            Matrix4x4::createRotationYRads(degToRad(m_yaw)) * 
            Matrix4x4::createRotationZRads(degToRad(m_roll));
    }
};

} // namespace Math

} // namespace Arns