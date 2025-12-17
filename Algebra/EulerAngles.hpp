/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"

namespace Arns
{

namespace Math
{

class EulerAngles
{
public:
    float m_pitch;
    float m_yaw;
    float m_roll;

    EulerAngles() : m_pitch(0), m_yaw(0), m_roll(0) {}

    EulerAngles(float pitch, float yaw, float roll) : m_pitch(pitch), m_yaw(yaw), m_roll(roll) {}

    float length() const
    {
        return sqrt(m_pitch * m_pitch + m_yaw * m_yaw + m_roll * m_roll);
    }

    float lengthSquared() const
    {
        return m_pitch * m_pitch + m_yaw * m_yaw + m_roll * m_roll;
    }

    EulerAngles& normalize()
    {
        m_pitch = wrapValue(m_pitch, -180.f, 180.f);
        m_yaw   = wrapValue(m_yaw,   -180.f, 180.f);
        m_roll  = wrapValue(m_roll,  -180.f, 180.f);
        return *this;
    }
    
    bool isNormalized() const
    {
        return m_pitch >= -180 && m_pitch <= 180 && m_yaw >= -180 && m_yaw <= 180 && m_roll >= -180 && m_roll <= 180;
    }

    EulerAngles& clampPitch(float min, float max)
    {
        m_pitch = fmin(max, fmax(min, m_pitch));
        return *this;
    }

    EulerAngles& clampYaw(float min, float max)
    {
        m_yaw = fmin(max, fmax(min, m_yaw));
        return *this;
    }

    EulerAngles& clampRoll(float min, float max)
    {
        m_roll = fmin(max, fmax(min, m_roll));
        return *this;
    }

    EulerAngles lerp(const EulerAngles &other, float t) const
    {
        return EulerAngles(m_pitch + (other.m_pitch - m_pitch) * t, m_yaw + (other.m_yaw - m_yaw) * t, m_roll + (other.m_roll - m_roll) * t).normalize();
    }

    float distance(const EulerAngles &other) const
    {
        return (*this - other).length();
    }

    EulerAngles operator+(const EulerAngles &other) const
    {
        return EulerAngles(m_pitch + other.m_pitch, m_yaw + other.m_yaw, m_roll + other.m_roll).normalize();
    }

    EulerAngles operator-(const EulerAngles &other) const
    {
        return EulerAngles(m_pitch - other.m_pitch, m_yaw - other.m_yaw, m_roll - other.m_roll).normalize();
    }

    EulerAngles operator*(float scalar) const
    {
        return EulerAngles(m_pitch * scalar, m_yaw * scalar, m_roll * scalar).normalize();
    }

    EulerAngles operator/(float scalar) const
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

    EulerAngles& operator*=(float scalar)
    {
        m_pitch *= scalar;
        m_yaw *= scalar;
        m_roll *= scalar;
        return this->normalize();
    }

    EulerAngles& operator/=(float scalar)
    {
        m_pitch /= scalar;
        m_yaw /= scalar;
        m_roll /= scalar;
        return this->normalize();
    }

    bool operator==(const EulerAngles &other) const
    {
        return approximatelyEqual(m_pitch, other.m_pitch) && approximatelyEqual(m_yaw, other.m_yaw) && approximatelyEqual(m_roll, other.m_roll);
    }

    bool operator!=(const EulerAngles &other) const
    {
        return !approximatelyEqual(m_pitch, other.m_pitch) || !approximatelyEqual(m_yaw, other.m_yaw) || !approximatelyEqual(m_roll, other.m_roll);
    }
};

} // namespace Math

} // namespace Arns