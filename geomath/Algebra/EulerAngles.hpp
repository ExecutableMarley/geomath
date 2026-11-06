/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Formatting/Formatting.hpp"
#include "Algebra/Angle.hpp"
#include "Geometry/Vector3D.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"

namespace Arns
{

namespace geomath
{

class EulerAngles
{
public:
    Radians m_pitch;
    Radians m_yaw;
    Radians m_roll;

    EulerAngles() : m_pitch(Radians(0)), m_yaw(Radians(0)), m_roll(Radians(0)) {}

    EulerAngles(Radians pitch, Radians yaw, Radians roll) : m_pitch(pitch), m_yaw(yaw), m_roll(roll) {}

    // --- Geometric Properties ---

    real_t distance(const EulerAngles& other) const
    {
        const EulerAngles delta = shortestDeltaTo(other);
        return sqrt(delta.m_pitch.value * delta.m_pitch.value +
            delta.m_yaw.value * delta.m_yaw.value +
            delta.m_roll.value * delta.m_roll.value);
    }

    // --- State Queries ---

    bool isNormalized(Radians minVal = Radians(-PI), Radians maxVal = Radians(PI)) const
    {
        return m_pitch >= minVal && m_pitch <= maxVal &&
            m_yaw >= minVal && m_yaw <= maxVal &&
            m_roll >= minVal && m_roll <= maxVal;
    }

    // --- Transform / modification ---

    EulerAngles& normalize(Radians minVal = Radians(-PI), Radians maxVal = Radians(PI))
    {
        m_pitch.value = wrapValue(m_pitch.value, minVal.value, maxVal.value);
        m_yaw.value   = wrapValue(m_yaw.value,   minVal.value, maxVal.value);
        m_roll.value  = wrapValue(m_roll.value,  minVal.value, maxVal.value);
        return *this;
    }

    EulerAngles& clampPitch(Radians min, Radians max)
    {
        m_pitch.value = fmin(max.value, fmax(min.value, m_pitch.value));
        return *this;
    }

    EulerAngles& clampYaw(Radians min, Radians max)
    {
        m_yaw.value = fmin(max.value, fmax(min.value, m_yaw.value));
        return *this;
    }

    EulerAngles& clampRoll(Radians min, Radians max)
    {
        m_roll.value = fmin(max.value, fmax(min.value, m_roll.value));
        return *this;
    }

    EulerAngles lerp(const EulerAngles &other, real_t t) const
    {
        return EulerAngles(m_pitch + (other.m_pitch - m_pitch) * t, 
            m_yaw + (other.m_yaw - m_yaw) * t, 
            m_roll + (other.m_roll - m_roll) * t).normalize();
    }

    EulerAngles shortestDeltaTo(const EulerAngles& target) const
    {
        return EulerAngles(
            angleDifference(m_pitch, target.m_pitch),
            angleDifference(m_yaw, target.m_yaw),
            angleDifference(m_roll, target.m_roll));
    }

    // --- Operators ---

    EulerAngles operator+(const EulerAngles &other) const
    {
        return EulerAngles(Radians(m_pitch + other.m_pitch), Radians(m_yaw + other.m_yaw), Radians(m_roll + other.m_roll)).normalize();
    }

    EulerAngles operator-(const EulerAngles &other) const
    {
        return EulerAngles(Radians(m_pitch - other.m_pitch), Radians(m_yaw - other.m_yaw), Radians(m_roll - other.m_roll)).normalize();
    }

    EulerAngles operator*(real_t scalar) const
    {
        return EulerAngles(Radians(m_pitch * scalar), Radians(m_yaw * scalar), Radians(m_roll * scalar)).normalize();
    }

    EulerAngles operator/(real_t scalar) const
    {
        return EulerAngles(Radians(m_pitch / scalar), Radians(m_yaw / scalar), Radians(m_roll / scalar)).normalize();
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
        return approximatelyEqual(m_pitch.value, other.m_pitch.value) && approximatelyEqual(m_yaw.value, other.m_yaw.value) 
            && approximatelyEqual(m_roll.value, other.m_roll.value);
    }

    bool operator!=(const EulerAngles &other) const
    {
        return !approximatelyEqual(m_pitch.value, other.m_pitch.value) || !approximatelyEqual(m_yaw.value, other.m_yaw.value) 
            || !approximatelyEqual(m_roll.value, other.m_roll.value);
    }

    // --- Conversion between EulerAngles and Matrix ---

    Matrix3x3 toMatrix3x3() const
    {
        return Matrix3x3::createRotationX(m_pitch) * 
            Matrix3x3::createRotationY(m_yaw) * 
            Matrix3x3::createRotationZ(m_roll);
    }

    Matrix4x4 toMatrix4x4() const
    {
        return Matrix4x4::createRotationX(m_pitch) * 
            Matrix4x4::createRotationY(m_yaw) * 
            Matrix4x4::createRotationZ(m_roll);
    }

    // --- Stream Output ---

    friend std::ostream& operator<<(std::ostream& stream, const EulerAngles& angles)
    {
        return stream << "[" << angles.m_pitch.value << ", " << angles.m_yaw.value << ", " << angles.m_roll.value << "]";
    }
};

} // namespace geomath

} // namespace Arns

template <>
struct std::formatter<geomath::EulerAngles>
{
    int precision = 6;
    bool hasPrecision = false;

    constexpr auto parse(std::format_parse_context& context)
    {
        return geomath::parse_optional_float_format(context, precision, hasPrecision);
    }

    template <typename FormatContext>
    auto format(const geomath::EulerAngles& angles, FormatContext& context) const
    {
        if (hasPrecision)
            return std::format_to(context.out(), "[{:.{}f}, {:.{}f}, {:.{}f}]", angles.m_pitch.value, precision, angles.m_yaw.value, precision, angles.m_roll.value, precision);

        return std::format_to(context.out(), "[{}, {}, {}]", angles.m_pitch.value, angles.m_yaw.value, angles.m_roll.value);
    }
};