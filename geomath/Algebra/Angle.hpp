/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"

namespace Arns::geomath
{



struct Radians;

struct Degrees {
    real_t value;
    explicit constexpr Degrees(real_t v) : value(v) {}
    
    //Implicit Degree -> Radian
    operator Radians() const; 
};

struct Radians
{
    real_t value;
    explicit constexpr Radians(real_t v) : value(v) {}

    constexpr Radians operator+(Radians other) const { return Radians(value + other.value); }
    constexpr Radians operator-(Radians other) const { return Radians(value - other.value); }
    constexpr Radians operator*(real_t scalar) const { return Radians(value * scalar); }
    constexpr Radians operator/(real_t scalar) const { return Radians(value / scalar); }

    constexpr Radians& operator+=(Radians other) { value += other.value; return *this; }
    constexpr Radians& operator-=(Radians other) { value -= other.value; return *this; }
    constexpr Radians& operator*=(real_t scalar) { value *= scalar; return *this; }
    constexpr Radians& operator/=(real_t scalar) { value /= scalar; return *this; }

    constexpr bool operator==(Radians other) const { return value == other.value; }
    constexpr bool operator!=(Radians other) const { return value != other.value; }
    constexpr bool operator<(Radians other) const { return value < other.value; }
    constexpr bool operator<=(Radians other) const { return value <= other.value; }
    constexpr bool operator>(Radians other) const { return value > other.value; }
    constexpr bool operator>=(Radians other) const { return value >= other.value; }

    constexpr Radians operator-() const { return Radians(-value); }
    constexpr Radians abs() const { return Radians(std::abs(value)); }

    // Explicit! Degree -> Radian
    explicit operator Degrees() const;

    static constexpr Radians zero()   { return Radians(real_t(0)); }
    static constexpr Radians halfPi() { return Radians(PI * real_t(0.5)); }
    static constexpr Radians pi()     { return Radians(PI); }
    static constexpr Radians twoPi()  { return Radians(PI * real_t(2)); }
};

constexpr Radians operator*(real_t scalar, Radians angle) { return angle * scalar; }

inline Radians angleDifference(Radians from, Radians to)
{
    const real_t fullTurn = real_t(2) * PI;
    return Radians(std::remainder(to.value - from.value, fullTurn));
}

/*
inline real_t sin(Radians rad) { return std::sin(rad.value); }
inline real_t cos(Radians rad) { return std::cos(rad.value); }
inline real_t tan(Radians rad) { return std::tan(rad.value); }
*/

// lerp, inverseLerp, remap defined in CommonMath.hpp

namespace literals
{
    // --- QOL Literals for convenience ---
    constexpr Degrees operator""_deg(long double v)        { return Degrees(static_cast<real_t>(v)); }
    constexpr Degrees operator""_deg(unsigned long long v) { return Degrees(static_cast<real_t>(v)); }

    constexpr Radians operator""_rad(long double v)        { return Radians(static_cast<real_t>(v)); }
    constexpr Radians operator""_rad(unsigned long long v) { return Radians(static_cast<real_t>(v)); }
}

inline Degrees::operator Radians() const { return Radians(degToRad(value)); }
inline Radians::operator Degrees() const { return Degrees(radToDeg(value)); }



}