/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <cmath>
#include <algorithm>
#include <numbers>
#include <stdexcept>

namespace Arns
{

namespace Math
{


constexpr float FloatRelEpsilon = 1e-5f;
constexpr float FloatAbsEpsilon = 1e-6f;
constexpr double DoubleRelEpsilon = 1e-10;
constexpr double DoubleAbsEpsilon = 1e-12;

using real_t = float;

constexpr real_t PI = std::numbers::pi_v<real_t>;

constexpr real_t T_MAX = std::numeric_limits<real_t>::max();
constexpr real_t T_MIN = std::numeric_limits<real_t>::min();

inline bool approximatelyZero(float value, float absEpsilon = FloatAbsEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyZero(double value, double absEpsilon = DoubleAbsEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyEqual(float a, float b, 
                                float absEpsilon = FloatAbsEpsilon,
                                float relEpsilon = FloatRelEpsilon)
{
    float diff = std::fabs(a - b);
    if (diff <= absEpsilon)
        return true;

    return diff <= relEpsilon * std::max(std::fabs(a), std::fabs(b));
}

inline bool approximatelyEqual(double a, double b, 
                                double absEpsilon = DoubleAbsEpsilon,
                                double relEpsilon = DoubleRelEpsilon)
{
    double diff = std::fabs(a - b);
    if (diff <= absEpsilon)
        return true;

    return diff <= relEpsilon * std::max(std::fabs(a), std::fabs(b));
}

inline bool approximatelyGreater(float a, float b,
                                    float absEpsilon = FloatAbsEpsilon,
                                    float relEpsilon = FloatRelEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a > b;
}

inline bool approximatelyGreater(double a, double b,
                                    double absEpsilon = DoubleRelEpsilon,
                                    double relEpsilon = DoubleAbsEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a > b;
}

inline bool approximatelyLess(float a, float b,
                                float absEpsilon = FloatAbsEpsilon,
                                float relEpsilon = FloatRelEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a < b;
}

inline bool approximatelyLess(double a, double b,
                                double absEpsilon = DoubleRelEpsilon,
                                double relEpsilon = DoubleAbsEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a < b;
}

inline bool approximatelyZeroAbs(float value, float absEpsilon = FloatAbsEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyZeroAbs(double value, double absEpsilon = DoubleAbsEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyEqualAbs(float a, float b, float absEpsilon = FloatAbsEpsilon)
{
    return approximatelyZero(a - b, absEpsilon);
}

inline bool approximatelyEqualAbs(double a, double b, double absEpsilon = DoubleAbsEpsilon)
{
    return approximatelyZero(a - b, absEpsilon);
}

inline bool approximatelyGreaterAbs(float a, float b, float absEpsilon = FloatAbsEpsilon)
{
    return (a - b) > absEpsilon;
}

inline bool approximatelyGreaterAbs(double a, double b, double absEpsilon = DoubleAbsEpsilon)
{
    return (a - b) > absEpsilon;
}

inline bool approximatelyLessAbs(float a, float b, float epsilon = FloatAbsEpsilon)
{
    return (b - a) > epsilon;
}

inline bool approximatelyLessAbs(double a, double b, double epsilon = DoubleAbsEpsilon)
{
    return (b - a) > epsilon;
}

// Included in C++17
template <class T>
inline T clamp(T value, T minVal, T maxVal)
{
    return value < minVal ? minVal : (value > maxVal ? maxVal : value);
    return std::max(minVal, std::min(value, maxVal));
}

// Included in C++20
template <class T>
inline T lerp(T a, T b, T t)
{
    return a + (b - a) * t;
}

template <class T>
inline T inverseLerp(T a, T b, T value)
{
    return (value - a) / (b - a);
}

template <class T>
inline T remap(T value, T min1, T max1, T min2, T max2)
{
    return lerp(min2, max2, inverseLerp(min1, max1, value));
}

template <class T>
inline int sign(T value)
{
    return (value > T(0)) - (value < T(0));
}

template <class T>
inline T wrapValue(T value, T min, T max)
{
    const T range = max - min;
    if (range == T(0))
        return min;

    if constexpr (std::is_integral_v<T>)
    {
        value = (value - min) % range;
        if (value < 0)
            value += range;
    }
    else if constexpr (std::is_floating_point_v<T>)
    {
        value = std::fmod(value - min, range);
        if (value < T(0))
            value += range;
    }
    else
    {
        static_assert(std::is_arithmetic_v<T>, "wrapValue requires arithmetic types.");
    }
    return value + min;
}

template <typename T>
inline bool inInterval(const T& x, const T& minVal, const T& maxVal)
{
    return (minVal <= x && x <= maxVal);
}

template <typename T>
inline bool inIntervalExclusive(const T& x, const T& minVal, const T& maxVal)
{
    return (minVal < x && x < maxVal);
}

template <typename T>
inline bool intervalsOverlap(const T& minA, const T& maxA, const T& minB, const T& maxB)
{
    return !(maxA < minB || maxB < minA);
}

template <typename T>
inline T degToRad(T degrees)
{
    return degrees * (PI / real_t{180.0});
}

template <typename T>
inline T radToDeg(T radians)
{
    return radians * (real_t{180.0} / PI);
}

template <typename T>
inline T sinDeg(T degrees)
{
    return sin(degToRad(degrees));
}

template <typename T>
inline T cosDeg(T degrees)
{
    return cos(degToRad(degrees));
}

template <typename T>
inline T tanDeg(T degrees)
{
    return tan(degToRad(degrees));
}

template <typename T>
inline void sinCos(T radians, T &sine, T &cosine)
{
    sine = sin(radians);
    cosine = cos(radians);
}

template <typename T>
inline void sinCosDeg(T degrees, T &sine, T &cosine)
{
    sinCos(degToRad(degrees), sine, cosine);
}


// Check for GCC/Clang builtins
#if defined(__GNUC__) || defined(__clang__)
    #define HAS_BUILTIN_OVERFLOW 1
#else
    #define HAS_BUILTIN_OVERFLOW 0
#endif

// Overflow aware arithmetic

// Safe signed multiplication
template <std::integral T>
inline T safe_mul(T a, T b) 
{
#if HAS_BUILTIN_OVERFLOW
    T result;
    if (__builtin_mul_overflow(a, b, &result)) {
        throw std::overflow_error("Multiplication overflow");
    }
    return result;
#else
    // Software Fallback
    if (a == 0 || b == 0) return 0;
    if (a > 0)
    {
        if (b > 0) { if (a > std::numeric_limits<T>::max() / b) throw std::overflow_error("Multiplication Overflow"); }
        else { if (b < std::numeric_limits<T>::min() / a) throw std::overflow_error("Multiplication Overflow"); }
    } 
    else
    {
        if (b > 0) { if (a < std::numeric_limits<T>::min() / b) throw std::overflow_error("Multiplication Overflow"); }
        else {
            if (a == std::numeric_limits<T>::min() || b == std::numeric_limits<T>::min()) throw std::overflow_error("Multiplication Overflow");
            if (-a > std::numeric_limits<T>::max() / (-b)) throw std::overflow_error("Multiplication Overflow");
        }
    }
    return a * b;
#endif
}

// Safe signed addition
template <std::integral T>
inline T safe_add(T a, T b) 
{
#if HAS_BUILTIN_OVERFLOW
    T result;
    if (__builtin_add_overflow(a, b, &result)) {
        throw std::overflow_error("Addition overflow");
    }
    return result;
#else
    // Software Fallback
    if ((b > 0 && a > std::numeric_limits<T>::max() - b) ||
        (b < 0 && a < std::numeric_limits<T>::min() - b)) {
        throw std::overflow_error("Addition overflow");
    }
    return a + b;
#endif
}

// Safe signed subtraction
template <std::integral T>
inline T safe_sub(T a, T b) 
{
#if HAS_BUILTIN_OVERFLOW
    T result;
    if (__builtin_sub_overflow(a, b, &result)) {
        throw std::overflow_error("Subtraction overflow");
    }
    return result;
#else
    // Software Fallback
    if ((b > 0 && a < std::numeric_limits<T>::min() + b) ||
        (b < 0 && a > std::numeric_limits<T>::max() + b)) {
        throw std::overflow_error("Subtraction overflow");
    }
    return a - b;
#endif
}


} // namespace Math

} // namespace Arns