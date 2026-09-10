/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <cmath>
#include <algorithm>
#include <concepts>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>

namespace Arns
{

namespace geomath
{

/*
constexpr float FloatRelEpsilon = 1e-5f;
constexpr float FloatAbsEpsilon = 1e-6f;
constexpr double DoubleRelEpsilon = 1e-10;
constexpr double DoubleAbsEpsilon = 1e-12;
constexpr long double LongDoubleRelEpsilon = 1e-15L;
constexpr long double LongDoubleAbsEpsilon = 1e-18L;
*/

#ifndef GEOMATH_REAL_TYPE
#define GEOMATH_REAL_TYPE float
#endif

using real_t = GEOMATH_REAL_TYPE;

static_assert(std::is_same_v<real_t, float> ||
              std::is_same_v<real_t, double> ||
              std::is_same_v<real_t, long double>,
              "GEOMATH_REAL_TYPE must be float, double, or long double");

template <typename T>
struct RealTraits;

template <>
struct RealTraits<float>
{
    static constexpr float relEpsilon = 1e-5f;
    static constexpr float absEpsilon = 1e-6f;
};

template <>
struct RealTraits<double>
{
    static constexpr double relEpsilon = 1e-10;
    static constexpr double absEpsilon = 1e-12;
};

template <>
struct RealTraits<long double>
{
    static constexpr long double relEpsilon = 1e-15L;
    static constexpr long double absEpsilon = 1e-18L;
};

constexpr real_t RelEpsilon = RealTraits<real_t>::relEpsilon;
constexpr real_t AbsEpsilon = RealTraits<real_t>::absEpsilon;

constexpr real_t PI = std::numbers::pi_v<real_t>;

constexpr real_t T_MAX = std::numeric_limits<real_t>::max();
constexpr real_t T_MIN = std::numeric_limits<real_t>::min();

inline bool approximatelyZero(float value, float absEpsilon = RealTraits<float>::absEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyZero(double value, double absEpsilon = RealTraits<double>::absEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyZero(long double value, long double absEpsilon = RealTraits<long double>::absEpsilon)
{
    return std::fabs(value) < absEpsilon;
}

inline bool approximatelyEqual(float a, float b, 
                                float absEpsilon = RealTraits<float>::absEpsilon,
                                float relEpsilon = RealTraits<float>::relEpsilon)
{
    float diff = std::fabs(a - b);
    if (diff <= absEpsilon)
        return true;

    return diff <= relEpsilon * std::max(std::fabs(a), std::fabs(b));
}

inline bool approximatelyEqual(double a, double b, 
                                double absEpsilon = RealTraits<double>::absEpsilon,
                                double relEpsilon = RealTraits<double>::relEpsilon)
{
    double diff = std::fabs(a - b);
    if (diff <= absEpsilon)
        return true;

    return diff <= relEpsilon * std::max(std::fabs(a), std::fabs(b));
}

inline bool approximatelyEqual(long double a, long double b,
                                long double absEpsilon = RealTraits<long double>::absEpsilon,
                                long double relEpsilon = RealTraits<long double>::relEpsilon)
{
    long double diff = std::fabs(a - b);
    if (diff <= absEpsilon)
        return true;

    return diff <= relEpsilon * std::max(std::fabs(a), std::fabs(b));
}

inline bool approximatelyGreater(float a, float b,
                                    float absEpsilon = RealTraits<float>::absEpsilon,
                                    float relEpsilon = RealTraits<float>::relEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a > b;
}

inline bool approximatelyGreater(double a, double b,
                                    double absEpsilon = RealTraits<double>::absEpsilon,
                                    double relEpsilon = RealTraits<double>::relEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a > b;
}

inline bool approximatelyGreater(long double a, long double b,
                                  long double absEpsilon = RealTraits<long double>::absEpsilon,
                                  long double relEpsilon = RealTraits<long double>::relEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a > b;
}

inline bool approximatelyLess(float a, float b,
                                float absEpsilon = RealTraits<float>::absEpsilon,
                                float relEpsilon = RealTraits<float>::relEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a < b;
}

inline bool approximatelyLess(double a, double b,
                                double absEpsilon = RealTraits<double>::absEpsilon,
                                double relEpsilon = RealTraits<double>::relEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a < b;
}

inline bool approximatelyLess(long double a, long double b,
                               long double absEpsilon = RealTraits<long double>::absEpsilon,
                               long double relEpsilon = RealTraits<long double>::relEpsilon)
{
    if (approximatelyEqual(a, b, relEpsilon, absEpsilon))
    {
        return false;
    }
    return a < b;
}

inline bool approximatelyZeroAbs(float value, float absEpsilon = RealTraits<float>::absEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyZeroAbs(double value, double absEpsilon = RealTraits<double>::absEpsilon)
{
    return fabs(value) < absEpsilon;
}

inline bool approximatelyZeroAbs(long double value, long double absEpsilon = RealTraits<long double>::absEpsilon)
{
    return std::fabs(value) < absEpsilon;
}

inline bool approximatelyEqualAbs(float a, float b, float absEpsilon = RealTraits<float>::absEpsilon)
{
    return approximatelyZero(a - b, absEpsilon);
}

inline bool approximatelyEqualAbs(double a, double b, double absEpsilon = RealTraits<double>::absEpsilon)
{
    return approximatelyZero(a - b, absEpsilon);
}

inline bool approximatelyEqualAbs(long double a, long double b, long double absEpsilon = RealTraits<long double>::absEpsilon)
{
    return approximatelyZero(a - b, absEpsilon);
}

inline bool approximatelyGreaterAbs(float a, float b, float absEpsilon = RealTraits<float>::absEpsilon)
{
    return (a - b) > absEpsilon;
}

inline bool approximatelyGreaterAbs(double a, double b, double absEpsilon = RealTraits<double>::absEpsilon)
{
    return (a - b) > absEpsilon;
}

inline bool approximatelyGreaterAbs(long double a, long double b, long double absEpsilon = RealTraits<long double>::absEpsilon)
{
    return (a - b) > absEpsilon;
}

inline bool approximatelyLessAbs(float a, float b, float epsilon = RealTraits<float>::absEpsilon)
{
    return (b - a) > epsilon;
}

inline bool approximatelyLessAbs(double a, double b, double epsilon = RealTraits<double>::absEpsilon)
{
    return (b - a) > epsilon;
}

inline bool approximatelyLessAbs(long double a, long double b, long double epsilon = RealTraits<long double>::absEpsilon)
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

namespace geomath = Arns::geomath;