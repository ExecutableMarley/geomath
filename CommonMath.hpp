/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <cmath>
#include <algorithm>

namespace Arns
{

namespace Math
{


#define PI 3.14159265359f

constexpr float FloatRelEpsilon = 1e-5f;
constexpr float FloatAbsEpsilon = 1e-6f;
constexpr double DoubleRelEpsilon = 1e-10;
constexpr double DoubleAbsEpsilon = 1e-12;

using real_t = float;

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
inline float clamp(float value, float minVal, float maxVal)
{
    return value < minVal ? minVal : (value > maxVal ? maxVal : value);
    return std::max(minVal, std::min(value, maxVal));
}

// Included in C++20
inline float lerp(float a, float b, float t)
{
    return a + (b - a) * t;
}

inline float inverseLerp(float a, float b, float value)
{
    return (value - a) / (b - a);
}

inline float remap(float value, float min1, float max1, float min2, float max2)
{
    return lerp(min2, max2, inverseLerp(min1, max1, value));
}

inline int sign(float value)
{
    return (value > 0) - (value < 0);
}

inline int wrapValue(int value, int min, int max)
{
    const int range = max - min;
    if (range == 0)
        return min;
    value = (value - min) % range;
    if (value < 0) 
        value += range;
    
    return value + min;
}

inline float wrapValue(float value, float min, float max)
{
    const float range = max - min;
    if (range == 0) 
        return min;
    value = fmod(value - min, range);
    if (value < 0) 
        value += range;
    
    return value + min;
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

inline float degToRad(float degrees)
{
    return degrees * (PI / 180.0f);
}

inline float radToDeg(float radians)
{
    return radians * (180.0f / PI);
}

inline float sinDeg(float degrees)
{
    return sin(degToRad(degrees));
}

inline float cosDeg(float degrees)
{
    return cos(degToRad(degrees));
}

inline float tanDeg(float degrees)
{
    return tan(degToRad(degrees));
}

inline void sinCos(float radians, float &sine, float &cosine)
{
    sine = sin(radians);
    cosine = cos(radians);
}

inline void sinCosDeg(float degrees, float &sine, float &cosine)
{
    sinCos(degToRad(degrees), sine, cosine);
}


} // namespace Math

} // namespace Arns