#pragma once

#include <math.h>
#include <algorithm>

namespace Utility
{

namespace Math
{


#define PI 3.14159265359f
#define EPSILON 0.00001f

bool approximatelyZero(float value, float epsilon = EPSILON)
{
    return fabs(value) < epsilon;
}

bool approximatelyEqual(float a, float b, float epsilon = EPSILON)
{
    return approximatelyZero(a - b, epsilon);
}

bool approximatelyGreater(float a, float b, float epsilon = EPSILON)
{
    return (a - b) > epsilon;
}

bool approximatelyLess(float a, float b, float epsilon = EPSILON)
{
    return (b - a) > epsilon;
}

// Included in C++17
float clamp(float value, float minVal, float maxVal)
{
    return std::max(minVal, std::min(value, maxVal));
}

// Included in C++20
float lerp(float a, float b, float t)
{
    return a + (b - a) * t;
}

float inverseLerp(float a, float b, float value)
{
    return (value - a) / (b - a);
}

float remap(float value, float min1, float max1, float min2, float max2)
{
    return lerp(min2, max2, inverseLerp(min1, max1, value));
}

int sign(float value)
{
    return (value > 0) - (value < 0);
}

float wrapValue(float value, float min, float max)
{
    const float range = max - min;
    while (value < min)
        value += range;
    while (value >= max)
        value -= range;
    return value;
}

float degToRad(float degrees)
{
    return degrees * (PI / 180.0f);
}

float radToDeg(float radians)
{
    return radians * (180.0f / PI);
}

float sinDeg(float degrees)
{
    return sin(degToRad(degrees));
}

float cosDeg(float degrees)
{
    return cos(degToRad(degrees));
}

float tanDeg(float degrees)
{
    return tan(degToRad(degrees));
}

void sinCos(float radians, float &sine, float &cosine)
{
    sine = sin(radians);
    cosine = cos(radians);
}

void sinCosDeg(float degrees, float &sine, float &cosine)
{
    sinCos(degToRad(degrees), sine, cosine);
}


} // namespace Math

} // namespace Utility