#include "third_party/doctest.h"

#include "CommonMath.hpp"

using namespace Arns::Math;

// ============================================================================
// Approximations
// ============================================================================

TEST_CASE("approximatelyZero(float)")
{
    constexpr float eps = FloatAbsEpsilon;

    CHECK(approximatelyZero(0.0f));
    CHECK(approximatelyZero(+eps * 0.5f));
    CHECK(approximatelyZero(-eps * 0.5f));

    CHECK_FALSE(approximatelyZero(+eps * 1.1f));
    CHECK_FALSE(approximatelyZero(-eps * 1.1f));

    CHECK_FALSE(approximatelyZero(1.0f));
    CHECK_FALSE(approximatelyZero(-1.0f));
}

TEST_CASE("approximatelyZero(double)")
{
    constexpr double eps = DoubleAbsEpsilon;

    CHECK(approximatelyZero(0.0));
    CHECK(approximatelyZero(+eps * 0.5));
    CHECK(approximatelyZero(-eps * 0.5));

    CHECK_FALSE(approximatelyZero(+eps * 1.1));
    CHECK_FALSE(approximatelyZero(-eps * 1.1));

    CHECK_FALSE(approximatelyZero(1.0));
    CHECK_FALSE(approximatelyZero(-1.0));
}

TEST_CASE("approximatelyEqual - absolute epsilon near zero (float)")
{
    float a = 0.0f;
    float b = FloatAbsEpsilon * 0.9f;

    CHECK(approximatelyEqual(a, b));
    CHECK(approximatelyEqual(b, a));
}

TEST_CASE("approximatelyEqual - absolute epsilon near zero (double)")
{
    double a = 0.0;
    double b = DoubleAbsEpsilon * 0.9;

    CHECK(approximatelyEqual(a, b));
    CHECK(approximatelyEqual(b, a));
}

TEST_CASE("approximatelyEqual - relative epsilon large magnitude (float)")
{
    float a = 1e6f;
    float b = a + a * FloatRelEpsilon * 0.9f;

    CHECK(approximatelyEqual(a, b));
    CHECK(approximatelyEqual(b, a));
}

TEST_CASE("approximatelyEqual - relative epsilon large magnitude (double)")
{
    double a = 1e12;
    double b = a + a * DoubleRelEpsilon * 0.9;

    CHECK(approximatelyEqual(a, b));
    CHECK(approximatelyEqual(b, a));
}

TEST_CASE("approximatelyEqual - clearly unequal")
{
    CHECK_FALSE(approximatelyEqual(1.0f, 1.1f));
    CHECK_FALSE(approximatelyEqual(1.0, 1.1));
}

// ============================================================================
// 
// ============================================================================