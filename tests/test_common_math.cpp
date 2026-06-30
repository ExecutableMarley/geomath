#include "third_party/doctest.h"

#include "CommonMath.hpp"
#include <stdexcept>
#include <cstdint>

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

TEST_SUITE("Safe Integer Arithmetic")
{
    TEST_CASE("safe_add")
    {
        CHECK(safe_add<int64_t>(0, 0) == 0);
        CHECK(safe_add<int64_t>(1, 2) == 3);
        CHECK(safe_add<int64_t>(-1, -2) == -3);
        CHECK(safe_add<int64_t>(5, -3) == 2);
        CHECK(safe_add<int64_t>(-5, 3) == -2);

        CHECK(safe_add<int64_t>(INT64_MAX, 0) == INT64_MAX);
        CHECK(safe_add<int64_t>(INT64_MIN, 0) == INT64_MIN);

        CHECK(safe_add<int64_t>(INT64_MAX - 1, 1) == INT64_MAX);
        CHECK(safe_add<int64_t>(INT64_MIN + 1, -1) == INT64_MIN);

        CHECK_THROWS_AS(safe_add<int64_t>(INT64_MAX, 1), std::overflow_error);
        CHECK_THROWS_AS(safe_add<int64_t>(INT64_MAX - 10, 11), std::overflow_error);

        CHECK_THROWS_AS(safe_add<int64_t>(INT64_MIN, -1), std::overflow_error);
        CHECK_THROWS_AS(safe_add<int64_t>(INT64_MIN + 10, -11), std::overflow_error);
    }

    TEST_CASE("safe_sub")
    {
        CHECK(safe_sub<int64_t>(0, 0) == 0);
        CHECK(safe_sub<int64_t>(5, 3) == 2);
        CHECK(safe_sub<int64_t>(-5, -3) == -2);
        CHECK(safe_sub<int64_t>(5, -3) == 8);
        CHECK(safe_sub<int64_t>(-5, 3) == -8);

        CHECK(safe_sub<int64_t>(INT64_MAX, 0) == INT64_MAX);
        CHECK(safe_sub<int64_t>(INT64_MIN, 0) == INT64_MIN);

        CHECK(safe_sub<int64_t>(INT64_MAX, INT64_MAX) == 0);
        CHECK(safe_sub<int64_t>(INT64_MIN, INT64_MIN) == 0);

        CHECK_THROWS_AS(safe_sub<int64_t>(INT64_MAX, -1), std::overflow_error);

        CHECK_THROWS_AS(safe_sub<int64_t>(INT64_MIN, 1), std::overflow_error);
        CHECK_THROWS_AS(safe_sub<int64_t>(INT64_MAX, INT64_MIN), std::overflow_error);
        CHECK_THROWS_AS(safe_sub<int64_t>(INT64_MIN, INT64_MAX), std::overflow_error);
    }

    TEST_CASE("safe_mul")
    {
        CHECK(safe_mul<int64_t>(0, 0) == 0);
        CHECK(safe_mul<int64_t>(0, 123456) == 0);
        CHECK(safe_mul<int64_t>(123456, 0) == 0);

        CHECK(safe_mul<int64_t>(2, 3) == 6);
        CHECK(safe_mul<int64_t>(-2, 3) == -6);
        CHECK(safe_mul<int64_t>(2, -3) == -6);
        CHECK(safe_mul<int64_t>(-2, -3) == 6);

        CHECK(safe_mul<int64_t>(1, INT64_MAX) == INT64_MAX);
        CHECK(safe_mul<int64_t>(INT64_MAX, 1) == INT64_MAX);

        CHECK(safe_mul<int64_t>(-1, INT64_MAX) == -INT64_MAX);
        CHECK(safe_mul<int64_t>(INT64_MAX, -1) == -INT64_MAX);

        CHECK(safe_mul<int64_t>(INT64_MIN, 1) == INT64_MIN);
        CHECK(safe_mul<int64_t>(1, INT64_MIN) == INT64_MIN);

        CHECK(safe_mul<int64_t>(INT64_MAX / 2, 2) == (INT64_MAX / 2) * 2);
        CHECK(safe_mul<int64_t>(INT64_MIN / 2, 2) == (INT64_MIN / 2) * 2);

        CHECK_THROWS_AS(safe_mul<int64_t>(INT64_MAX, 2), std::overflow_error);
        CHECK_THROWS_AS(safe_mul<int64_t>(2, INT64_MAX), std::overflow_error);

        CHECK_THROWS_AS(safe_mul<int64_t>(INT64_MIN, -1), std::overflow_error);
        CHECK_THROWS_AS(safe_mul<int64_t>(-1, INT64_MIN), std::overflow_error);

        CHECK_THROWS_AS(safe_mul<int64_t>(INT64_MIN, 2), std::overflow_error);
        CHECK_THROWS_AS(safe_mul<int64_t>(2, INT64_MIN), std::overflow_error);

        CHECK_THROWS_AS(safe_mul<int64_t>(INT64_MAX, INT64_MAX), std::overflow_error);
        CHECK_THROWS_AS(safe_mul<int64_t>(INT64_MIN, INT64_MIN), std::overflow_error);
    }
}