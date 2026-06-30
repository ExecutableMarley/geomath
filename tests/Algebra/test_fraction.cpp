#include "../third_party/doctest.h"
#include "Algebra/Fraction.hpp"

using namespace Arns::Math;

TEST_SUITE("Fraction")
{
    TEST_CASE("Default constructor")
    {
        Fraction f;

        CHECK(f.numerator == 0);
        CHECK(f.denominator == 1);
        CHECK(f.decimal() == doctest::Approx(0.0));
    }

    TEST_CASE("Construction and simplification")
    {
        CHECK(Fraction(2, 4) == Fraction(1, 2));
        CHECK(Fraction(10, 5) == Fraction(2));
        CHECK(Fraction(-6, 8) == Fraction(-3, 4));
        CHECK(Fraction(6, -8) == Fraction(-3, 4));
        CHECK(Fraction(-6, -8) == Fraction(3, 4));
        CHECK(Fraction(0, 5) == Fraction(0));
    }

    TEST_CASE("Construction with zero denominator throws")
    {
        CHECK_THROWS_AS(Fraction(1, 0), std::invalid_argument);
    }

    TEST_CASE("Decimal conversion")
    {
        CHECK(Fraction(1, 2).decimal() == doctest::Approx(0.5));
        CHECK(Fraction(-3, 4).decimal() == doctest::Approx(-0.75));
        CHECK(Fraction(5).decimal() == doctest::Approx(5.0));
    }

    TEST_CASE("Addition")
    {
        CHECK(Fraction(1, 2) + Fraction(1, 3) == Fraction(5, 6));
        CHECK(Fraction(3, 4) + Fraction(1, 4) == Fraction(1));
        CHECK(Fraction(-1, 2) + Fraction(1, 2) == Fraction(0));
    }

    TEST_CASE("Subtraction")
    {
        CHECK(Fraction(3, 4) - Fraction(1, 4) == Fraction(1, 2));
        CHECK(Fraction(1, 2) - Fraction(3, 2) == Fraction(-1));
        CHECK(Fraction(5, 7) - Fraction(5, 7) == Fraction(0));
    }

    TEST_CASE("Multiplication")
    {
        CHECK(Fraction(2, 3) * Fraction(3, 4) == Fraction(1, 2));
        CHECK(Fraction(-2, 3) * Fraction(3, 5) == Fraction(-2, 5));
        CHECK(Fraction(0) * Fraction(7, 8) == Fraction(0));
    }

    TEST_CASE("Division")
    {
        CHECK(Fraction(2, 3) / Fraction(4, 5) == Fraction(5, 6));
        CHECK(Fraction(-2, 3) / Fraction(4, 5) == Fraction(-5, 6));
        CHECK(Fraction(5, 7) / Fraction(5, 7) == Fraction(1));
    }

    TEST_CASE("Compound assignment operators")
    {
        Fraction f;

        f = Fraction(1, 2);
        f += Fraction(1, 3);
        CHECK(f == Fraction(5, 6));

        f = Fraction(3, 4);
        f -= Fraction(1, 4);
        CHECK(f == Fraction(1, 2));

        f = Fraction(2, 3);
        f *= Fraction(3, 4);
        CHECK(f == Fraction(1, 2));

        f = Fraction(2, 3);
        f /= Fraction(4, 5);
        CHECK(f == Fraction(5, 6));
    }

    TEST_CASE("createSimplified")
    {
        Fraction f(8, 12);
        Fraction s = f.createSimplified();

        CHECK(f.numerator == 2);
        CHECK(f.denominator == 3);

        CHECK(s == Fraction(2, 3));
    }

    TEST_CASE("fromDecimal with fixed precision")
    {
        CHECK(Fraction::fromFixedDecimal(0.5, 100) == Fraction(1, 2));
        CHECK(Fraction::fromFixedDecimal(0.25, 100) == Fraction(1, 4));
        CHECK(Fraction::fromFixedDecimal(1.75, 100) == Fraction(7, 4));
    }

    TEST_CASE("fromDecimal using continued fractions")
    {
        CHECK(Fraction::fromDecimal(0.5) == Fraction(1, 2));
        CHECK(Fraction::fromDecimal(0.25) == Fraction(1, 4));
        CHECK(Fraction::fromDecimal(1.75) == Fraction(7, 4));
        CHECK(Fraction::fromDecimal(-0.75) == Fraction(-3, 4));
        CHECK(Fraction::fromDecimal(2.0) == Fraction(2));
    }
}