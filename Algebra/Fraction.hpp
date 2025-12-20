#pragma once

#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "../CommonMath.hpp"

namespace Arns
{

namespace Math
{

class Fraction
{
public:
    int64_t numerator;
    int64_t denominator;

    Fraction() : numerator(0), denominator(1) {}

    Fraction(int64_t numerator, int64_t denominator = 1) : numerator(numerator), denominator(denominator)
    {
        if (denominator == 0)
        {
            throw std::invalid_argument("Denominator cannot be zero.");
        }
        simplify();
    }

    double decimal()
    {
        return (double)numerator / (double)denominator;
    }

    Fraction& simplify()
    {
        if (denominator < 0)
        {
            numerator   = -numerator;
            denominator = -denominator;
        }

        int g = gcd(numerator, denominator);
        if (g != 1)
        {
            this->numerator   /= g;
            this->denominator /= g;
        }
        return *this;
    }

    Fraction createSimplified() const
    {
        Fraction tmp = *this;
        tmp.simplify();
        return tmp;
    }

    int64_t gcd(int64_t a, int64_t b)
    {
        while (b != 0)
        {
            int64_t tmp = b;
            b = a % b;
            a = tmp;
        }
        return a;
    }

    //Todo: integer overflow

    Fraction operator+(const Fraction& other) const
    {
        const auto a = this->numerator * other.denominator;
        const auto b = other.numerator * this->denominator;
        return Fraction(a + b, this->denominator * other.denominator);
    }

    Fraction operator-(const Fraction& other) const
    {
        const auto a = this->numerator * other.denominator;
        const auto b = other.numerator * this->denominator;
        return Fraction(a - b, this->denominator * other.denominator);
    }

    Fraction operator*(const Fraction& other) const
    {
        return Fraction(this->numerator * other.numerator, this->denominator * other.denominator);
    }

    Fraction operator/(const Fraction& other) const
    {
        return Fraction(this->numerator * other.denominator, this->denominator * other.numerator);
    }

    Fraction& operator+=(const Fraction& other)
    {
        *this = *this + other;
        return *this;
    }

    Fraction& operator-=(const Fraction& other)
    {
        *this = *this - other;
        return *this;
    }

    Fraction& operator*=(const Fraction& other)
    {
        *this = *this * other;
        return *this;
    }

    Fraction& operator/=(const Fraction& other)
    {
        *this = *this / other;
        return *this;
    }

    bool operator==(const Fraction& other) const
    {
        auto a = this->createSimplified();
        auto b = other.createSimplified();
        return a.numerator == b.numerator && a.denominator == b.denominator;
    }

    Fraction operator!=(const Fraction& other) const
    {
        auto a = this->createSimplified();
        auto b = other.createSimplified();
        return a.numerator != b.numerator || a.denominator != b.denominator;
    }

    //Todo: Improve
    static Fraction fromDecimal(double decimal, int maxPrecision)
    {
        int64_t num = floor(decimal * maxPrecision);
        return Fraction(num, (int64_t)maxPrecision).simplify();
    }

    //Todo: Can this fail from normal input?
    static Fraction fromDecimal(double value, double tolerance = 1e-9)
    {
        bool negative = value < 0;
        value = std::abs(value);

        if (std::abs(value - std::round(value)) <= tolerance)
        return Fraction(negative ? -static_cast<int64_t>(std::llround(value))
                                 : static_cast<int64_t>(std::llround(value)),
                        (int64_t)1);

        double x = value;
        int64_t a = static_cast<int64_t>(std::floor(x));

        int64_t h0 = 0, h1 = 1;
        int64_t k0 = 1, k1 = 0;

        int64_t h2 = a, k2 = 1;

        constexpr size_t MAX_ITER = 100;
        for (size_t i = 0; i < MAX_ITER; ++i)
        {
            double approx = double(h2) / k2;
            /*
            if (std::abs(value - approx) <= tolerance * std::max(1.0, value))
                break;
                */
            if (approximatelyEqual(value, approx, tolerance, tolerance))
                break;

            double frac = x - a;
            if (approximatelyZero(frac))
                break;

            x = 1.0 / frac;
            a = static_cast<int64_t>(std::floor(x));

            h0 = h1; h1 = h2;
            k0 = k1; k1 = k2;

            //Overflow check
            if (std::abs(a) > (INT64_MAX - std::abs(h0)) / std::abs(h1))
                break;

            if (std::abs(a) > (INT64_MAX - std::abs(k0)) / std::abs(k1))
                break;

            h2 = a * h1 + h0;
            k2 = a * k1 + k0;
        }
        //Probably already simplified?
        return Fraction(negative ? -h2 : h2, k2).simplify();
    }

    //Todo: Force specific denominator
};

} // namespace Math

} // namespace Arns