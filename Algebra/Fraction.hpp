#pragma once

#include <cmath>
#include <cstdint>
#include <stdexcept>

#define DefaultPrecision 5

namespace Utility
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

    Fraction(float decimal, int maxPrecision)
    {
        int64_t num = floor(decimal * maxPrecision);
        this->numerator = num;
        this->denominator = maxPrecision;
        this->simplify();
    }

    Fraction(double decimal, int maxPrecision)
    {
        int64_t num = floor(decimal * maxPrecision);
        this->numerator = num;
        this->denominator = maxPrecision;
        this->simplify();
    }

    float decimal()
    {
        return (float)numerator / (float)denominator;
    }

    void simplify()
    {
        int g = gcd(numerator, denominator);
        if (g != 1)
        {
            this->numerator   /= g;
            this->denominator /= g;
        }
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




    static Fraction fromDecimal(float decimal)
    {
        return Fraction(decimal, 10);
    }

    //Todo: Force specific denominator
};

} // namespace Math

} // namespace Utility