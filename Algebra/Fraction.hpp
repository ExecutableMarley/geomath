#pragma once

#include <cmath>

#define DefaultPrecision 5

namespace Utility
{

namespace Math
{

//Todo: Currently not able to represent double precision
//Consider using long instead of int

class Fraction
{
public:
    int numerator;
    int denominator;

    Fraction() : numerator(0), denominator(1) {}

    //Prevent denominator == 0
    Fraction(int numerator, int denominator = 1) : numerator(numerator), denominator(denominator) {}

    Fraction(float decimal, int maxPrecision)
    {
        int num = floor(decimal * maxPrecision);
        this->numerator = num;
        this->denominator = maxPrecision;
        this->simplify();
    }

    Fraction(double decimal, int maxPrecision)
    {
        int num = floor(decimal * maxPrecision);
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

    int gcd(int a, int b)
    {
        while (b != 0)
        {
            int tmp = b;
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