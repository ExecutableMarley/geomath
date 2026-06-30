#pragma once

#include "../CommonMath.hpp"

#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>

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

    double decimal() const
    {
        return (double)numerator / (double)denominator;
    }

    Fraction& simplify()
    {
        if (numerator == 0)
        {
            denominator = 1;
            return *this;
        }

        if (denominator < 0)
        {
            numerator   = -numerator;
            denominator = -denominator;
        }

        int g = std::gcd(numerator, denominator);
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

    Fraction operator+(const Fraction& other) const
    {
        int64_t g = std::gcd(this->denominator, other.denominator);

        int64_t d1 = this->denominator / g;
        int64_t d2 = other.denominator / g;

        int64_t a = safe_mul<int64_t>(this->numerator, d2);
        int64_t b = safe_mul<int64_t>(other.numerator, d1);

        int64_t num = safe_add(a, b);
        int64_t den = safe_mul(this->denominator, d2);

        return Fraction(num, den);
    }

    Fraction operator-(const Fraction& other) const
    {
        int64_t g = std::gcd(this->denominator, other.denominator);

        int64_t d1 = this->denominator / g;
        int64_t d2 = other.denominator / g;

        int64_t a = safe_mul(this->numerator, d2);
        int64_t b = safe_mul(other.numerator, d1);

        int64_t num = safe_sub(a, b);
        int64_t den = safe_mul(this->denominator, d2);

        return Fraction(num, den);
    }

    Fraction operator*(const Fraction &other) const
    {
        int64_t g1 = std::gcd(this->numerator, other.denominator);
        int64_t g2 = std::gcd(other.numerator, this->denominator);

        int64_t num = safe_mul(this->numerator / g1, other.numerator / g2);
        int64_t den = safe_mul(this->denominator / g2, other.denominator / g1);

        return Fraction(num, den);
    }

    Fraction operator/(const Fraction &other) const
    {
        if (other.numerator == 0)
            throw std::domain_error("Division by zero fraction");

        // (num1 / den1) * (den2 / num2)
        int64_t g1 = std::gcd(this->numerator, other.numerator);
        int64_t g2 = std::gcd(other.denominator, this->denominator);

        int64_t num = safe_mul(this->numerator / g1, other.denominator / g2);
        int64_t den = safe_mul(this->denominator / g2, other.numerator / g1);

        if (den < 0)
        {
            // Edge case, will overflow
            if (num == std::numeric_limits<int64_t>::min())
                throw std::overflow_error("Overflow");
            num = -num;
            den = -den;
        }

        return Fraction(num, den);
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

    static Fraction fromFixedDecimal(double decimal, int maxPrecision)
    {
        int64_t num = floor(decimal * maxPrecision);
        return Fraction(num, (int64_t)maxPrecision).simplify();
    }

    static Fraction fromDecimal(double value, double tolerance = 1e-9)
    {
        // NaN Inf protection
        if (std::isnan(value) || std::isinf(value))
        {
            return Fraction(0, 1); 
        }

        bool negative = value < 0;
        value = std::abs(value);

        // Fast exit for close integers
        if (std::abs(value - std::round(value)) <= tolerance)
        {
            int64_t intPart = static_cast<int64_t>(std::llround(value));
            return Fraction(negative ? -intPart : intPart, 1);
        }

        // Continued Fraction algorithm
        double x = value;
        int64_t h0 = 0, h1 = 1; // Numerators
        int64_t k0 = 1, k1 = 0; // Denominators

        constexpr size_t MAX_ITER = 100;
        for (size_t i = 0; i < MAX_ITER; ++i)
        {
            double floor_x = std::floor(x);
            
            // Overflow Safety check
            if (floor_x > static_cast<double>(std::numeric_limits<int64_t>::max())) {
                break;
            }
            int64_t a = static_cast<int64_t>(floor_x);

            // Double-check overflow
            if (h1 > 0 && a > (std::numeric_limits<int64_t>::max() - h0) / h1) break;
            if (k1 > 0 && a > (std::numeric_limits<int64_t>::max() - k0) / k1) break;

            int64_t h2 = a * h1 + h0;
            int64_t k2 = a * k1 + k0;

            // Update history
            h0 = h1; h1 = h2;
            k0 = k1; k1 = k2;

            // Check if current approximation is accurate enough
            double approx = static_cast<double>(h1) / k1;
            if (approximatelyEqualAbs(value, approx, tolerance)) {
                break;
            }

            // Prepare for next iteration
            double frac = x - floor_x;
            if (approximatelyZero(frac))
            {
                break;
            }
            x = 1.0 / frac;
        }
        // h1 and k1 hold the most recent safe convergent
        return Fraction(negative ? -h1 : h1, k1).simplify();
    }
};

} // namespace Math

} // namespace Arns