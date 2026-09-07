/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "IMatrix.hpp"
#include "Matrix.hpp"

#include <math.h>
#include <format>
#include <stdexcept>


namespace Arns
{

namespace geomath
{

// Unary operators

inline Matrix operator-(const IMatrix& matrix)
{
    Matrix result(matrix.rows(), matrix.columns());
    for (int i = 0; i < matrix.rows(); i++)
        for (int j = 0; j < matrix.columns(); j++)
            result(i, j) = -matrix(i, j);
    return result;
}

// Binary operators

// [*]

inline Matrix operator*(const IMatrix& lhs, const IMatrix& rhs)
{
    if (lhs.columns() != rhs.rows())
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");

    Matrix result(lhs.rows(), rhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < rhs.columns(); j++)
            for (int k = 0; k < lhs.columns(); k++)
                result(i, j) += lhs(i, k) * rhs(k, j);

    return result;
}

inline Matrix operator*(const IMatrix& lhs, real_t scalar)
{
    Matrix result(lhs.rows(), lhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            result(i, j) = lhs(i, j) * scalar;

    return result;
}

inline Matrix operator*(real_t scalar, const IMatrix& rhs)
{
    return rhs * scalar;
}

// [/]

inline Matrix operator/(const IMatrix& lhs, real_t scalar)
{
    Matrix result(lhs.rows(), lhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            result(i, j) = lhs(i, j) / scalar;

    return result;
}

// [+]

inline Matrix operator+(const IMatrix& lhs, const IMatrix& rhs)
{
    if (lhs.rows() != rhs.rows() || lhs.columns() != rhs.columns())
        throw std::invalid_argument("Matrix dimensions do not match for addition");

    Matrix result(lhs.rows(), lhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            result(i, j) = lhs(i, j) + rhs(i, j);
    return result;
}

inline Matrix operator+(const IMatrix& lhs, real_t scalar)
{
    Matrix result(lhs.rows(), lhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            result(i, j) = lhs(i, j) + scalar;

    return result;
}

inline Matrix operator+(real_t scalar, const IMatrix& rhs)
{
    return rhs + scalar;
    
}

// [-]

inline Matrix operator-(const IMatrix& lhs, const IMatrix& rhs)
{
    if (lhs.rows() != rhs.rows() || lhs.columns() != rhs.columns())
        throw std::invalid_argument("Matrix dimensions do not match for subtraction");

    Matrix result(lhs.rows(), lhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            result(i, j) = lhs(i, j) - rhs(i, j);
    return result;
}

inline Matrix operator-(const IMatrix& lhs, real_t scalar)
{
    Matrix result(lhs.rows(), lhs.columns());

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            result(i, j) = lhs(i, j) - scalar;

    return result;
}

// [==]

inline bool operator==(const IMatrix& lhs, const IMatrix& rhs)
{
    if (lhs.rows() != rhs.rows() || lhs.columns() != rhs.columns())
        return false;

    for (int i = 0; i < lhs.rows(); i++)
        for (int j = 0; j < lhs.columns(); j++)
            if (!approximatelyEqual(lhs(i, j), rhs(i, j)))
                return false;

    return true;
}

inline bool operator!=(const IMatrix& lhs, const IMatrix& rhs)
{
    return !(lhs == rhs);
}

inline std::ostream& operator<<(std::ostream &os, const geomath::IMatrix &matrix)
{
    if (matrix.rows() == 0 || matrix.columns() == 0) {
        return os << "[] (Empty Matrix)";
    }

    os << "[";

    for (size_t i = 0; i < matrix.rows(); ++i) {
        os << "[ ";
        for (size_t j = 0; j < matrix.columns(); ++j) {
            os << std::setw(8) << matrix(i, j);
            if (j < matrix.columns() - 1) {
                os << " ";
            }
        }
        os << " ]";

        if (i < matrix.rows() - 1) {
            os << ",";
        }
    }

    os << "]";

    return os;
}

} // namespace geomath

} // namespace Arns


template <>
struct std::formatter<geomath::IMatrix>
{
    constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const geomath::IMatrix& matrix, FormatContext& ctx) const
    {
        std::string result = "[";
        for (size_t i = 0; i < matrix.rows(); ++i)
        {
            result += "[";
            for (size_t j = 0; j < matrix.columns(); ++j)
            {
                result += std::to_string(matrix(i, j));
                if (j < matrix.columns() - 1)
                    result += ", ";
            }
            result += "]";
            if (i < matrix.rows() - 1)
                result += ", ";
        }
        result += "]";
        return std::format_to(ctx.out(), "{}", result);
    }
};