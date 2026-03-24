/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <stdexcept>

#include "CommonMath.hpp"
#include "IMatrix.hpp"
#include "Matrix.hpp"

namespace Arns
{

namespace Math
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


} // namespace Math

} // namespace Arns