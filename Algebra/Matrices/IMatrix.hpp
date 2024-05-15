/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <stdexcept>

namespace Utility
{

namespace Math
{

class IMatrix
{
public:
    virtual ~IMatrix() = default;

    virtual size_t rows() const = 0;
    virtual size_t columns() const = 0;

    virtual float& operator()(size_t row, size_t column) = 0;
    virtual const float& operator()(size_t row, size_t column) const = 0;


    IMatrix& operator+=(const IMatrix& other)
    {
        if (rows() != other.rows() || columns() != other.columns())
            throw std::invalid_argument("Matrix dimensions do not match for addition");

        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < columns(); j++)
                (*this)(i, j) += other(i, j);
        return *this;
    }

    IMatrix& operator+=(float scalar)
    {
        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < columns(); j++)
                (*this)(i, j) += scalar;
        return *this;
    }

    IMatrix& operator-=(const IMatrix& other)
    {
        if (rows() != other.rows() || columns() != other.columns())
            throw std::invalid_argument("Matrix dimensions do not match for subtraction");

        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < columns(); j++)
                (*this)(i, j) -= other(i, j);
        return *this;
    }

    IMatrix& operator-=(float scalar)
    {
        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < columns(); j++)
                (*this)(i, j) -= scalar;
        return *this;
    }

    IMatrix& operator*=(float scalar)
    {
        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < columns(); j++)
                (*this)(i, j) *= scalar;
        return *this;
    }

    IMatrix& operator/=(float scalar)
    {
        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < columns(); j++)
                (*this)(i, j) /= scalar;
        return *this;
    }
};


} // namespace Math

} // namespace Utility