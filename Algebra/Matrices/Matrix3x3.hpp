/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <stdexcept>

#include "IMatrix.hpp"

namespace Utility
{

namespace Math
{

class Matrix3x3 : public IMatrix
{
protected:
    float m_data[3][3];
public:

    Matrix3x3(float initValue = 0.0f)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m_data[i][j] = initValue;
    }

    Matrix3x3(float m00, float m01, float m02,
              float m10, float m11, float m12,
              float m20, float m21, float m22)
    {
        m_data[0][0] = m00;
        m_data[0][1] = m01;
        m_data[0][2] = m02;

        m_data[1][0] = m10;
        m_data[1][1] = m11;
        m_data[1][2] = m12;

        m_data[2][0] = m20;
        m_data[2][1] = m21;
        m_data[2][2] = m22;
    }

    Matrix3x3(const Matrix3x3 &other)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m_data[i][j] = other.m_data[i][j];
    }

    Matrix3x3(const IMatrix &other)
    {
        if (other.rows() != 3 || other.columns() != 3)
            throw std::invalid_argument("Matrix must be 3x3.");
            
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m_data[i][j] = other(i, j);
    }

    // IMatrix Interface

    size_t rows()    const override { return 3; }
    size_t columns() const override { return 3; }

    float& operator()(size_t row, size_t column) override
    {
        return m_data[row][column];
    }

    const float& operator()(size_t row, size_t column) const override
    {
        return m_data[row][column];
    }

    //

    Matrix3x3 transpose() const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[j][i];
        return result;
    }

    // Operators

    Matrix3x3 operator-() const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = -m_data[i][j];
        return result;
    } 

    //[+]

    Matrix3x3 operator+(const IMatrix &other) const
    {
        if (other.rows() != 3 || other.columns() != 3)
            throw std::invalid_argument("Matrix dimensions do not match for addition");

        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[i][j] + other(i, j);
        return result;
    }

    Matrix3x3 operator+(float scalar) const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[i][j] + scalar;
        return result;
    }

    //[-]

    Matrix3x3 operator-(const IMatrix &other) const
    {
        if (other.rows() != 3 || other.columns() != 3)
            throw std::invalid_argument("Matrix dimensions do not match for subtraction");

        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[i][j] - other(i, j);
        return result;
    }

    Matrix3x3 operator-(float scalar) const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[i][j] - scalar;
        return result;
    }   

    //[*]

    Matrix3x3 operator*(const Matrix3x3 &other) const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                result(i, j) = 0;
                for (int k = 0; k < 3; k++)
                {
                    result(i, j) += m_data[i][k] * other.m_data[k][j];
                }
            }
        }
        return result;
    }

    Matrix3x3 operator*(float scalar) const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[i][j] * scalar;
        return result;
    }

    //[/]

    Matrix3x3 operator/(float scalar) const
    {
        Matrix3x3 result;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                result(i, j) = m_data[i][j] / scalar;
        return result;
    }
};

} // namespace Math

} // namespace Utility