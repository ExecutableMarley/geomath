/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <stdexcept>

#include "CommonMath.hpp"
#include "IMatrix.hpp"
#include "Geometry/Vector3D.hpp"

namespace Arns
{

namespace Math
{

class Matrix4x4 : public IMatrix
{
protected:
    real_t m_data[4][4];
public:

    Matrix4x4(real_t initValue = 0.0f)
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                m_data[i][j] = initValue;
    }

    Matrix4x4(real_t m00, real_t m01, real_t m02, real_t m03,
              real_t m10, real_t m11, real_t m12, real_t m13,
              real_t m20, real_t m21, real_t m22, real_t m23,
              real_t m30, real_t m31, real_t m32, real_t m33)
    {
        m_data[0][0] = m00;
        m_data[0][1] = m01;
        m_data[0][2] = m02;
        m_data[0][3] = m03;

        m_data[1][0] = m10;
        m_data[1][1] = m11;
        m_data[1][2] = m12;
        m_data[1][3] = m13;

        m_data[2][0] = m20;
        m_data[2][1] = m21;
        m_data[2][2] = m22;
        m_data[2][3] = m23;

        m_data[3][0] = m30;
        m_data[3][1] = m31;
        m_data[3][2] = m32;
        m_data[3][3] = m33;
    }

    Matrix4x4(const IMatrix& other)
    {
        if (other.rows() != 4 || other.columns() != 4)
            throw std::invalid_argument("Matrix must be 4x4.");

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                m_data[i][j] = other(i, j);
    }

    // IMatrix Interface

    size_t rows()    const override { return 4; }
    size_t columns() const override { return 4; }

    real_t &operator()(size_t row, size_t column) override
    {
        return m_data[row][column];
    }

    const real_t &operator()(size_t row, size_t column) const override
    {
        return m_data[row][column];
    }

    //

    Matrix4x4 transpose() const
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[j][i];
        return result;
    }

    //

    Vector3D transform(const Vector3D& vec) const
    {
        return Vector3D(
            m_data[0][0] * vec.x + m_data[0][1] * vec.y + m_data[0][2] * 1,
            m_data[1][0] * vec.x + m_data[1][1] * vec.y + m_data[1][2] * 1,
            m_data[2][0] * vec.x + m_data[2][1] * vec.y + m_data[2][2] * 1);
    }

    real_t* operator[](size_t row)
    {
        return m_data[row];
    }

    // Operators

    Matrix4x4 operator-()
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = -m_data[i][j];
        return result;
    }

    //[+]

    Matrix4x4 operator+(const IMatrix& other) const
    {
        if (other.rows() != 4 || other.columns() != 4)
            throw std::invalid_argument("Matrix dimensions do not match for addition");

        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][j] + other(i, j);
        return result;
    }

    Matrix4x4 operator+(real_t scalar) const
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][j] + scalar;
        return result;
    }

    //[-]

    Matrix4x4 operator-(const IMatrix& other) const
    {
        if (other.rows() != 4 || other.columns() != 4)
            throw std::invalid_argument("Matrix dimensions do not match for subtraction");

        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][j] - other(i, j);
        return result;
    }

    Matrix4x4 operator-(real_t scalar) const
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][j] - scalar;
        return result;
    }

    //[*]

    Matrix4x4 operator*(const Matrix4x4 &other) const
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][0] * other(0, j) +
                               m_data[i][1] * other(1, j) +
                               m_data[i][2] * other(2, j) +
                               m_data[i][3] * other(3, j);
        return result;
    }

    Matrix4x4 operator*(real_t scalar) const
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][j] * scalar;
        return result;
    }

    //[/]

    Matrix4x4 operator/(real_t scalar) const
    {
        Matrix4x4 result;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result(i, j) = m_data[i][j] / scalar;
        return result;
    }
};


} // namespace Math

} // namespace Arns