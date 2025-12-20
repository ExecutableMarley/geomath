/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <stdexcept>

#include "IMatrix.hpp"
#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Geometry/Vector3D.hpp"

namespace Arns
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

    real_t determinant() const
    {
        const real_t a = m_data[0][0], b = m_data[0][1], c = m_data[0][2];
        const real_t d = m_data[1][0], e = m_data[1][1], f = m_data[1][2];
        const real_t g = m_data[2][0], h = m_data[2][1], i = m_data[2][2];

        return a * (e * i - f * h) -
            b * (d * i - f * g) +
            c * (d * h - e * g);
    }

    bool isInvertible() const
    {
        return !approximatelyZero(determinant());
    }

    bool inverse(Matrix3x3& out) const
    {
        const real_t a = m_data[0][0], b = m_data[0][1], c = m_data[0][2];
        const real_t d = m_data[1][0], e = m_data[1][1], f = m_data[1][2];
        const real_t g = m_data[2][0], h = m_data[2][1], i = m_data[2][2];

        real_t det =
            a * (e * i - f * h) -
            b * (d * i - f * g) +
            c * (d * h - e * g);

        if (approximatelyZero(det))
            return false;

        real_t invDet = real_t{1.0} / det;

        out.m_data[0][0] =  (e * i - f * h) * invDet;
        out.m_data[0][1] =  (c * h - b * i) * invDet;
        out.m_data[0][2] =  (b * f - c * e) * invDet;

        out.m_data[1][0] =  (f * g - d * i) * invDet;
        out.m_data[1][1] =  (a * i - c * g) * invDet;
        out.m_data[1][2] =  (c * d - a * f) * invDet;

        out.m_data[2][0] =  (d * h - e * g) * invDet;
        out.m_data[2][1] =  (b * g - a * h) * invDet;
        out.m_data[2][2] =  (a * e - b * d) * invDet;

        return true;
    }

    //

    Vector2D transform(const Vector2D& vec) const
    {
        return Vector2D(
            m_data[0][0] * vec.x + m_data[0][1] * vec.y + m_data[0][2] * 1,
            m_data[1][0] * vec.x + m_data[1][1] * vec.y + m_data[1][2] * 1);
    }

    float* operator[](size_t row)
    {
        return m_data[row];
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

    Vector3D operator*(const Vector3D& other) const
    {
        return Vector3D(
            m_data[0][0] * other.x + m_data[0][1] * other.y + m_data[0][2] * other.z,
            m_data[1][0] * other.x + m_data[1][1] * other.y + m_data[1][2] * other.z,
            m_data[2][0] * other.x + m_data[2][1] * other.y + m_data[2][2] * other.z);
    }

    Vector2D operator*(const Vector2D& other)
    {
        real_t resultX = m_data[0][0] * other.x + m_data[0][1] * other.y + m_data[0][2] * 1.0f;
        real_t resultY = m_data[1][0] * other.x + m_data[1][1] * other.y + m_data[1][2] * 1.0f;
        real_t resultW = m_data[2][0] * other.x + m_data[2][1] * other.y + m_data[2][2] * 1.0f;

        // Normalization
        if (!approximatelyZero(resultW) && !approximatelyEqual(resultW, 1.0f))
            return Vector2D(resultX / resultW, resultY / resultW);
        return Vector2D(resultX, resultY);
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

    //[Static]

    static Matrix3x3 createTranslation(const Vector2D& translation)
    {
        return Matrix3x3(
            1, 0, translation.x,
            0, 1, translation.y,
            0, 0, 1.f);
    }

    static Matrix3x3 createScale(real_t sx, real_t sy)
    {
        return Matrix3x3(
            sx,   0.0f, 0.0f,
            0.0f, sy,   0.0f,
            0.0f, 0.0f, 1.0f
        );
    }

    static Matrix3x3 createRotationRads(real_t angle)
    {
        float sin, cos;
        sinCos(angle, sin, cos);
        
        return Matrix3x3(
            cos, -sin, 0.f,
            sin,  cos, 0.f,
            0.f,  0.f, 1.f);
    }

    static Matrix3x3 createRotationDegs(real_t angle)
    {
        float sin, cos;
        sinCosDeg(angle, sin, cos);
        
        return Matrix3x3(
            cos, -sin, 0.f,
            sin,  cos, 0.f,
            0.f,  0.f, 1.f);
    }
};

} // namespace Math

} // namespace Arns