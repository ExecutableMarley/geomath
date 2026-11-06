/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "IMatrix.hpp"
#include "Geometry/Vector3D.hpp"
#include "Geometry/Vector4D.hpp"
#include <stdexcept>

namespace Arns
{

namespace geomath
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

    real_t determinant() const
    {
        const real_t a00 = m_data[0][0];
        const real_t a01 = m_data[0][1];
        const real_t a02 = m_data[0][2];
        const real_t a03 = m_data[0][3];

        const real_t a10 = m_data[1][0];
        const real_t a11 = m_data[1][1];
        const real_t a12 = m_data[1][2];
        const real_t a13 = m_data[1][3];

        const real_t a20 = m_data[2][0];
        const real_t a21 = m_data[2][1];
        const real_t a22 = m_data[2][2];
        const real_t a23 = m_data[2][3];

        const real_t a30 = m_data[3][0];
        const real_t a31 = m_data[3][1];
        const real_t a32 = m_data[3][2];
        const real_t a33 = m_data[3][3];

        const real_t b00 = a00 * a11 - a01 * a10;
        const real_t b01 = a00 * a12 - a02 * a10;
        const real_t b02 = a00 * a13 - a03 * a10;
        const real_t b03 = a01 * a12 - a02 * a11;
        const real_t b04 = a01 * a13 - a03 * a11;
        const real_t b05 = a02 * a13 - a03 * a12;

        const real_t b06 = a20 * a31 - a21 * a30;
        const real_t b07 = a20 * a32 - a22 * a30;
        const real_t b08 = a20 * a33 - a23 * a30;
        const real_t b09 = a21 * a32 - a22 * a31;
        const real_t b10 = a21 * a33 - a23 * a31;
        const real_t b11 = a22 * a33 - a23 * a32;

        return b00 * b11 -
               b01 * b10 +
               b02 * b09 +
               b03 * b08 -
               b04 * b07 +
               b05 * b06;
    }

    bool isInvertible() const
    {
        return !approximatelyZero(determinant());
    }

    bool inverse(Matrix4x4 &out) const
    {
        const real_t a00 = m_data[0][0];
        const real_t a01 = m_data[0][1];
        const real_t a02 = m_data[0][2];
        const real_t a03 = m_data[0][3];

        const real_t a10 = m_data[1][0];
        const real_t a11 = m_data[1][1];
        const real_t a12 = m_data[1][2];
        const real_t a13 = m_data[1][3];

        const real_t a20 = m_data[2][0];
        const real_t a21 = m_data[2][1];
        const real_t a22 = m_data[2][2];
        const real_t a23 = m_data[2][3];

        const real_t a30 = m_data[3][0];
        const real_t a31 = m_data[3][1];
        const real_t a32 = m_data[3][2];
        const real_t a33 = m_data[3][3];

        const real_t b00 = a00 * a11 - a01 * a10;
        const real_t b01 = a00 * a12 - a02 * a10;
        const real_t b02 = a00 * a13 - a03 * a10;
        const real_t b03 = a01 * a12 - a02 * a11;
        const real_t b04 = a01 * a13 - a03 * a11;
        const real_t b05 = a02 * a13 - a03 * a12;

        const real_t b06 = a20 * a31 - a21 * a30;
        const real_t b07 = a20 * a32 - a22 * a30;
        const real_t b08 = a20 * a33 - a23 * a30;
        const real_t b09 = a21 * a32 - a22 * a31;
        const real_t b10 = a21 * a33 - a23 * a31;
        const real_t b11 = a22 * a33 - a23 * a32;

        const real_t det =
            b00 * b11 -
            b01 * b10 +
            b02 * b09 +
            b03 * b08 -
            b04 * b07 +
            b05 * b06;

        if (approximatelyZero(det))
            return false;

        const real_t invDet = real_t{1} / det;

        out(0, 0) = (a11 * b11 - a12 * b10 + a13 * b09)  * invDet;
        out(0, 1) = (-a01 * b11 + a02 * b10 - a03 * b09) * invDet;
        out(0, 2) = (a31 * b05 - a32 * b04 + a33 * b03)  * invDet;
        out(0, 3) = (-a21 * b05 + a22 * b04 - a23 * b03) * invDet;

        out(1, 0) = (-a10 * b11 + a12 * b08 - a13 * b07) * invDet;
        out(1, 1) = (a00 * b11 - a02 * b08 + a03 * b07)  * invDet;
        out(1, 2) = (-a30 * b05 + a32 * b02 - a33 * b01) * invDet;
        out(1, 3) = (a20 * b05 - a22 * b02 + a23 * b01)  * invDet;

        out(2, 0) = (a10 * b10 - a11 * b08 + a13 * b06)  * invDet;
        out(2, 1) = (-a00 * b10 + a01 * b08 - a03 * b06) * invDet;
        out(2, 2) = (a30 * b04 - a31 * b02 + a33 * b00)  * invDet;
        out(2, 3) = (-a20 * b04 + a21 * b02 - a23 * b00) * invDet;

        out(3, 0) = (-a10 * b09 + a11 * b07 - a12 * b06) * invDet;
        out(3, 1) = (a00 * b09 - a01 * b07 + a02 * b06)  * invDet;
        out(3, 2) = (-a30 * b03 + a31 * b01 - a32 * b00) * invDet;
        out(3, 3) = (a20 * b03 - a21 * b01 + a22 * b00)  * invDet;

        return true;
    }

    //

    Vector3D transformPoint(const Vector3D& vec) const
    {
        real_t resultX = m_data[0][0] * vec.x + m_data[0][1] * vec.y + m_data[0][2] * vec.z + m_data[0][3];
        real_t resultY = m_data[1][0] * vec.x + m_data[1][1] * vec.y + m_data[1][2] * vec.z + m_data[1][3];
        real_t resultZ = m_data[2][0] * vec.x + m_data[2][1] * vec.y + m_data[2][2] * vec.z + m_data[2][3];
        real_t resultW = m_data[3][0] * vec.x + m_data[3][1] * vec.y + m_data[3][2] * vec.z + m_data[3][3];

        // Normalization
        if (!approximatelyZero(resultW) && !approximatelyEqual(resultW, real_t{1}))
            return Vector3D(resultX / resultW, resultY / resultW, resultZ / resultW);
        return Vector3D(resultX, resultY, resultZ);
    }

    Vector3D transformVector(const Vector3D& vec) const
    {
        return Vector3D(
            m_data[0][0] * vec.x + m_data[0][1] * vec.y + m_data[0][2] * vec.z,
            m_data[1][0] * vec.x + m_data[1][1] * vec.y + m_data[1][2] * vec.z,
            m_data[2][0] * vec.x + m_data[2][1] * vec.y + m_data[2][2] * vec.z);
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

    Vector4D operator*(const Vector4D& vec) const
    {
        return Vector4D(
            m_data[0][0] * vec.x + m_data[0][1] * vec.y + m_data[0][2] * vec.z + m_data[0][3] * vec.w,
            m_data[1][0] * vec.x + m_data[1][1] * vec.y + m_data[1][2] * vec.z + m_data[1][3] * vec.w,
            m_data[2][0] * vec.x + m_data[2][1] * vec.y + m_data[2][2] * vec.z + m_data[2][3] * vec.w,
            m_data[3][0] * vec.x + m_data[3][1] * vec.y + m_data[3][2] * vec.z + m_data[3][3] * vec.w);
    }

    Vector3D operator*(const Vector3D& vec) const
    {
        return Vector3D(
            m_data[0][0] * vec.x + m_data[0][1] * vec.y + m_data[0][2] * vec.z,
            m_data[1][0] * vec.x + m_data[1][1] * vec.y + m_data[1][2] * vec.z,
            m_data[2][0] * vec.x + m_data[2][1] * vec.y + m_data[2][2] * vec.z);
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

    //[Static]  

    static Matrix4x4 identity()
    {
        return Matrix4x4(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1);
    }

    // 3D homogeneous transforms

    static Matrix4x4 createTranslation3D(real_t tx, real_t ty, real_t tz)
    {
        return Matrix4x4(
            1.0, 0.0, 0.0, tx,
            0.0, 1.0, 0.0, ty,
            0.0, 0.0, 1.0, tz,
            0.0, 0.0, 0.0, 1.0);
    }

    static Matrix4x4 createScale3D(real_t sx, real_t sy, real_t sz)
    {
        return Matrix4x4(
            sx , 0.0, 0.0, 0.0,
            0.0, sy , 0.0, 0.0,
            0.0, 0.0, sz , 0.0,
            0.0, 0.0, 0.0, 1.0
        );
    }

    static Matrix4x4 createRotation3D(Radians angle, const Vector3D& axis)
    {
        real_t sin, cos;
        sinCos(angle.value, sin, cos);
        real_t oneMinusCos = 1.0 - cos;

        return Matrix4x4(
            cos + axis.x * axis.x * oneMinusCos,
            axis.x * axis.y * oneMinusCos - axis.z * sin,
            axis.x * axis.z * oneMinusCos + axis.y * sin,
            0.0,

            axis.y * axis.x * oneMinusCos + axis.z * sin,
            cos + axis.y * axis.y * oneMinusCos,
            axis.y * axis.z * oneMinusCos - axis.x * sin,
            0.0,

            axis.z * axis.x * oneMinusCos - axis.y * sin,
            axis.z * axis.y * oneMinusCos + axis.x * sin,
            cos + axis.z * axis.z * oneMinusCos,
            0.0,

            0.0, 0.0, 0.0, 1.0
        );
    }

    static Matrix4x4 createRotationX(Radians angle)
    {
        real_t sin, cos;
        sinCos(angle.value, sin, cos);

        return Matrix4x4(
            1.0, 0.0,  0.0, 0.0,
            0.0, cos, -sin, 0.0,
            0.0, sin,  cos, 0.0,
            0.0, 0.0,  0.0, 1.0
        );
    }

    static Matrix4x4 createRotationY(Radians angle)
    {
        real_t sin, cos;
        sinCos(angle.value, sin, cos);

        return Matrix4x4(
            cos , 0.0, sin, 0.0,
            0.0 , 1.0, 0.0, 0.0,
            -sin, 0.0, cos, 0.0,
            0.0 , 0.0, 0.0, 1.0
        );
    }

    static Matrix4x4 createRotationZ(Radians angle)
    {
        real_t sin, cos;
        sinCos(angle.value, sin, cos);

        return Matrix4x4(
            cos, -sin, 0.0, 0.0,
            sin,  cos, 0.0, 0.0,
            0.0,  0.0, 1.0, 0.0,
            0.0,  0.0, 0.0, 1.0
        );
    }
};

} // namespace Math

} // namespace Arns