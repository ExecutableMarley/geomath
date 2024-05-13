/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "../CommonMath.hpp"
#include "Geometry/Vector3D.hpp"

namespace Utility
{

namespace Math
{
//Matrix3x3 Specialization
class RotationMatrix
{
    float m[3][3];

    RotationMatrix()
    {
        m[0][0] = 1.0f;
        m[0][1] = 0.0f;
        m[0][2] = 0.0f;

        m[1][0] = 0.0f;
        m[1][1] = 1.0f;
        m[1][2] = 0.0f;

        m[2][0] = 0.0f;
        m[2][1] = 0.0f;
        m[2][2] = 1.0f;
    }

    RotationMatrix(float m00, float m01, float m02,
                   float m10, float m11, float m12,
                   float m20, float m21, float m22)
    {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;

        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;

        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
    }

    RotationMatrix(const RotationMatrix &other)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m[i][j] = other.m[i][j];
    }

    Vector3D transform(Vector3D &vec3) const
    {
        return Vector3D(
            vec3.x * m[0][0] + vec3.y * m[0][1] + vec3.z * m[0][2],
            vec3.x * m[1][0] + vec3.y * m[1][1] + vec3.z * m[1][2],
            vec3.x * m[2][0] + vec3.y * m[2][1] + vec3.z * m[2][2]);
    }

    // Static functions
};

//Todo: Implement the following functions:
//Create from Euler Angles
//Create from Quaternion
//Rotate Around Axis by Angle
//Transform Vector3D


} // namespace Math

} // namespace Utility