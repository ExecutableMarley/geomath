/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "../CommonMath.hpp"
#include "Geometry/Vector3D.hpp"

#include "Matrices/Matrix3x3.hpp"

namespace Utility
{

namespace Math
{

class RotationMatrix : public Matrix3x3
{
public:
    using Matrix3x3::Matrix3x3;

    RotationMatrix()
    {
        m_data[0][0] = 1.0f;
        m_data[0][1] = 0.0f;
        m_data[0][2] = 0.0f;

        m_data[1][0] = 0.0f;
        m_data[1][1] = 1.0f;
        m_data[1][2] = 0.0f;

        m_data[2][0] = 0.0f;
        m_data[2][1] = 0.0f;
        m_data[2][2] = 1.0f;
    }

    Vector3D transform(Vector3D &vec3) const
    {
        return Vector3D(
            vec3.x * m_data[0][0] + vec3.y * m_data[0][1] + vec3.z * m_data[0][2],
            vec3.x * m_data[1][0] + vec3.y * m_data[1][1] + vec3.z * m_data[1][2],
            vec3.x * m_data[2][0] + vec3.y * m_data[2][1] + vec3.z * m_data[2][2]);
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