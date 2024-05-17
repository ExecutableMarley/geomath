/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Vector2D.hpp"
#include "Vector3D.hpp"
#include "Vector4D.hpp"

#include "Matrices/Matrix4x4.hpp"

namespace Utility
{

namespace Math
{

class ViewMatrix : public Matrix4x4
{
public:
    using Matrix4x4::Matrix4x4;

    ViewMatrix()
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i == j)
                    m_data[i][j] = 1;
                else
                    m_data[i][j] = 0;
    }

    Vector4D transform(const Vector4D &vec4) const
    {
        return Vector4D(
            vec4.x * m_data[0][0] + vec4.y * m_data[0][1] + vec4.z * m_data[0][2] + vec4.w * m_data[0][3],
            vec4.x * m_data[1][0] + vec4.y * m_data[1][1] + vec4.z * m_data[1][2] + vec4.w * m_data[1][3],
            vec4.x * m_data[2][0] + vec4.y * m_data[2][1] + vec4.z * m_data[2][2] + vec4.w * m_data[2][3],
            vec4.x * m_data[3][0] + vec4.y * m_data[3][1] + vec4.z * m_data[3][2] + vec4.w * m_data[3][3]);
    }

    Vector4D transform(const Vector3D &vec3) const
    {
        return Vector4D(
            vec3.x * m_data[0][0] + vec3.y * m_data[0][1] + vec3.z * m_data[0][2] + m_data[0][3],
            vec3.x * m_data[1][0] + vec3.y * m_data[1][1] + vec3.z * m_data[1][2] + m_data[1][3],
            vec3.x * m_data[2][0] + vec3.y * m_data[2][1] + vec3.z * m_data[2][2] + m_data[2][3],
            vec3.x * m_data[3][0] + vec3.y * m_data[3][1] + vec3.z * m_data[3][2] + m_data[3][3]);
    }

    bool isInFront(const Vector3D &point) const
    {
        return point.x * m_data[0][2] + point.y * m_data[1][2] + point.z * m_data[2][2] + m_data[3][2] > 0;
    }

    bool worldToScreen(const Vector3D &world, Vector3D &screen, int screenWidth, int screenHeight) const
    {
        const Vector4D clip = transform(world);
        if (clip.w <= 0)
            return false;
        screen.x = (clip.x / clip.w + 1.0f) * 0.5f * screenWidth;
        screen.y = (1.0f - clip.y / clip.w) * 0.5f * screenHeight;
        screen.z = clip.z / clip.w;
        return true;
    }

    bool worldToScreen(const Vector3D &world, Vector2D &screen, int screenWidth, int screenHeight) const
    {
        const Vector4D clip = transform(world);
        if (clip.w <= 0)
            return false;
        screen.x = (clip.x / clip.w + 1.0f) * 0.5f * screenWidth;
        screen.y = (1.0f - clip.y / clip.w) * 0.5f * screenHeight;
        return true;
    }
};

} // namespace Math

} // namespace Utility