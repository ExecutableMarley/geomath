/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Vector2D.hpp"
#include "Vector3D.hpp"
#include "Vector4D.hpp"

namespace Utility
{

namespace Math
{
// Matrix4x4 Specialization
class ViewMatrix
{
    float m[4][4];

    ViewMatrix()
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i == j)
                    m[i][j] = 1;
                else
                    m[i][j] = 0;
    }

    ViewMatrix(float m00, float m01, float m02, float m03,
               float m10, float m11, float m12, float m13,
               float m20, float m21, float m22, float m23,
               float m30, float m31, float m32, float m33)
    {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;
        m[0][3] = m03;
        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;
        m[1][3] = m13;
        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
        m[2][3] = m23;
        m[3][0] = m30;
        m[3][1] = m31;
        m[3][2] = m32;
        m[3][3] = m33;
    }

    ViewMatrix(const ViewMatrix &other)
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                m[i][j] = other.m[i][j];
    }

    Vector4D transform(const Vector4D &vec4) const
    {
        return Vector4D(
            vec4.x * m[0][0] + vec4.y * m[0][1] + vec4.z * m[0][2] + vec4.w * m[0][3],
            vec4.x * m[1][0] + vec4.y * m[1][1] + vec4.z * m[1][2] + vec4.w * m[1][3],
            vec4.x * m[2][0] + vec4.y * m[2][1] + vec4.z * m[2][2] + vec4.w * m[2][3],
            vec4.x * m[3][0] + vec4.y * m[3][1] + vec4.z * m[3][2] + vec4.w * m[3][3]);
    }

    Vector4D transform(const Vector3D &vec3) const
    {
        return Vector4D(
            vec3.x * m[0][0] + vec3.y * m[0][1] + vec3.z * m[0][2] + m[0][3],
            vec3.x * m[1][0] + vec3.y * m[1][1] + vec3.z * m[1][2] + m[1][3],
            vec3.x * m[2][0] + vec3.y * m[2][1] + vec3.z * m[2][2] + m[2][3],
            vec3.x * m[3][0] + vec3.y * m[3][1] + vec3.z * m[3][2] + m[3][3]);
    }

    bool isInFront(const Vector3D &point) const
    {
        return point.x * m[0][2] + point.y * m[1][2] + point.z * m[2][2] + m[3][2] > 0;
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