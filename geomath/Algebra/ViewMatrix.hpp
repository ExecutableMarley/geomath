/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Geometry/Vector2D.hpp"
#include "Geometry/Vector3D.hpp"
#include "Geometry/Vector4D.hpp"

#include "Matrices/Matrix3x3.hpp"
#include "Matrices/Matrix4x4.hpp"

namespace Arns
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
        screen.x = (clip.x / clip.w + 1) * real_t{0.5} * screenWidth;
        screen.y = (1 - clip.y / clip.w) * real_t{0.5} * screenHeight;
        screen.z = clip.z / clip.w;
        return true;
    }

    bool worldToScreen(const Vector3D &world, Vector2D &screen, int screenWidth, int screenHeight) const
    {
        const Vector4D clip = transform(world);
        if (clip.w <= 0)
            return false;
        screen.x = (clip.x / clip.w + 1) * real_t{0.5} * screenWidth;
        screen.y = (1 - clip.y / clip.w) * real_t{0.5} * screenHeight;
        return true;
    }

    static ViewMatrix LookAt(
        const Vector3D &eye,
        const Vector3D &target,
        const Vector3D &up = Vector3D(0, 1, 0))
    {
        const Vector3D forward = (target - eye).normalize();
        const Vector3D right = forward.cross(up).normalize();
        const Vector3D cameraUp = right.cross(forward).normalize();

        ViewMatrix view;

        view.m_data[0][0] = right.x;
        view.m_data[0][1] = right.y;
        view.m_data[0][2] = right.z;
        view.m_data[0][3] = -right.dot(eye);

        view.m_data[1][0] = cameraUp.x;
        view.m_data[1][1] = cameraUp.y;
        view.m_data[1][2] = cameraUp.z;
        view.m_data[1][3] = -cameraUp.dot(eye);

        view.m_data[2][0] = forward.x;
        view.m_data[2][1] = forward.y;
        view.m_data[2][2] = forward.z;
        view.m_data[2][3] = -forward.dot(eye);

        view.m_data[3][0] = 0;
        view.m_data[3][1] = 0;
        view.m_data[3][2] = 0;
        view.m_data[3][3] = 1;

        return view;
    }

    static ViewMatrix Perspective(
        real_t fovRad,
        real_t aspect,
        real_t nearPlane,
        real_t farPlane)
    {
        const real_t f =
            real_t(1) / std::tan(fovRad * real_t(0.5));

        ViewMatrix result{};

        result.m_data[0][0] = f / aspect;
        result.m_data[0][1] = 0;
        result.m_data[0][2] = 0;
        result.m_data[0][3] = 0;

        result.m_data[1][0] = 0;
        result.m_data[1][1] = f;
        result.m_data[1][2] = 0;
        result.m_data[1][3] = 0;

        result.m_data[2][0] = 0;
        result.m_data[2][1] = 0;
        result.m_data[2][2] =
            farPlane / (farPlane - nearPlane);
        result.m_data[2][3] =
            -(nearPlane * farPlane) /
            (farPlane - nearPlane);

        result.m_data[3][0] = 0;
        result.m_data[3][1] = 0;
        result.m_data[3][2] = 1;
        result.m_data[3][3] = 0;

        return result;
    }

    static ViewMatrix Orthographic(
        real_t left,
        real_t right,
        real_t bottom,
        real_t top,
        real_t nearPlane,
        real_t farPlane)
    {
        ViewMatrix result{};

        result.m_data[0][0] =
            real_t(2) / (right - left);

        result.m_data[0][3] =
            -(right + left) /
            (right - left);

        result.m_data[1][1] =
            real_t(2) / (top - bottom);

        result.m_data[1][3] =
            -(top + bottom) /
            (top - bottom);

        result.m_data[2][2] =
            real_t(1) /
            (farPlane - nearPlane);

        result.m_data[2][3] =
            -nearPlane /
            (farPlane - nearPlane);

        result.m_data[3][3] = 1;

        return result;
    }

    static ViewMatrix FromEuler(const Vector3D &eye, real_t pitchRad, real_t yawRad, real_t rollRad = 0)
    {
        // Calculate forward vector from Yaw and Pitch
        Vector3D forward;
        forward.x = cosf(pitchRad) * cosf(yawRad);
        forward.y = sinf(pitchRad);
        forward.z = cosf(pitchRad) * sinf(yawRad);
        forward.normalize();

        return LookAt(eye, eye + forward, Vector3D::unitY());
    }
};

} // namespace Math

} // namespace Arns