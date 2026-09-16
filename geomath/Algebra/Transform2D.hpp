#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"

namespace Arns::geomath
{


class Transform2D 
{
public:
    Transform2D()
        : m_position(0, 0), m_rotationRads(0), m_scale(1), m_dirty(true) {}

    Transform2D(const Vector2D& position, real_t rotationRads = 0, real_t scale = 1)
        : m_position(position), m_rotationRads(rotationRads), m_scale(scale), m_dirty(true) {}

    // --- Setters (Mark cached matrices dirty) ---

    void setPosition(const Vector2D& position) { m_position = position; m_dirty = true; }
    void setRotation(real_t radians) { m_rotationRads = radians; m_dirty = true; }
    void setScale(real_t scale) { m_scale = scale; m_dirty = true; }

    // --- Getters ---

    const Vector2D& getPosition() const { return m_position; }
    real_t getRotation() const { return m_rotationRads; }
    real_t getScale() const { return m_scale; }

    // --- Point and Vector Transformations ---
    
    // Transforms a point from local to world space
    Vector2D transformPoint(const Vector2D& point) const 
    {
        updateMatrices();
        return this->m_matrix.transformPoint(point);
    }

    // Transforms a direction/vector
    Vector2D transformVector(const Vector2D& vec) const 
    {
        updateMatrices();
        return this->m_matrix.transformVector(vec);
    }

    // Transforms a point from world space back to local mesh space
    Vector2D inverseTransformPoint(const Vector2D& worldPoint) const 
    {
        updateMatrices();
        return this->m_invMatrix.transformPoint(worldPoint);
    }

    // Transforms a vector from world space back to local mesh space
    Vector2D inverseTransformVector(const Vector2D& worldVec) const 
    {
        updateMatrices();
        return this->m_invMatrix.transformVector(worldVec);
    }

    // --- Transform / modification ---

    Transform2D& reset()
    {
        this->m_position     = Vector2D::zero();
        this->m_rotationRads = real_t{0};
        this->m_scale        = real_t{1};
        m_dirty = true;
        return *this;
    }

    Transform2D& translate(const Vector2D &translation) 
    {
        m_position += translation; 
        m_dirty = true; 
        return *this; 
    }

    Transform2D& rotate(real_t angle) 
    {
        m_rotationRads += angle; 
        m_dirty = true; 
        return *this; 
    }

    Transform2D& scale(real_t scaleFactor) 
    {
        m_scale *= scaleFactor;
        m_dirty = true; 
        return *this; 
    }

    Transform2D operator*(const Transform2D& rhs) const
    {
        return Transform2D(
            transformPoint(rhs.m_position),
            m_rotationRads + rhs.m_rotationRads,
            m_scale * rhs.m_scale);
    }

    Transform2D& operator *=(const Transform2D& rhs)
    {
        return (*this) = (*this) * rhs;
    }

    // Access underlying matrices

    const Matrix3x3& getMatrix()        const { updateMatrices(); return m_matrix; }
    const Matrix3x3& getInverseMatrix() const { updateMatrices(); return m_invMatrix; }

    // --- Constants ---

    static Transform2D identity() { return Transform2D(); }

private:
    void updateMatrices() const 
    {
        if (!m_dirty) return;

        real_t sinAngle, cosAngle;
        sinCos(m_rotationRads, sinAngle, cosAngle);

        // Pre-calculate rotation * scale terms
        real_t m00 = cosAngle * m_scale;
        real_t m01 = -sinAngle * m_scale;
        real_t m10 = sinAngle * m_scale;
        real_t m11 = cosAngle * m_scale;

        // Forward Matrix (SRT: Scale, Rotate, Translate)
        m_matrix = Matrix3x3(
            m00, m01, m_position.x,
            m10, m11, m_position.y,
            0, 0, 1
        );

        // Fast Inverse (T^-1 * R^-1 * S^-1) computed analytically skipping 3x3 inversion
        real_t invScale = (m_scale != real_t{0}) ? (real_t{1} / m_scale) : real_t{0};

        real_t inv00 = cosAngle * invScale;
        real_t inv01 = sinAngle * invScale;
        real_t inv10 = -sinAngle * invScale;
        real_t inv11 = cosAngle * invScale;

        // Apply inverse translation: -R^-1 * S^-1 * T
        real_t inv02 = -(inv00 * m_position.x + inv01 * m_position.y);
        real_t inv12 = -(inv10 * m_position.x + inv11 * m_position.y);

        m_invMatrix = Matrix3x3(
            inv00, inv01, inv02,
            inv10, inv11, inv12,
            0,  0,  1
        );

        m_dirty = false;
    }

    Vector2D m_position;
    real_t m_rotationRads;
    real_t m_scale;

    mutable Matrix3x3 m_matrix;
    mutable Matrix3x3 m_invMatrix;
    mutable bool m_dirty = true;
};


}