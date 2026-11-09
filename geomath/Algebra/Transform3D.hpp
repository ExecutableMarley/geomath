/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "Algebra/EulerAngles.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"

namespace Arns::geomath
{

class Transform3D
{
public:
	Transform3D()
		: m_position(0, 0, 0), m_rotation(), m_scale(1), m_dirty(true) {}

	Transform3D(const Vector3D& position, const EulerAngles& rotation = EulerAngles(), const real_t scale = 1)
		: m_position(position), m_rotation(rotation), m_scale(scale), m_dirty(true) {}

    // --- Setters (Mark cached matrices dirty) ---
        
	void setPosition(const Vector3D& position) { m_position = position; m_dirty = true; }
	void setRotation(const EulerAngles& rotation) { m_rotation = rotation; m_dirty = true; }
	void setScale(const real_t scale) { m_scale = scale; m_dirty = true; }

    // --- Getters ---

	const Vector3D& getPosition() const { return m_position; }
	const EulerAngles& getRotation() const { return m_rotation; }
	const real_t getScale() const { return m_scale; }

    // --- Point and Vector Transformations ---

	Vector3D transformPoint(const Vector3D& point) const
	{
		updateMatrices();
		return m_matrix.transformPoint(point);
	}

	Vector3D transformVector(const Vector3D& vector) const
	{
		updateMatrices();
		return m_matrix.transformVector(vector);
	}

	Vector3D inverseTransformPoint(const Vector3D& worldPoint) const
	{
		updateMatrices();
		return m_invMatrix.transformPoint(worldPoint);
	}

	Vector3D inverseTransformVector(const Vector3D& worldVector) const
	{
		updateMatrices();
		return m_invMatrix.transformVector(worldVector);
	}

    // --- Transform / modification ---

	Transform3D& reset()
	{
		m_position = Vector3D(0, 0, 0);
		m_rotation = EulerAngles();
		m_scale = 1;
		m_dirty = true;
		return *this;
	}

	Transform3D& translate(const Vector3D& translation)
	{
		m_position += translation;
		m_dirty = true;
		return *this;
	}

	Transform3D& rotate(const EulerAngles& rotation)
	{
		m_rotation += rotation;
		m_dirty = true;
		return *this;
	}

	Transform3D& rotate(Radians pitch, Radians yaw, Radians roll)
	{
		return rotate(EulerAngles(pitch, yaw, roll));
	}

	Transform3D& scale(real_t scaleFactor)
	{
		m_scale *= scaleFactor;
        m_dirty  = true;
        return *this;
	}

	Transform3D operator*(const Transform3D& rhs) const
	{
		return Transform3D(
			transformPoint(rhs.m_position),
			m_rotation + rhs.m_rotation,
			m_scale * rhs.m_scale);
	}

	Transform3D& operator*=(const Transform3D& rhs)
	{
		return *this = *this * rhs;
	}

    // Access underlying matrices

	const Matrix4x4& getMatrix() const { updateMatrices(); return m_matrix; }
	const Matrix4x4& getInverseMatrix() const { updateMatrices(); return m_invMatrix; }

    // --- Constants ---

	static Transform3D identity() { return Transform3D(); }

private:
	static real_t inverseScale(real_t value)
	{
		return approximatelyZero(value) ? real_t{0} : real_t{1} / value;
	}

	void updateMatrices() const
	{
		if (!m_dirty) return;

		const Matrix4x4 translation = Matrix4x4::createTranslation3D(m_position.x, m_position.y, m_position.z);
		const Matrix4x4 rotation = m_rotation.toMatrix4x4();
		const Matrix4x4 scaleMatrix = Matrix4x4::createScale3D(m_scale, m_scale, m_scale);
		m_matrix = translation * rotation * scaleMatrix;

		const Matrix4x4 inverseRotation = rotation.transpose();
		const Matrix4x4 inverseScaleMatrix = Matrix4x4::createScale3D(
			inverseScale(m_scale), inverseScale(m_scale), inverseScale(m_scale));
		const Matrix4x4 inverseTranslation = Matrix4x4::createTranslation3D(-m_position.x, -m_position.y, -m_position.z);
		m_invMatrix = inverseScaleMatrix * inverseRotation * inverseTranslation;
		m_dirty = false;
	}

	Vector3D m_position;
	EulerAngles m_rotation;
	real_t m_scale;

	mutable Matrix4x4 m_matrix;
	mutable Matrix4x4 m_invMatrix;
	mutable bool m_dirty = true;
};

}
