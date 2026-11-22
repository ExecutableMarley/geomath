/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "BBox3D.hpp"
#include "Variants/EShapeType3D.hpp"
#include "Variants/ShapeConcepts3D.hpp"

namespace Arns
{

namespace geomath
{

class Cylinder3D
{
public:
    Vector3D m_startPoint;
    Vector3D m_endPoint;
    real_t m_radius;

    Cylinder3D() : m_startPoint(), m_endPoint(), m_radius(0) {}

    Cylinder3D(const Vector3D &startPoint, const Vector3D &endPoint, real_t radius) : m_startPoint(startPoint), m_endPoint(endPoint), m_radius(radius) {}

    ShapeType3D type() const
    {
        return ShapeType3D::SHAPE3D_CYLINDER;
    }

    real_t volume() const
    {
        return PI * m_radius * m_radius * (m_endPoint - m_startPoint).length();
    }

    real_t surfaceArea() const
    {
        return 2.0 * PI * m_radius * (m_endPoint - m_startPoint).length() + 2.0 * PI * m_radius * m_radius;
    }

    Vector3D centroid() const
    {
        return (m_startPoint + m_endPoint) * 0.5;
    }

    bool contains(const Vector3D &point) const
    {
        // Todo: Calculating the closest parameter t should be done in a separate function
        const Vector3D delta = m_endPoint - m_startPoint;
        real_t t = (point - m_startPoint).dot(delta) / (delta).lengthSquared();

        Vector3D closestPoint = m_startPoint + delta * t;
        if (t >= 0.0 && t <= 1.0 && (point - closestPoint).lengthSquared() <= m_radius * m_radius)
        {
            return true;
        }
        return false;
    }

    Cylinder3D& translate(const Vector3D &translation)
    {
        m_startPoint += translation;
        m_endPoint += translation;
        return *this;
    }

    Cylinder3D& scale(real_t scaleFactor)
    {
        Vector3D centroid = this->centroid();
        m_startPoint = centroid + scaleFactor * (m_startPoint - centroid);
        m_endPoint   = centroid + scaleFactor * (m_endPoint - centroid);
        m_radius    *= scaleFactor;
        return *this;
    }

    Cylinder3D& rotateAround(const Matrix3x3 &rotationMatrix, const Vector3D &point)
    {
        m_startPoint = rotationMatrix * (m_startPoint - point) + point;
        m_endPoint   = rotationMatrix * (m_endPoint - point) + point;
        return *this;
    }

    Cylinder3D& rotate(const Matrix3x3 &rotationMatrix)
    {
        return this->rotateAround(rotationMatrix, this->centroid());
    }

    BBox3D boundingBox() const
    {
        const Vector3D min = Vector3D::min(m_startPoint, m_endPoint) - Vector3D(m_radius, m_radius, m_radius);
        const Vector3D max = Vector3D::max(m_startPoint, m_endPoint) + Vector3D(m_radius, m_radius, m_radius);
        return BBox3D(min, max);
    }
};

static_assert(Translatable3D<Cylinder3D>);
static_assert(UniformlyScalable3D<Cylinder3D>);
static_assert(Rotatable3D<Cylinder3D>);
//static_assert(Transformable3D<Cylinder3D>);
static_assert(PointContainable3D<Cylinder3D>);
static_assert(Bounded3D<Cylinder3D>);

} // namespace Math

} // namespace Arns