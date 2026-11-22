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

class Sphere3D
{
public:
    Vector3D m_center;
    real_t m_radius;

    Sphere3D() : m_center(), m_radius(0) {}

    Sphere3D(const Vector3D &center, real_t radius) : m_center(center), m_radius(radius) {}

    ShapeType3D type() const 
    {
        return ShapeType3D::SHAPE3D_SPHERE;
    }

    real_t volume() const 
    {
        return 4.0 / 3.0 * PI * m_radius * m_radius * m_radius;
    }

    real_t surfaceArea() const 
    {
        return 4.0 * PI * m_radius * m_radius;
    }

    Vector3D centroid() const 
    {
        return m_center;
    }

    bool contains(const Vector3D &point) const 
    {
        return (point - m_center).lengthSquared() <= m_radius * m_radius;
    }

    bool contains(const Sphere3D &sphere) const
    {
        return (m_center - sphere.m_center).lengthSquared() + sphere.m_radius * sphere.m_radius <= m_radius * m_radius;
    }

    Sphere3D& translate(const Vector3D &translation) 
    {
        m_center += translation;
        return *this;
    }

    Sphere3D& scale(real_t scaleFactor)
    {
        m_radius *= scaleFactor;
        return *this;
    }

    Sphere3D& rotateAround(const Matrix3x3 &rotationMatrix, const Vector3D &point)
    {
        m_center = rotationMatrix * (m_center - point) + point;
        return *this;
    }

    Sphere3D& rotate(const Matrix3x3 &rotationMatrix)
    {
        return *this;
    }

    bool intersects(const Sphere3D &sphere) const
    {
        return (m_center - sphere.m_center).lengthSquared() <= (m_radius + sphere.m_radius) * (m_radius + sphere.m_radius);
    }

    BBox3D boundingBox() const 
    {
        return BBox3D(m_center - Vector3D(m_radius, m_radius, m_radius), m_center + Vector3D(m_radius, m_radius, m_radius));
    }
};

static_assert(Translatable3D<Sphere3D>);
static_assert(UniformlyScalable3D<Sphere3D>);
static_assert(Rotatable3D<Sphere3D>);
//static_assert(Transformable3D<Sphere3D>);
static_assert(PointContainable3D<Sphere3D>);
static_assert(Bounded3D<Sphere3D>);

} // namespace Math

} // namespace Arns