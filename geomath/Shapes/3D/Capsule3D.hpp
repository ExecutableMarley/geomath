/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "BBox3D.hpp"
#include "Variants/EShapeType3D.hpp"
#include "Variants/ShapeConcepts3D.hpp"
#include "Algorithms3D/Algorithms3D.hpp"

namespace Arns
{

namespace geomath
{

class Capsule3D
{
public:
    Vector3D m_startPoint;
    Vector3D m_endPoint;
    real_t m_radius;

    Capsule3D() : m_startPoint(), m_endPoint(), m_radius(0) {}
    
    Capsule3D(const Vector3D &startPoint, const Vector3D &endPoint, real_t radius) : m_startPoint(startPoint), m_endPoint(endPoint), m_radius(radius) {}

    ShapeType3D type() const 
    {
        return ShapeType3D::SHAPE3D_CAPSULE;
    }

    real_t volume() const 
    {
        return PI * m_radius * m_radius * (4.0f / 3.0f * m_radius + (m_endPoint - m_startPoint).length());
    }

    real_t surfaceArea() const 
    {
        return 2.0f * PI * m_radius * (2.0f * m_radius + (m_endPoint - m_startPoint).length());
    }

    Vector3D centroid() const 
    {
        return (m_startPoint + m_endPoint) * 0.5f;
    }

    Capsule3D& translate(const Vector3D &translation) 
    {
        m_startPoint += translation;
        m_endPoint += translation;
        return *this;
    }

    Capsule3D& scale(real_t scaleFactor)
    {
        Vector3D centroid = this->centroid();
        m_startPoint = centroid + scaleFactor * (m_startPoint - centroid);
        m_endPoint   = centroid + scaleFactor * (m_endPoint - centroid);
        m_radius    *= scaleFactor;
        return *this;
    }

    Capsule3D& rotate(Radians xAngle, Radians yAngle, Radians zAngle)
    {
        Vector3D centroid = this->centroid();
        m_startPoint = m_startPoint.rotateAround(xAngle, yAngle, zAngle, centroid);
        m_endPoint   = m_endPoint.rotateAround(xAngle, yAngle, zAngle, centroid);
        return *this;
    }

    /*
    Capsule& rotate(real_t angle, const Vector3D &axis)
    {
        m_startPoint = m_startPoint.rotate(angle, axis);
        m_endPoint   = m_endPoint.rotate(angle, axis);
        return *this;
    }*/

    /*
    bool contains(const Vector3D &point) const 
    {
        return false;
    }*/

    
    bool contains(const Vector3D &point) const 
    {
        return detail::distancePointToSegment(m_startPoint, m_endPoint, point) < m_radius;
    }
    

    BBox3D boundingBox() const 
    {
        const Vector3D min = Vector3D::min(m_startPoint, m_endPoint) - Vector3D(m_radius, m_radius, m_radius);
        const Vector3D max = Vector3D::max(m_startPoint, m_endPoint) + Vector3D(m_radius, m_radius, m_radius);
        return BBox3D(min, max);
    }
};

static_assert(Translatable3D<Capsule3D>);
static_assert(UniformlyScalable3D<Capsule3D>);
//static_assert(Rotatable3D<Capsule3D>);
//static_assert(Transformable3D<Capsule3D>);
static_assert(PointContainable3D<Capsule3D>);

} // namespace Math

} // namespace Arns