/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <array>
#include <vector>
#include <memory>
#include <stdexcept>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Shapes/2D/Algorithms2D/Algorithms2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Arns
{

namespace Math
{

class Triangle2D : public IPolygonalShape2D
{
public:
    std::array<Vector2D, 3> m_vertices;

    // --- Constructors ---

    Triangle2D() :  m_vertices{ Vector2D{}, Vector2D{}, Vector2D{} } {}

    Triangle2D(const Vector2D &a, const Vector2D &b, const Vector2D &c) : m_vertices{ a, b, c } {}

    ShapeType2D type() const override { return SHAPE2D_TRIANGLE; }
    static constexpr ShapeType2D shapeType = SHAPE2D_TRIANGLE;

    Vector2D& a() { return m_vertices[0]; }
    Vector2D& b() { return m_vertices[1]; }
    Vector2D& c() { return m_vertices[2]; }

    const Vector2D& a() const { return m_vertices[0]; }
    const Vector2D& b() const { return m_vertices[1]; }
    const Vector2D& c() const { return m_vertices[2]; }

    size_t vertexCount() const override
    {
        return 3;
    }

    const Vector2D& vertexAt(size_t index) const
    {
        return m_vertices[index];
    }

    Vector2D& vertexAt(size_t index)
    {
        return m_vertices[index];
    }

    const std::span<const Vector2D> vertices() const override
    {
        return m_vertices;
    }

    real_t area() const override
    {
        return 0.5f * fabs((a().x - c().x) * (b().y - c().y) - (b().x - c().x) * (a().y - c().y));
    }

    real_t perimeter() const override
    {
        return (a() - b()).length() + (b() - c()).length() + (c() - a()).length();
    }

    Vector2D centroid() const override
    {
        return (a() + b() + c()) / real_t{3};
    }

    // --- Transform / modification ---

    Triangle2D& translate(const Vector2D &translation) override
    {
        for (auto& v : m_vertices)
            v += translation;
        return *this;
    }

    Triangle2D& rotate(real_t angle)
    {
        const Vector2D c = centroid();
        for (auto& v : m_vertices)
            v.rotateAround(angle, c);
        return *this;
    }

    Triangle2D& rotate(real_t angle, const Vector2D& point)
    {
        for (auto& v : m_vertices)
            v.rotateAround(angle, point);
        return *this;
    }

    Triangle2D& scale(real_t factor)
    {
        const Vector2D c = centroid();
        for (auto& v : m_vertices)
            v = c + (v - c) * factor;
        return *this;
    }

    Triangle2D copy() const 
    {
        return *this;
    }

    std::unique_ptr<IFiniteShape2D> clone() const override
    {
        return std::make_unique<Triangle2D>(*this);        
    }

    bool contains(const Vector2D &point) const override
    {
        const real_t areaABC = area();
        const auto& m_a = this->a();
        const auto& m_b = this->b();
        const auto& m_c = this->c();
        const real_t areaPBC = real_t{0.5} * std::abs((m_b.x - point.x) * (m_c.y - point.y) - (m_c.x - point.x) * (m_b.y - point.y));
        const real_t areaPCA = real_t{0.5} * std::abs((m_c.x - point.x) * (m_a.y - point.y) - (m_a.x - point.x) * (m_c.y - point.y));
        const real_t areaPAB = real_t{0.5} * std::abs((m_a.x - point.x) * (m_b.y - point.y) - (m_b.x - point.x) * (m_a.y - point.y));
        return approximatelyEqual(areaABC, areaPBC + areaPCA + areaPAB);
    }

    bool contains(const Segment2D& segment) const
    {
        return isSegmentInsideShape(segment, *this);
    }

    bool contains(const Triangle2D &triangle) const
    {
        return contains(triangle.a()) && contains(triangle.b()) && contains(triangle.c());
    }

    bool contains(const IFiniteShape2D& shape) const
    {
        return isShapeInsideConvexPolygon(shape, m_vertices);
    }

    BBox2D boundingBox() const override
    {
        return BBox2D(Vector2D::min(a(), b(), c()), Vector2D::max(a(), b(), c()));
    }

    Vector2D& operator[](size_t index)
    {
        return vertexAt(index);
    }

    const Vector2D& operator[](size_t index) const override
    {
        return vertexAt(index);
    }

    // --- Operators ---

    bool operator==(const Triangle2D& other) const
    {
        return this->a() == other.a() && 
            this->b() == other.b() && 
            this->c() == other.c();
    }

    bool operator!=(const Triangle2D& other) const
    {
        return !(*this == other);
    }
};

} // namespace Math

} // namespace Arns