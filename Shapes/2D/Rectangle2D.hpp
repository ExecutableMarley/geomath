/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>
#include <memory>
#include <stdexcept>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Arns
{

namespace Math
{

class Rectangle2D : public IFiniteShape2D, public IPolygonalShape2D
{
public:
    std::array<Vector2D, 4> m_vertices;

    Rectangle2D() : m_vertices{ Vector2D{}, Vector2D{}, Vector2D{} } {}

    Rectangle2D(const Vector2D &a, const Vector2D &b, const Vector2D &c, const Vector2D &d) : m_vertices{ a, b, c, d } {}

    Rectangle2D(const Vector2D &pos, real_t width, real_t height) : m_vertices{ pos, pos + Vector2D(width, 0), pos + Vector2D(width, height), pos + Vector2D(0, height) } {}

    ShapeType2D type() const override { return SHAPE2D_RECTANGLE; }
    static constexpr ShapeType2D shapeType = SHAPE2D_RECTANGLE;

    Vector2D& a() { return m_vertices[0]; }
    Vector2D& b() { return m_vertices[1]; }
    Vector2D& c() { return m_vertices[2]; }
    Vector2D& d() { return m_vertices[3]; }

    const Vector2D& a() const { return m_vertices[0]; }
    const Vector2D& b() const { return m_vertices[1]; }
    const Vector2D& c() const { return m_vertices[2]; }
    const Vector2D& d() const { return m_vertices[3]; }

    size_t vertexCount() const
    {
        return 4;
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

    real_t width() const
    {
        return (b() - a()).length();
    }

    real_t height() const
    {
        return (d() - a()).length();
    }

    real_t area() const override
    {
        return width() * height();
    }

    real_t perimeter() const override
    {
        return 2.0f * (width() + height());
    }

    Vector2D centroid() const override
    {
        return (a() + b() + c() + d()) / 4.0f;
    }

    Rectangle2D& translate(const Vector2D &translation) override
    {
        for (auto& v : m_vertices)
            v += translation;
        return *this;
    }

    Rectangle2D& rotate(real_t angle)
    {
        const Vector2D c = centroid();
        for (auto& v : m_vertices)
            v.rotateAround(angle, c);
        return *this;
    }

    Rectangle2D& rotate(real_t angle, const Vector2D& point)
    {
        for (auto& v : m_vertices)
            v.rotateAround(angle, point);
        return *this;
    }

    Rectangle2D& scale(real_t factor)
    {
        const Vector2D c = centroid();
        for (auto& v : m_vertices)
            v = c + (v - c) * factor;
        return *this;
    }

    Rectangle2D copy() const
    {
        return *this;
    }

    std::unique_ptr<IFiniteShape2D> clone() const
    {
        return std::make_unique<Rectangle2D>(*this);        
    }

    bool contains(const Vector2D &point) const override
    {
        const auto& m_a = this->a();
        const auto& m_b = this->b();
        const auto& m_c = this->c();
        const auto& m_d = this->d();

        const Vector2D ab = m_b - m_a;
        const Vector2D ap = point - m_a;
        const Vector2D bc = m_c - m_b;
        const Vector2D bp = point - m_b;
        const Vector2D cd = m_d - m_c;
        const Vector2D cp = point - m_c;
        const Vector2D da = m_a - m_d;
        const Vector2D dp = point - m_d;

        return ab.cross(ap) >= 0 && bc.cross(bp) >= 0 && cd.cross(cp) >= 0 && da.cross(dp) >= 0;
    }

    bool contains(const Rectangle2D &rectangle) const
    {
        return contains(rectangle.a()) && contains(rectangle.b()) && contains(rectangle.c()) && contains(rectangle.d());
    }

    BBox2D boundingBox() const override
    {
        return BBox2D(Vector2D::min(a(), b(), c(), d()), Vector2D::max(a(), b(), c(), d()));
    }

    Vector2D& operator[](size_t index)
    {
        return vertexAt(index);
    }

    const Vector2D& operator[](size_t index) const
    {
        return vertexAt(index);
    }
};

} // namespace Math

} // namespace Arns