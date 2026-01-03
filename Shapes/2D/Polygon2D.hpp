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
#include "Shapes/2D/Algorithms2D/Algorithms2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"

namespace Arns
{

namespace Math
{

class Polygon2D : public IFiniteShape2D
{
public:
    std::vector<Vector2D> m_vertices;

    Polygon2D() : m_vertices() {}

    Polygon2D(const std::vector<Vector2D> &vertices) : m_vertices(vertices) {}

    Polygon2D(std::vector<Vector2D>&& vertices) noexcept : m_vertices(std::move(vertices)) {}

    ShapeType2D type() const
    {
        return SHAPE2D_POLYGON;
    }

    static constexpr ShapeType2D shapeType = SHAPE2D_POLYGON;

    size_t vertexCount() const
    {
        return m_vertices.size();
    }

    Vector2D vertexAt(size_t index) const
    {
        return m_vertices[index];
    }

    Vector2D wrappedVertexAt(std::ptrdiff_t index) const
    {
        size_t n = m_vertices.size();
        if (n == 0)
            throw std::out_of_range("Cannot access element in empty vector");
        
        std::ptrdiff_t wrapped = index % static_cast<std::ptrdiff_t>(n);
        if (wrapped < 0)
            wrapped += n;
        return m_vertices[wrapped];
    }

    std::vector<Vector2D> getVertices() const
    {
        return m_vertices;
    }

    using iterator = std::vector<Vector2D>::iterator;
    using const_iterator = std::vector<Vector2D>::const_iterator;

    iterator begin() { return m_vertices.begin(); }
    iterator end() { return m_vertices.end(); }

    const_iterator begin() const { return m_vertices.begin(); }
    const_iterator end() const { return m_vertices.end(); }

    bool isConvex() const
    {
        const size_t n = m_vertices.size();
        if (n < 3) 
            return false;

        int sign = 0;

        for (size_t i = 0; i < n; ++i)
        {
            const Vector2D& a = m_vertices[i];
            const Vector2D& b = m_vertices[(i + 1) % n];
            const Vector2D& c = m_vertices[(i + 2) % n];

            // Compute edge vectors
            const Vector2D ab = b - a;
            const Vector2D bc = c - b;

            real_t cross = ab.cross(bc);

            if (approximatelyZero(cross))
            {
                // Collinear points
                continue;
            }

            int currentSign = (cross > 0) ? 1 : -1;

            if (sign == 0)
            {
                // First turn
                sign = currentSign;
            }
            else if (currentSign != sign)
            {
                // Turning direction changes -> concave
                return false;
            }
        }
        return true;
    }

    real_t area() const
    {
        double area = 0;
        for (int i = 0; i < m_vertices.size(); i++)
        {
            const Vector2D& a = m_vertices[i];
            const Vector2D& b = m_vertices[(i + 1 == m_vertices.size()) ? 0 : i + 1];
            area += static_cast<double>(a.x) * static_cast<double>(b.y)
              - static_cast<double>(b.x) * static_cast<double>(a.y);
        }

        return area / 2.0f;
    }

    real_t perimeter() const
    {
        real_t perimeter = 0;

        for (int i = 0; i < m_vertices.size(); i++)
        {
            perimeter += (m_vertices[i] - m_vertices[(i + 1) % m_vertices.size()]).length();
        }

        return perimeter;
    }

    Vector2D centroid() const
    {
        Vector2D centroid;

        for (int i = 0; i < m_vertices.size(); i++)
        {
            centroid += m_vertices[i];
        }

        return centroid / m_vertices.size();
    }

    Polygon2D& translate(const Vector2D &translation)
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i] += translation;
        return *this;
    }

    Polygon2D& rotate(real_t angle)
    {
        const Vector2D centroid = this->centroid();
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i].rotateAround(angle, centroid);
        return *this;
    }

    Polygon2D& rotate(real_t angle, const Vector2D& point)
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i].rotateAround(angle, point);
        return *this;
    }

    Polygon2D& scale(real_t factor)
    {
        const Vector2D centroid = this->centroid();
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i] = centroid + (m_vertices[i] - centroid) * factor;
        return *this;
    }

    Polygon2D copy() const
    {
        return *this;
    }

    std::unique_ptr<IFiniteShape2D> clone() const
    {
        return std::make_unique<Polygon2D>(*this);        
    }

    /*This only works for "simple" (non self intersecting) polygons*/
    bool contains(const Vector2D &point) const
    {
        bool contains = false;
        for (int i = 0, j = m_vertices.size() - 1; i < m_vertices.size(); j = i++)
        {
            if (isPointOnSegment(point, m_vertices[i], m_vertices[j]))
                return true;

            if (((m_vertices[i].y > point.y) != (m_vertices[j].y > point.y)) &&
                (point.x < (m_vertices[j].x - m_vertices[i].x) * (point.y - m_vertices[i].y) / (m_vertices[j].y - m_vertices[i].y) + m_vertices[i].x))
            {
                contains = !contains;
            }
        }

        return contains;
    }

    bool strictlyContains(const Vector2D &point) const
    {
        bool contains = false;
        for (int i = 0, j = m_vertices.size() - 1; i < m_vertices.size(); j = i++)
        {
            if (((m_vertices[i].y > point.y) != (m_vertices[j].y > point.y)) &&
                (point.x < (m_vertices[j].x - m_vertices[i].x) * (point.y - m_vertices[i].y) / (m_vertices[j].y - m_vertices[i].y) + m_vertices[i].x))
            {
                contains = !contains;
            }
        }

        return contains;
    }

    BBox2D boundingBox() const
    {
        if (m_vertices.empty())
            return BBox2D();

        Vector2D min = m_vertices[0];
        Vector2D max = m_vertices[0];

        for (int i = 1; i < m_vertices.size(); i++)
        {
            min.x = fmin(min.x, m_vertices[i].x);
            min.y = fmin(min.y, m_vertices[i].y);
            max.x = fmax(max.x, m_vertices[i].x);
            max.y = fmax(max.y, m_vertices[i].y);
        }

        return BBox2D(min, max);
    }

    Vector2D& operator[](size_t index)
    {
        return m_vertices[index];
    }

    const Vector2D& operator[](size_t index) const
    {
        return m_vertices[index];
    }
};


} // namespace Math

} // namespace Arns