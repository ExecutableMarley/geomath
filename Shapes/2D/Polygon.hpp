/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>

#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"

namespace Utility
{

namespace Math
{

class Polygon
{
public:
    std::vector<Vector2D> m_vertices;

    Polygon() : m_vertices() {}

    Polygon(const std::vector<Vector2D> &vertices) : m_vertices(vertices) {}

    size_t vertexCount() const
    {
        return m_vertices.size();
    }

    std::vector<Vector2D> getVertices() const
    {
        return m_vertices;
    }

    Vector2D getVertexAt(size_t index) const
    {
        return m_vertices[index];
    }

    float area() const
    {
        float area = 0;

        for (int i = 0; i < m_vertices.size(); i++)
        {
            area += m_vertices[i].x * m_vertices[(i + 1) % m_vertices.size()].y - m_vertices[(i + 1) % m_vertices.size()].x * m_vertices[i].y;
        }

        return area / 2;
    }

    float perimeter() const
    {
        float perimeter = 0;

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

    Polygon& translate(const Vector2D &translation)
    {
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i] += translation;
        return *this;
    }

    Polygon& rotate(float angle)
    {
        const Vector2D centroid = this->centroid();
        for (int i = 0; i < m_vertices.size(); i++)
            m_vertices[i].rotateAround(angle, centroid);
        return *this;
    }

    BBox2D boundingBox() const
    {
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
};


} // namespace Math

} // namespace Utility