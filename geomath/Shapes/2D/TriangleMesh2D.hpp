#pragma once

#include <math.h>
#include <stdexcept>
#include <vector>
#include <span>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "Interfaces/IFiniteShape2D.hpp"
#include "Triangle2D.hpp"

namespace Arns
{

namespace geomath
{

struct TriangleIndices
{
    TriangleIndices(int v0, int v1, int v2) : v0(v0), v1(v1), v2(v2) {}

    int v0, v1, v2;

    unsigned int operator[](size_t index) const
    {
        switch (index)
        {
        case 0: return v0;
        case 1: return v1;
        case 2: return v2;
        default:
            return 0; //Throw exception
        }
    }
};

// Todo: Should have the same structure as TriangleMesh3D

class TriangleMesh2D
{
public:
    std::vector<Vector2D> m_vertices;
    std::vector<TriangleIndices> m_triangles;

    // --- Constructors ---

    TriangleMesh2D() : m_vertices(), m_triangles() {};

    TriangleMesh2D(const std::span<const Vector2D>& vertices, const std::span<const TriangleIndices>& triangles) :
        m_vertices(vertices.begin(), vertices.end()), m_triangles(triangles.begin(), triangles.end()) {}

    TriangleMesh2D(
        std::span<const Vector2D> vertices,
        std::initializer_list<TriangleIndices> triangles)
        : m_vertices(vertices.begin(), vertices.end()),
          m_triangles(triangles)
    {
    }

    bool validate() const
    {
        // Check indices are within bounds
        for (size_t i = 0; i < m_triangles.size(); i += 1)
        {
            if (m_triangles[i].v0 >= m_vertices.size() || m_triangles[i].v1 >= m_vertices.size() || m_triangles[i].v2 >= m_vertices.size())
                return false;
        }
        return true;
    }

    void addTriangle(const Triangle2D& triangle);

    void addTriangle(const Vector2D& a, const Vector2D& b, const Vector2D& c);

    Triangle2D getTriangle(size_t index) const
    {
        if (index >= m_triangles.size())
            throw std::out_of_range("Triangle index out of range");

        const TriangleIndices& triIndices = m_triangles[index];
        return Triangle2D(m_vertices[triIndices.v0], m_vertices[triIndices.v1], m_vertices[triIndices.v2]);
    }
};



} // namespace Math

} // namespace Arns