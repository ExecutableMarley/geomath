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
#include "Shapes/Common/TriangleIndices.hpp"

namespace Arns
{

namespace geomath
{


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

    // --- Data Accessors & Mutators ---

    Triangle2D getTriangle(size_t index) const
    {
        if (index >= m_triangles.size())
            throw std::out_of_range("Triangle index out of range");

        const TriangleIndices& triIndices = m_triangles[index];
        return Triangle2D(m_vertices[triIndices.v0], m_vertices[triIndices.v1], m_vertices[triIndices.v2]);
    }

    // --- Geometric Properties ---

    // area

    // --- Transform / modification ---

    // translate, scale, rotate, transform

    // --- Spatial Queries ---

    // --- Operators ---

    TriangleIndices& operator[](size_t index)
    {
        return m_triangles[index];
    }
};

class TriangleMeshBuilder2D 
{
};



} // namespace Math

} // namespace Arns