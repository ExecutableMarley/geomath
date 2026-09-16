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

struct TriangleMeshData
{
    std::vector<Vector2D> m_vertices;
    std::vector<TriangleIndices> m_triangles;

    TriangleMeshData(const std::span<const Vector2D>& vertices, const std::span<const TriangleIndices>& triangles) :
        m_vertices(vertices.begin(), vertices.end()), m_triangles(triangles.begin(), triangles.end()) {}

    TriangleMeshData(std::span<const Vector2D> vertices,
        std::initializer_list<TriangleIndices> triangles)
        : m_vertices(vertices.begin(), vertices.end()),
          m_triangles(triangles)
    {}

    bool isValid() const 
    {
        const size_t vertexCount = m_vertices.size();
        for (const auto& t : m_triangles) 
        {
            if (static_cast<size_t>(t.v0) >= vertexCount ||
                static_cast<size_t>(t.v1) >= vertexCount ||
                static_cast<size_t>(t.v2) >= vertexCount) 
            {
                return false;
            }

            if (t.v0 == t.v1 || t.v1 == t.v2 || t.v2 == t.v0)
                return false;

            const Vector2D& v0 = m_vertices[t.v0];
            const Vector2D& v1 = m_vertices[t.v1];
            const Vector2D& v2 = m_vertices[t.v2];

            if (!isCCW(v0, v1, v2)) return false;
        }
        return true;
    }

    bool validate()
    {
        const size_t vertexCount = m_vertices.size();
        for (auto& t : m_triangles) 
        {
            if (static_cast<size_t>(t.v0) >= vertexCount ||
                static_cast<size_t>(t.v1) >= vertexCount ||
                static_cast<size_t>(t.v2) >= vertexCount) 
            {
                return false;
            }

            if (t.v0 == t.v1 || t.v1 == t.v2 || t.v2 == t.v0)
                return false;

            const Vector2D& v0 = m_vertices[t.v0];
            const Vector2D& v1 = m_vertices[t.v1];
            const Vector2D& v2 = m_vertices[t.v2];

            real_t signedArea = orient2D(v0, v1, v2);
            if (approximatelyZero(signedArea))
                return false;

            if (approximatelyLessAbs(signedArea, 0))
            {
                std::swap(t.v1, t.v2);
            }
        }
        return true;
    }
};

class TriangleMesh2D
{
private:
    TriangleMeshData m_data;

    //BBox2D
    //cached area
    //Vector2D centroid
public:

    // --- Constructors ---

    TriangleMesh2D() = default;

    explicit TriangleMesh2D(TriangleMeshData data) 
        : m_data(std::move(data)) {}

    TriangleMesh2D(const std::span<const Vector2D>& vertices, const std::span<const TriangleIndices>& triangles) :
        m_data(vertices, triangles) {}

    TriangleMesh2D(
        std::span<const Vector2D> vertices,
        std::initializer_list<TriangleIndices> triangles) :
            m_data(vertices, triangles) {}

    TriangleMeshData releaseData() &&
    {
        return std::move(m_data);
    }

    // --- Data Accessors & Mutators ---

    bool isValid() const { return m_data.isValid(); }

    const std::vector<Vector2D>& vertices()         const { return m_data.m_vertices; }
    const std::vector<TriangleIndices>& triangles() const { return m_data.m_triangles; }

    Triangle2D getTriangle(size_t index) const
    {
        if (index >= m_data.m_triangles.size())
            throw std::out_of_range("Triangle index out of range");

        const TriangleIndices& triIndices = m_data.m_triangles[index];
        return Triangle2D(m_data.m_vertices[triIndices.v0], m_data.m_vertices[triIndices.v1], m_data.m_vertices[triIndices.v2]);
    }

    // --- Geometric Properties ---

    // area

    // --- Transform / modification ---

    // translate, scale, rotate, transform

    // --- Spatial Queries ---

    // --- Operators ---

    const TriangleIndices& operator[](size_t index) const
    {
        return m_data.m_triangles[index];
    }
};

class TriangleMeshBuilder2D 
{
private:
    TriangleMeshData m_data;
public:

    TriangleMeshBuilder2D() = default;

    explicit TriangleMeshBuilder2D(TriangleMeshData data) 
        : m_data(std::move(data)) {}

    // --- Data Accessors & Mutators ---

    TriangleMeshBuilder2D& addVertex(Vector2D v) {
        m_data.m_vertices.push_back(std::move(v));
        return *this;
    }

    TriangleMeshBuilder2D& addVertex(real_t x, real_t y) {
        m_data.m_vertices.push_back({x, y});
        return *this;
    }

    TriangleMeshBuilder2D& addTriangle(TriangleIndices t) {
        m_data.m_triangles.push_back(std::move(t));
        return *this;
    }

    TriangleMeshBuilder2D& addTriangle(uint32_t i0, uint32_t i1, uint32_t i2) {
        m_data.m_triangles.push_back({i0, i1, i2});
        return *this;
    }

    TriangleMeshBuilder2D& reserve(size_t numVertices, size_t numTriangles) {
        m_data.m_vertices.reserve(numVertices);
        m_data.m_triangles.reserve(numTriangles);
        return *this;
    }

    // Move out (consuming/invalidating m_data)
    TriangleMeshData build() && 
    {
        return std::move(m_data);
    }

    // ivalue reference
    /*TriangleMeshData build() const & 
    {
        return m_data;
    }*/

    // ivalue reference
    TriangleMeshData copyData() const &
    {
        return m_data;
    }

    const TriangleMeshData& data() const & {
        return m_data;
    }
};



} // namespace Math

} // namespace Arns