/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Algebra/Transform2D.hpp"
#include "Shapes/Common/TriangleIndices.hpp"
#include "Shapes/2D/BBox2D.hpp"
#include "Shapes/2D/Triangle2D.hpp"

#include <math.h>
#include <stdexcept>
#include <vector>
#include <span>
#include <initializer_list>

namespace Arns::geomath
{


struct TriangleMeshData2D
{
    std::vector<Vector2D> m_vertices;
    std::vector<TriangleIndices> m_triangles;

    TriangleMeshData2D() = default;

    TriangleMeshData2D(std::span<const Vector2D> vertices, std::span<const TriangleIndices> triangles) :
        m_vertices(vertices.begin(), vertices.end()), m_triangles(triangles.begin(), triangles.end()) {}

    TriangleMeshData2D(std::span<const Vector2D> vertices, std::initializer_list<TriangleIndices> triangles)
        : m_vertices(vertices.begin(), vertices.end()), m_triangles(triangles) {}

    TriangleMeshData2D(std::initializer_list<Vector2D> vertices, std::initializer_list<TriangleIndices> triangles)
        : m_vertices(vertices), m_triangles(triangles) {}

    void reserve(size_t numVertices, size_t numTriangles) 
    {
        m_vertices.reserve(numVertices);
        m_triangles.reserve(numTriangles);
    }

    bool isValid() const 
    {
        const size_t vertexCount = m_vertices.size();
        for (const TriangleIndices& t : m_triangles) 
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
    // Base Members

    TriangleMeshData2D m_data;
    Transform2D m_transform;

    // Derived local data from Mesh

    BBox2D   m_localBounds;
    real_t   m_area;
    Vector2D m_centroid;
    
    // Cached world data

    mutable BBox2D m_worldBounds;
    mutable bool m_worldBoundsDirty = true;


    void _invalidate()
    {
        m_worldBoundsDirty = true;
    }

    void _derivedData()
    {
        if (this->m_data.m_triangles.empty())
        {
            m_localBounds = BBox2D();
            m_area = 0;
            m_centroid = Vector2D();
            return;
        }

        m_localBounds = BBox2D::fromPoint(m_data.m_vertices[m_data.m_triangles[0].v0]);
        real_t totalArea = real_t{0};
        Vector2D weightedCentroid = Vector2D::zero();

        for (const auto& t : m_data.m_triangles)
        {
            const Vector2D& v0 = m_data.m_vertices[t.v0];
            const Vector2D& v1 = m_data.m_vertices[t.v1];
            const Vector2D& v2 = m_data.m_vertices[t.v2];

            m_localBounds.encapsulate(v0);
            m_localBounds.encapsulate(v1);
            m_localBounds.encapsulate(v2);

            const real_t triangleArea = cross(v1 - v0, v2 - v0) * real_t{0.5};
            const Vector2D triangleCentroid = (v0 + v1 + v2) / real_t{3};

            weightedCentroid += triangleCentroid * triangleArea;
            totalArea        += triangleArea;
        }

        this->m_area = totalArea;

        if (!approximatelyZero(totalArea))
            this->m_centroid = weightedCentroid / totalArea;
    }

public:

    // --- Constructors ---

    TriangleMesh2D() = default;

    explicit TriangleMesh2D(TriangleMeshData2D data) 
        : m_data(std::move(data)) {
            this->_derivedData();
        }

    TriangleMesh2D(const std::span<const Vector2D>& vertices, const std::span<const TriangleIndices>& triangles) :
        m_data(vertices, triangles) {
            this->_derivedData();
        }

    TriangleMesh2D(
        std::span<const Vector2D> vertices,
        std::initializer_list<TriangleIndices> triangles) :
            m_data(vertices, triangles) {
                this->_derivedData();
            }

    TriangleMesh2D(
        std::initializer_list<Vector2D> vertices,
        std::initializer_list<TriangleIndices> triangles) :
            m_data(vertices, triangles) {
                this->_derivedData();
            }

    TriangleMeshData2D releaseData() &&
    {
        return std::move(m_data);
    }

    // --- Data Accessors & Mutators ---

    bool isValid() const { return m_data.isValid(); }

    const std::vector<Vector2D>& vertices()         const { return m_data.m_vertices; }
    const std::vector<TriangleIndices>& triangles() const { return m_data.m_triangles; }
    const TriangleMeshData2D& data()                const { return m_data;}

    Triangle2D getTriangle(size_t index) const
    {
        if (index >= m_data.m_triangles.size())
            throw std::out_of_range("Triangle index out of range");

        const TriangleIndices& triIndices = m_data.m_triangles[index];
        return Triangle2D(m_data.m_vertices[triIndices.v0], m_data.m_vertices[triIndices.v1], m_data.m_vertices[triIndices.v2]);
    }

    Triangle2D getTriangleWorld(size_t index) const
    {
        if (index >= m_data.m_triangles.size())
            throw std::out_of_range("Triangle index out of range");

        const TriangleIndices& t = m_data.m_triangles[index];
        const Vector2D& a = m_transform.transformPoint(m_data.m_vertices[t.v0]);
        const Vector2D& b = m_transform.transformPoint(m_data.m_vertices[t.v1]);
        const Vector2D& c = m_transform.transformPoint(m_data.m_vertices[t.v2]);
        return Triangle2D(a,b,c);
    }

    // --- Geometric Properties ---

    real_t localArea() const { return m_area; }
    real_t worldArea() const { return m_area * m_transform.getScale() * m_transform.getScale(); }

    Vector2D localCentroid() const { return m_centroid; }
    Vector2D worldCentroid() const { return m_transform.transformPoint(m_centroid); }

    BBox2D localBounds () const { return m_localBounds; }
    BBox2D worldBounds () const
    {
        if (m_worldBoundsDirty)
        {
            m_worldBounds = BBox2D::fromPoints(
                m_transform.transformPoint(m_localBounds.topLeft()),
                m_transform.transformPoint(m_localBounds.topRight()),
                m_transform.transformPoint(m_localBounds.bottomRight()),
                m_transform.transformPoint(m_localBounds.bottomLeft()));

            m_worldBoundsDirty = false;
        }
        return m_worldBounds;
    }

    // --- Transform / modification ---

    TriangleMesh2D& translate(const Vector2D& translation)
    {
        this->m_transform.translate(translation);
        _invalidate();
        return *this;
    }

    TriangleMesh2D& rotate(Radians rotate)
    {
        this->m_transform.rotate(rotate);
        _invalidate();
        return *this;
    }

    TriangleMesh2D& scale(real_t factor, const Vector2D& point)
    {
        const Vector2D position = m_transform.getPosition();
        m_transform.setScale(m_transform.getScale() * factor);
        m_transform.setPosition(point + (position - point) * factor);
        _invalidate();
        return *this;
    }

    TriangleMesh2D& scale(real_t factor)
    {
        return scale(factor, this->worldCentroid());
    }

    TriangleMesh2D& transform(Transform2D transform)
    {
        this->m_transform *= transform;
        _invalidate();
        return *this;
    }

    void bakeTransform()
    {
        for (auto& v : m_data.m_vertices)
        {
            v = m_transform.transformPoint(v);
        }
        m_transform.reset();
        _derivedData();
        _invalidate();
    }


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
    TriangleMeshData2D m_data;
public:

    TriangleMeshBuilder2D() = default;
    explicit TriangleMeshBuilder2D(TriangleMeshData2D data) : m_data(std::move(data)) {}

    // --- Data Accessors & Mutators ---

    TriangleMeshBuilder2D& addVertex(Vector2D v) {
        m_data.m_vertices.push_back(std::move(v));
        return *this;
    }

    TriangleMeshBuilder2D& addVertex(real_t x, real_t y) {
        m_data.m_vertices.emplace_back(x, y);
        return *this;
    }

    TriangleMeshBuilder2D& addTriangle(TriangleIndices t) {
        m_data.m_triangles.push_back(std::move(t));
        return *this;
    }

    TriangleMeshBuilder2D& addTriangle(uint32_t i0, uint32_t i1, uint32_t i2) {
        m_data.m_triangles.emplace_back(i0, i1, i2);
        return *this;
    }

    TriangleMeshBuilder2D& addVertices(const std::vector<Vector2D>& vertices) {
        m_data.m_vertices.insert(m_data.m_vertices.end(), vertices.begin(), vertices.end());
        return *this;
    }

    TriangleMeshBuilder2D& addTriangles(const std::vector<TriangleIndices>& triangles) {
        m_data.m_triangles.insert(m_data.m_triangles.end(), triangles.begin(), triangles.end());
        return *this;
    }

    TriangleMeshBuilder2D& merge(const TriangleMeshData2D& other) {
        uint32_t vertexOffset = static_cast<uint32_t>(m_data.m_vertices.size());
        
        // Reserve space to avoid reallocations
        m_data.m_vertices.reserve(m_data.m_vertices.size() + other.m_vertices.size());
        m_data.m_triangles.reserve(m_data.m_triangles.size() + other.m_triangles.size());

        // Append vertices
        m_data.m_vertices.insert(m_data.m_vertices.end(), other.m_vertices.begin(), other.m_vertices.end());

        // Append triangles with shifted indices
        for (const auto& t : other.m_triangles) {
            m_data.m_triangles.push_back({
                t.v0 + vertexOffset,
                t.v1 + vertexOffset,
                t.v2 + vertexOffset
            });
        }
        return *this;
    }

    TriangleMeshBuilder2D& reserve(size_t numVertices, size_t numTriangles) {
        m_data.reserve(numVertices, numTriangles);
        return *this;
    }

    // Removes triangles where the predicate returns true
    template <typename Predicate>
    TriangleMeshBuilder2D& filterTriangles(Predicate&& pred) {
        auto& tris = m_data.m_triangles;
        tris.erase(
            std::remove_if(tris.begin(), tris.end(), [&](const TriangleIndices& t) {
                return pred(m_data.m_vertices[t.v0], 
                            m_data.m_vertices[t.v1], 
                            m_data.m_vertices[t.v2]);
            }),
            tris.end()
        );
        return *this;
    }

    template <typename Func>
    TriangleMeshBuilder2D& transformVertices(Func&& func) {
        for (auto& v : m_data.m_vertices) {
            func(v);
        }
        return *this;
    }

    // Removes orphan vertices that are no longer referenced by any triangle
    TriangleMeshBuilder2D& removeUnusedVertices() {
        std::vector<bool> usedVertex(m_data.m_vertices.size(), false);
        for (const auto& t : m_data.m_triangles) {
            usedVertex[t.v0] = true;
            usedVertex[t.v1] = true;
            usedVertex[t.v2] = true;
        }

        std::vector<uint32_t> oldToNewIndex(m_data.m_vertices.size());
        std::vector<Vector2D> compactedVertices;
        compactedVertices.reserve(m_data.m_vertices.size());

        for (size_t i = 0; i < m_data.m_vertices.size(); ++i) {
            if (usedVertex[i]) {
                oldToNewIndex[i] = static_cast<uint32_t>(compactedVertices.size());
                compactedVertices.push_back(m_data.m_vertices[i]);
            }
        }

        // Remap triangle indices
        for (auto& t : m_data.m_triangles) {
            t.v0 = oldToNewIndex[t.v0];
            t.v1 = oldToNewIndex[t.v1];
            t.v2 = oldToNewIndex[t.v2];
        }

        m_data.m_vertices = std::move(compactedVertices);
        return *this;
    }

    //Optimize Method?

    // Move out (consuming/invalidating m_data)
    TriangleMeshData2D build() && 
    {
        return std::move(m_data);
    }

    // ivalue reference
    TriangleMeshData2D copyData() const &
    {
        return m_data;
    }

    const TriangleMeshData2D& data() const & {
        return m_data;
    }
};


} // namespace Arns::geomath