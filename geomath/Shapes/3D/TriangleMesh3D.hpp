/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector3D.hpp"
#include "Algebra/Transform3D.hpp"
#include "Shapes/Common/TriangleIndices.hpp"
#include "Shapes/3D/BBox3D.hpp"
#include "Shapes/3D/Triangle3D.hpp"

#include <cstdint>
#include <initializer_list>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Arns::geomath
{

struct TriangleMeshData3D
{
    std::vector<Vector3D> m_vertices;
    std::vector<TriangleIndices> m_triangles;

    TriangleMeshData3D() = default;

    TriangleMeshData3D(std::span<const Vector3D> vertices, std::span<const TriangleIndices> triangles)
        : m_vertices(vertices.begin(), vertices.end()), m_triangles(triangles.begin(), triangles.end()) {}

    TriangleMeshData3D(std::span<const Vector3D> vertices, std::initializer_list<TriangleIndices> triangles)
        : m_vertices(vertices.begin(), vertices.end()), m_triangles(triangles) {}

    TriangleMeshData3D(std::initializer_list<Vector3D> vertices, std::initializer_list<TriangleIndices> triangles)
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

            const Vector3D& v0 = m_vertices[t.v0];
            const Vector3D& v1 = m_vertices[t.v1];
            const Vector3D& v2 = m_vertices[t.v2];

            const Vector3D crossProduct = (v1 - v0).cross(v2 - v0);
            if (approximatelyZero(crossProduct.lengthSquared()))
                return false;

            //Todo: Winding order
        }
        return true;
    }

    bool validate() const { return isValid(); }
};

class TriangleMesh3D
{
private:
    // Base Members

    TriangleMeshData3D m_data;
    Transform3D m_transform;

    // Derived local data from Mesh

    BBox3D   m_localBounds;
    real_t   m_surfaceArea = 0;
    Vector3D m_centroid;

    // Cached world data

    mutable BBox3D m_worldBounds;
    mutable bool m_worldBoundsDirty = true;

    void _invalidate()
    {
        m_worldBoundsDirty = true;
    }

    void _derivedData()
    {
        if (this->m_data.m_triangles.empty())
        {
            m_localBounds = BBox3D();
            m_surfaceArea = 0;
            m_centroid = Vector3D();
            return;
        }

        m_localBounds = BBox3D::fromPoint(m_data.m_vertices[m_data.m_triangles[0].v0]);
        real_t totalArea = real_t{0};
        Vector3D weightedCentroid = Vector3D::zero();

        for (const TriangleIndices& indices : m_data.m_triangles)
        {
            const Vector3D& v0 = m_data.m_vertices[indices.v0];
            const Vector3D& v1 = m_data.m_vertices[indices.v1];
            const Vector3D& v2 = m_data.m_vertices[indices.v2];
            
            m_localBounds.encapsulate(v0);
            m_localBounds.encapsulate(v1);
            m_localBounds.encapsulate(v2);
            
            const real_t triangleArea = (v1 - v0).cross(v2 - v0).length() * real_t{0.5};
            const Vector3D triangleCentroid = (v0 + v1 + v2) / real_t{3};

            weightedCentroid += triangleCentroid * triangleArea;
            totalArea        += triangleArea;
        }

        this->m_surfaceArea = totalArea;

        if (!approximatelyZero(m_surfaceArea))
            m_centroid = weightedCentroid / m_surfaceArea;
    }

public:

    // --- Constructors ---

    TriangleMesh3D() = default;

    explicit TriangleMesh3D(TriangleMeshData3D data)
        : m_data(std::move(data)) {
            this->_derivedData();
        }

    TriangleMesh3D(std::span<const Vector3D> vertices, std::span<const TriangleIndices> triangles)
        : m_data(vertices, triangles) {
            this->_derivedData();
        }

    TriangleMesh3D(std::span<const Vector3D> vertices, std::initializer_list<TriangleIndices> triangles)
        : m_data(vertices, triangles) {
            this->_derivedData();
        }

    TriangleMesh3D(std::initializer_list<Vector3D> vertices, std::initializer_list<TriangleIndices> triangles)
        : m_data(vertices, triangles) {
            this->_derivedData();
        }

    TriangleMeshData3D releaseData() && { return std::move(m_data); }

    // --- Data Accessors & Mutators ---

    bool isValid() const { return m_data.isValid(); }

    const std::vector<Vector3D>& vertices()         const { return m_data.m_vertices; }
    const std::vector<TriangleIndices>& triangles() const { return m_data.m_triangles; }
    const TriangleMeshData3D& data()                const { return m_data; }

    Triangle3D getTriangle(size_t index) const
    {
        if (index >= m_data.m_triangles.size())
            throw std::out_of_range("Triangle index out of range");
        
        const TriangleIndices& indices = m_data.m_triangles[index];
        return Triangle3D(m_data.m_vertices[indices.v0], m_data.m_vertices[indices.v1], m_data.m_vertices[indices.v2]);
    }

    Triangle3D getTriangleWorld(size_t index) const
    {
        if (index >= m_data.m_triangles.size())
            throw std::out_of_range("Triangle index out of range");
        
        const TriangleIndices& indices = m_data.m_triangles[index];
        return Triangle3D(
            m_transform.transformPoint(m_data.m_vertices[indices.v0]),
            m_transform.transformPoint(m_data.m_vertices[indices.v1]),
            m_transform.transformPoint(m_data.m_vertices[indices.v2]));
    }

    real_t localSurfaceArea() const { return m_surfaceArea; }
    real_t worldSurfaceArea() const { return m_surfaceArea * m_transform.getScale() * m_transform.getScale(); }

    Vector3D localCentroid() const { return m_centroid; }
    Vector3D worldCentroid() const { return m_transform.transformPoint(m_centroid); }

    BBox3D localBounds() const { return m_localBounds; }
    BBox3D worldBounds() const
    {
        if (m_worldBoundsDirty)
        {
            if (m_data.m_triangles.empty() || !m_data.isValid())
            {
                m_worldBounds = BBox3D();
            }
            else
            {
                const Vector3D& minimum = m_localBounds.m_min;
                const Vector3D& maximum = m_localBounds.m_max;
                bool firstCorner = true;
                for (unsigned int corner = 0; corner < 8; ++corner)
                {
                    const Vector3D point(
                        (corner & 1) ? maximum.x : minimum.x,
                        (corner & 2) ? maximum.y : minimum.y,
                        (corner & 4) ? maximum.z : minimum.z);
                    const Vector3D worldPoint = m_transform.transformPoint(point);
                    if (firstCorner)
                    {
                        m_worldBounds = BBox3D(worldPoint, worldPoint);
                        firstCorner = false;
                    }
                    else
                    {
                        m_worldBounds.encapsulate(worldPoint);
                    }
                }
            }
            m_worldBoundsDirty = false;
        }
        return m_worldBounds;
    }

    // --- Transform / modification ---

    TriangleMesh3D& translate(const Vector3D& translation)
    {
        m_transform.translate(translation);
        _invalidate();
        return *this;
    }

    TriangleMesh3D& rotate(const EulerAngles& rotation)
    {
        m_transform.rotate(rotation);
        _invalidate();
        return *this;
    }

    TriangleMesh3D& rotate(Radians pitch, Radians yaw, Radians roll)
    {
        return rotate(EulerAngles(pitch, yaw, roll));
    }

    TriangleMesh3D& scale(real_t factor, const Vector3D& point)
    {
        const Vector3D position = m_transform.getPosition();
        m_transform.setScale(m_transform.getScale() * factor);
        m_transform.setPosition(point + (position - point) * factor);
        _invalidate();
        return *this;
    }

    TriangleMesh3D& scale(real_t factor)
    {
        return scale(factor, this->worldCentroid());
    }

    TriangleMesh3D& transform(const Transform3D& transform)
    {
        m_transform *= transform;
        _invalidate();
        return *this;
    }

    void bakeTransform()
    {
        for (Vector3D& vertex : m_data.m_vertices)
            vertex = m_transform.transformPoint(vertex);
        m_transform.reset();
        _derivedData();
        _invalidate();
    }

    // --- Spatial Queries ---

    // --- Operators ---

    const TriangleIndices& operator[](size_t index) const { return m_data.m_triangles[index]; }
};

class TriangleMeshBuilder3D
{
private:
    TriangleMeshData3D m_data;

public:
    TriangleMeshBuilder3D() = default;
    explicit TriangleMeshBuilder3D(TriangleMeshData3D data) : m_data(std::move(data)) {}

    // --- Data Accessors & Mutators ---

    TriangleMeshBuilder3D& addVertex(Vector3D vertex) {
        m_data.m_vertices.push_back(std::move(vertex));
        return *this;
    }

    TriangleMeshBuilder3D& addVertex(real_t x, real_t y, real_t z) {
        m_data.m_vertices.emplace_back(x, y, z);
        return *this;
    }

    TriangleMeshBuilder3D& addTriangle(TriangleIndices triangle) {
        m_data.m_triangles.push_back(std::move(triangle));
        return *this;
    }

    TriangleMeshBuilder3D& addTriangle(uint32_t i0, uint32_t i1, uint32_t i2) {
        m_data.m_triangles.emplace_back(i0, i1, i2);
        return *this;
    }

    TriangleMeshBuilder3D& addVertices(const std::vector<Vector3D>& vertices) {
        m_data.m_vertices.insert(m_data.m_vertices.end(), vertices.begin(), vertices.end());
        return *this;
    }

    TriangleMeshBuilder3D& addTriangles(const std::vector<TriangleIndices>& triangles) {
        m_data.m_triangles.insert(m_data.m_triangles.end(), triangles.begin(), triangles.end());
        return *this;
    }

    TriangleMeshBuilder3D& merge(const TriangleMeshData3D& other) {
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

    TriangleMeshBuilder3D& reserve(size_t numVertices, size_t numTriangles) {
        m_data.reserve(numVertices, numTriangles);
        return *this;
    }

    // Removes triangles where the predicate returns true
    template <typename Predicate>
    TriangleMeshBuilder3D& filterTriangles(Predicate&& pred) {
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
    TriangleMeshBuilder3D& transformVertices(Func&& func) {
        for (auto& v : m_data.m_vertices) {
            func(v);
        }
        return *this;
    }

    // Removes orphan vertices that are no longer referenced by any triangle
    TriangleMeshBuilder3D& removeUnusedVertices() {
        std::vector<bool> usedVertex(m_data.m_vertices.size(), false);
        for (const auto& t : m_data.m_triangles) {
            usedVertex[t.v0] = true;
            usedVertex[t.v1] = true;
            usedVertex[t.v2] = true;
        }

        std::vector<uint32_t> oldToNewIndex(m_data.m_vertices.size());
        std::vector<Vector3D> compactedVertices;
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
    TriangleMeshData3D build() && { return std::move(m_data); }

    // ivalue reference
    TriangleMeshData3D copyData() const & { return m_data; }
    
    const TriangleMeshData3D& data() const & { return m_data; }
};

}