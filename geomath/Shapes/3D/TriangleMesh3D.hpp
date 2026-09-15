/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>
#include <memory>

#include "Geometry/Vector3D.hpp"
#include "BBox3D.hpp"
#include "IShape3D.hpp"
#include "Triangle3D.hpp"
#include "Line3D.hpp"
#include "Shapes/Common/TriangleIndices.hpp"

#include "Intersections3D.hpp"

namespace Arns
{

namespace geomath
{

// Todo: Should have the same structure as TriangleMesh2D

class TriangleMesh3D
{
public:
    std::vector<Vector3D> m_vertices;
    std::vector<TriangleIndices> m_triangles;
    BBox3D m_boundingBox;

    TriangleMesh3D() : m_vertices(), m_triangles() {}

    TriangleMesh3D(const std::vector<Vector3D> &vertices, const std::vector<TriangleIndices> &indices) : m_vertices(vertices), m_triangles(indices) 
    {
        if (m_vertices.empty() || m_triangles.empty())
            return;
        // Calculate the bounds of the shapes
        m_boundingBox = BBox3D(vertices[0], vertices[0]);
        for (auto& vector : m_vertices)
        {
            m_boundingBox.encapsulate(vector);
        }
    }

    bool validate() const
    {
        // Check indices count is multiple of 3
        if (m_triangles.size() % 3 != 0)
            return false;

        // Check indices are within bounds
        for (size_t i = 0; i < m_triangles.size(); i += 3)
        {
            if (m_triangles[i].v0 >= m_vertices.size() || m_triangles[i].v1 >= m_vertices.size() || m_triangles[i].v2 >= m_vertices.size())
                return false;
        }
        return true;
    }

    unsigned int findClosestVertex(const Vector3D& point)
    {
        int iClosest = -1;
        float fClosest = std::numeric_limits<float>().max();
        for (int i = 0; i < m_vertices.size(); i++)
        {
            float dist = (point - m_vertices[i]).length();
            if (dist < fClosest)
            {
                fClosest = dist;
            }
        }
        return iClosest;
    }

    bool addTriangle(const Triangle3D& triangle)
    {
        return addTriangle(triangle.m_a, triangle.m_b, triangle.m_c);
    }

    bool addTriangle(const Vector3D& a, const Vector3D& b, const Vector3D& c)
    {
        //Todo: Merge vertex with existing one if close
        uint32_t offset = m_vertices.size();
        this->m_vertices.push_back(a);
        this->m_vertices.push_back(b);
        this->m_vertices.push_back(c);
        this->m_triangles.push_back({offset + 0, offset + 1, offset + 2});
    }

    bool addTriangles(const std::vector<Vector3D>& vertices, const std::vector<TriangleIndices>& indices)
    {
        // Check indices count is multiple of 3
        if (indices.size() % 3 != 0)
            return false;

        // Check indices are within bounds
        for (size_t i = 0; i < indices.size(); i += 3)
        {
            if (indices[i].v0 >= vertices.size() || indices[i].v1 >= vertices.size() || indices[i].v2 >= vertices.size())
                return false;
        }

        size_t offset = m_vertices.size();

        std::vector<TriangleIndices> newIndices;
        for (auto t : indices)
        {
            newIndices.push_back({t.v0 + offset, t.v1 + offset, t.v2 + offset});
        }
            //newIndices.push_back(indice + offset);

        m_vertices.insert(m_vertices.end(), vertices.begin(), vertices.end());
        m_triangles.insert(m_triangles.end(), newIndices.begin(), newIndices.end());

        return true;
    }

    bool addMesh(const TriangleMesh3D& mesh)
    {
        return addTriangles(mesh.m_vertices, mesh.m_triangles);
    }

    bool shorten(const BBox3D& bounds)
    {
        if (bounds.contains(this->m_boundingBox))
            return true;

        if (!this->m_boundingBox.contains(bounds))
            return false;

        std::vector<TriangleIndices> newIndices;
        for (int i = 0; i + 2 < m_triangles.size(); i = i + 3)
        {
            if (bounds.contains(m_vertices[m_triangles[i].v0])
                && bounds.contains(m_vertices[m_triangles[i].v1])
                && bounds.contains(m_vertices[m_triangles[i].v2]))
            {
                newIndices.push_back(m_triangles[i]);
                newIndices.push_back(m_triangles[i+1]);
                newIndices.push_back(m_triangles[i+2]);
            }
        }

        this->m_triangles = std::move(newIndices);
        this->m_boundingBox = bounds;
        return true;
    }

    /*
    void removeUnusedVertices()
    {
        std::vector<bool> used(m_vertices.size(), false);
        for (const auto& index : m_triangles)
        {
            used[index.v0] = true;
            used[index.v1] = true;
            used[index.v2] = true;
        }

        std::vector<Vector3D> newVertices;
        std::vector<TriangleIndices> remap(m_vertices.size(), 0);
        unsigned int newIndex = 0;

        for (size_t i = 0; i < m_vertices.size(); ++i)
        {
            if (used[i])
            {
                remap[i] = newIndex++;
                newVertices.push_back(m_vertices[i]);
            }
        }

        for (auto& index : m_triangles)
        {
            index = remap[index];
        }

        m_vertices = std::move(newVertices);
    }
    */

    Vector3D centroid() const
    {
        Vector3D centroid;
        for (size_t i = 0; i < m_triangles.size(); i += 3)
        {
            const Vector3D &a = m_vertices[m_triangles[i].v0];
            const Vector3D &b = m_vertices[m_triangles[i].v1];
            const Vector3D &c = m_vertices[m_triangles[i].v2];
            centroid += (a + b + c) / 3.0f;
        }
        return centroid / (m_triangles.size());
    }

    TriangleMesh3D& translate(const Vector3D &translation)
    {
        for (size_t i = 0; i < m_vertices.size(); ++i)
        {
            m_vertices[i] += translation;
        }
        m_boundingBox.translate(translation);
        return *this;
    }

    bool intersects(const Line3D &line) const
    {
        for (size_t i = 0; i < m_triangles.size(); i += 3)
        {
            const Vector3D &a = m_vertices[m_triangles[i].v0];
            const Vector3D &b = m_vertices[m_triangles[i].v1];
            const Vector3D &c = m_vertices[m_triangles[i].v2];

            if (geomath::intersects(line, Triangle3D(a, b, c)))
                return true;
        }
        return false;
    }

    // --- Data Accessors & Mutators ---

    // --- Geometric Properties ---

    // --- Transform / modification ---

    // --- Spatial Queries ---

    // --- Operators ---
};

} // namespace Math

} // namespace Arns