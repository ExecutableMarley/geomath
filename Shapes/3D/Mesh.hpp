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

#include "Intersections3D.hpp"

namespace Utility
{

namespace Math
{

class Mesh
{
public:
    std::vector<Vector3D> m_vertices;
    std::vector<unsigned int> m_indices;
    BBox3D m_boundingBox;

    Mesh() : m_vertices(), m_indices() {}

    Mesh(const std::vector<Vector3D> &vertices, const std::vector<unsigned int> &indices) : m_vertices(vertices), m_indices(indices) 
    {
        if (m_vertices.empty() || m_indices.empty())
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
        if (m_indices.size() % 3 != 0)
            return false;

        // Check indices are within bounds
        for (size_t i = 0; i < m_indices.size(); i += 3)
        {
            if (m_indices[i] >= m_vertices.size() || m_indices[i + 1] >= m_vertices.size() || m_indices[i + 2] >= m_vertices.size())
                return false;
        }
        return true;
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
        this->m_indices.push_back(offset + 0);
        this->m_indices.push_back(offset + 1);
        this->m_indices.push_back(offset + 2);
    }

    bool addTriangles(const std::vector<Vector3D>& vertices, const std::vector<unsigned int>& indices)
    {
        // Check indices count is multiple of 3
        if (indices.size() % 3 != 0)
            return false;

        // Check indices are within bounds
        for (size_t i = 0; i < indices.size(); i += 3)
        {
            if (indices[i] >= vertices.size() || indices[i + 1] >= vertices.size() || indices[i + 2] >= vertices.size())
                return false;
        }

        size_t offset = m_vertices.size();

        std::vector<unsigned int> newIndices;
        for (auto indice : indices)
            newIndices.push_back(indice + offset);

        m_vertices.insert(m_vertices.end(), vertices.begin(), vertices.end());
        m_indices.insert(m_indices.end(), newIndices.begin(), newIndices.end());

        return true;
    }

    bool addMesh(const Mesh& mesh)
    {
        return addTriangles(mesh.m_vertices, mesh.m_indices);
    }

    bool shorten(const BBox3D& bounds)
    {
        if (bounds.contains(this->m_boundingBox))
            return true;

        if (!this->m_boundingBox.contains(bounds))
            return false;

        std::vector<unsigned int> newIndices;
        for (int i = 0; i + 2 < m_indices.size(); i = i + 3)
        {
            if (bounds.contains(m_vertices[m_indices[i]])
                && bounds.contains(m_vertices[m_indices[i+1]])
                && bounds.contains(m_vertices[m_indices[i+2]]))
            {
                newIndices.push_back(m_indices[i]);
                newIndices.push_back(m_indices[i+1]);
                newIndices.push_back(m_indices[i+2]);
            }
        }

        this->m_indices = std::move(newIndices);
        this->m_boundingBox = bounds;
        return true;
    }

    void removeUnusedVertices()
    {
        std::vector<bool> used(m_vertices.size(), false);
        for (const auto& index : m_indices)
        {
            used[index] = true;
        }

        std::vector<Vector3D> newVertices;
        std::vector<unsigned int> remap(m_vertices.size(), 0);
        unsigned int newIndex = 0;

        for (size_t i = 0; i < m_vertices.size(); ++i)
        {
            if (used[i])
            {
                remap[i] = newIndex++;
                newVertices.push_back(m_vertices[i]);
            }
        }

        for (auto& index : m_indices)
        {
            index = remap[index];
        }

        m_vertices = std::move(newVertices);
    }

    Vector3D centroid() const
    {
        Vector3D centroid;
        for (size_t i = 0; i < m_indices.size(); i += 3)
        {
            const Vector3D &a = m_vertices[m_indices[i]];
            const Vector3D &b = m_vertices[m_indices[i + 1]];
            const Vector3D &c = m_vertices[m_indices[i + 2]];
            centroid += (a + b + c) / 3.0f;
        }
        return centroid / (m_indices.size() / 3);
    }

    Mesh& translate(const Vector3D &translation)
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
        for (size_t i = 0; i < m_indices.size(); i += 3)
        {
            const Vector3D &a = m_vertices[m_indices[i]];
            const Vector3D &b = m_vertices[m_indices[i + 1]];
            const Vector3D &c = m_vertices[m_indices[i + 2]];

            if (Utility::Math::intersects(line, Triangle3D(a, b, c)))
                return true;
        }
        return false;
    }
};

} // namespace Math

} // namespace Utility