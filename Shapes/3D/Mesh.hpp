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
        for (auto vector : m_vertices)
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