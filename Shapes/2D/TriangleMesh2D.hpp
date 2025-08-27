#pragma once

#include <math.h>
#include <stdexcept>
#include <vector>

#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "Triangle2D.hpp"
#include "IShape2D.hpp"

namespace Utility
{

namespace Math
{

struct TriangleIndices
{
    TriangleIndices(int v0, int v1, int v2) : v0(v0), v1(v1), v2(v2) {}

    int v0, v1, v2;

    unsigned int operator[](size_t index) const
    {
        switch (index)
        {
        case 1: return v0;
        case 2: return v1;
        case 3: return v2;
        default:
            return 0; //Throw exception
        }
    }
};

class TriangleMesh2D
{
public:
    std::vector <Vector2D> m_vertices;
    std::vector <TriangleIndices> m_triangles;

    TriangleMesh2D() : m_vertices(), m_triangles() {};

    TriangleMesh2D(const std::vector<Vector2D>& vertices, const std::vector<TriangleIndices>& triangles) :
        m_vertices(vertices), m_triangles(triangles) {}

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
};



} // namespace Math

} // namespace Utility