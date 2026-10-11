#pragma once

#include "../../TriangleMesh2D.hpp"
#include "Detail/EarcutTriangulation.hpp"

namespace Arns::geomath
{

//[Predicates]

//[DelaunayAlgorithm]


// Defined inside DelaunayTriangulation.cpp

bool isDelaunay(const TriangleMeshData2D& mesh);

inline bool isDelaunay(const TriangleMesh2D& mesh) { return isDelaunay(mesh.data()); }

TriangleMeshData2D delaunay(const std::span<const Vector2D> points);


//Given by EarcutTriangulation.hpp and earcut.hpp

inline TriangleMeshData2D triangulate(std::span<const Vector2D> polygonVertices)
{
    auto triangles = detail::earcut(polygonVertices);
    return TriangleMeshData2D(polygonVertices, triangles);
}

template <PolygonalShape shape>
inline TriangleMeshData2D triangulate(const shape& polygon)
{
    auto triangles = detail::earcut(polygon.vertices());
    return TriangleMeshData2D(polygon.vertices(), triangles);
}

//bool refineToDelaunay(TriangleMesh2D& mesh);



}