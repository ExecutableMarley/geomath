#pragma once

#include "../../TriangleMesh2D.hpp"


namespace Arns::geomath
{


bool isDelaunay(const TriangleMesh2D& mesh);


TriangleMesh2D delaunay(std::span<const Vector2D> points);

TriangleMesh2D triangulate(const Polygon2D& polygon);

bool refineToDelaunay(TriangleMesh2D& mesh);



}