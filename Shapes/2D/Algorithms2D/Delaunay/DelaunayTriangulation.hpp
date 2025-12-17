#pragma once

#include <vector>

#include "Geometry/Vector2D.hpp"
#include "Shapes/2D/TriangleMesh2D.hpp"

namespace Arns
{

namespace Math
{

//[Predicates]

//[DelaunayAlgorithm]

bool isDelaunay(const TriangleMesh2D& mesh);

TriangleMesh2D fastDelaunayTriangulation(const std::vector<Vector2D>& points);


} // namespace Math

} // namespace Arns