/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "../Line2D.hpp"
#include "../Ray2D.hpp"
#include "../Variants/ShapeConcepts.hpp"
#include "Containment2D.hpp"
#include "Distance2D.hpp"
#include "Intersection2D.hpp"
#include "CCD2D.hpp"

#include <math.h>
#include <vector>
#include <span>

namespace Arns
{

namespace geomath
{

// Forward declarations
class Ray2D;
class Segment2D;
class BBox2D;
class IBaseShape2D;
class IFiniteShape2D;
class IPolygonalShape2D;
class Circle2D;
class Triangle2D;
class Rectangle2D;
class Polygon2D;
class ConvexPolygon2D;

} // namespace Math

} // namespace Arns