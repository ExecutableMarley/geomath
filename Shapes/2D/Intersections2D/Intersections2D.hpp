/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector2D.hpp"
#include "../Line2D.hpp"
#include "../BBox2D.hpp"
#include "../Circle.hpp"
#include "../Polygon.hpp"

namespace Utility
{

namespace Math
{

//[Line-Line Intersection]

bool intersects(const Line2D& line1, const Line2D& line2, Vector2D* intersection = nullptr);


//[Line-Rectangle Intersection]

bool intersects(const Line2D& line, const BBox2D& rectangle, Vector2D* intersection = nullptr);


//[Line-Circle Intersection]

bool intersects(const Line2D& line, const Circle& circle, Vector2D* intersection = nullptr);


//[Line-Polygon Intersection]

bool intersects(const Line2D& line, const Polygon& polygon, Vector2D* intersection = nullptr);

} // namespace Math

} // namespace Utility