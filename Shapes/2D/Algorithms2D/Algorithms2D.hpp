/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector2D.hpp"
#include "../Line2D.hpp"

namespace Utility
{

namespace Math
{

// Distance based algorithms

float distancePointToLine(const Vector2D& point, const Vector2D& lineStart, const Vector2D& lineEnd, Vector2D* closestPoint = nullptr);

float distancePointToLine(const Vector2D& point, const Line2D& line, Vector2D* closestPoint = nullptr);

float distanceLineToLine(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);

float distanceLineToLine(const Line2D& line1, const Line2D& line2, Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr);


} // namespace Math

} // namespace Utility