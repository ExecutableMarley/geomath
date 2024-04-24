/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Vector3D.hpp"

namespace Utility
{

namespace Math
{

class Line3D;

// Distance based algorithms

float distancePointToLine(const Vector3D& point, const Vector3D& lineStart, const Vector3D& lineEnd, Vector3D* closestPoint = nullptr);

float distancePointToLine(const Vector3D& point, const Line3D& line, Vector3D* closestPoint = nullptr);

float distanceLineToLine(const Vector3D& s1, const Vector3D& s2, const Vector3D& k1, const Vector3D& k2, Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);

float distanceLineToLine(const Line3D& line1, const Line3D& line2, Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);

} // namespace Math

} // namespace Utility