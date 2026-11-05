/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Geometry/Vector3D.hpp"

namespace Arns::geomath
{

class Segment3D;

namespace detail
{

real_t distancePointToSegment(const Vector3D& point, const Vector3D& lineStart, const Vector3D& lineEnd, Vector3D* closestPoint = nullptr);
real_t distancePointToSegment(const Vector3D& point, const Segment3D& line, Vector3D* closestPoint = nullptr);

real_t distanceLineToLine(const Vector3D& s1, const Vector3D& s2, const Vector3D& k1, const Vector3D& k2,
                          Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);
real_t distanceLineToLine(const Segment3D& line1, const Segment3D& line2,
                          Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);

real_t distanceSegmentToSegment(const Vector3D& s1, const Vector3D& s2,
                                const Vector3D& k1, const Vector3D& k2,
                                Vector3D* closestPoint1 = nullptr, Vector3D* closestPoint2 = nullptr);

}

} // namespace Arns::geomath