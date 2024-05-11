/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector3D.hpp"
#include "../Line3D.hpp"
#include "../BBox3D.hpp"
#include "../Sphere.hpp"
#include "../Cylinder.hpp"
#include "../Capsule.hpp"
#include "../Triangle3D.hpp"

namespace Utility
{

namespace Math
{

bool intersects(const Line3D &line1, const Line3D &line2, Vector3D *intersection = nullptr);

bool intersects(const Line3D &line, const BBox3D &bbox, Vector3D *intersection = nullptr);

bool intersects(const Line3D &line, const Sphere &sphere, Vector3D *intersection = nullptr);

bool intersects(const Line3D &line, const Cylinder &cylinder, Vector3D *intersection = nullptr);

bool intersects(const Line3D &line, const Capsule &capsule, Vector3D *intersection = nullptr);

bool intersects(const Line3D &line, const Triangle3D &triangle, Vector3D *intersection = nullptr);

} // namespace Math

} // namespace Utility