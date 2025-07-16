/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"

#include "Geometry/Vector2D.hpp"
#include "Geometry/Vector3D.hpp"
#include "Geometry/Vector4D.hpp"

#include "Algebra/EulerAngles.hpp"
#include "Algebra/Quaternion.hpp"
#include "Algebra/RotationMatrix.hpp"
#include "Algebra/ViewMatrix.hpp"

// 2D
#include "Shapes/2D/Algorithms2D/Algorithms2D.hpp"

#include "Shapes/2D/Ray2D.hpp"
#include "Shapes/2D/Line2D.hpp"
#include "Shapes/2D/BBox2D.hpp"
#include "Shapes/2D/IShape2D.hpp"
#include "Shapes/2D/Circle2D.hpp"
#include "Shapes/2D/Triangle2D.hpp"
#include "Shapes/2D/Rectangle2D.hpp"
#include "Shapes/2D/Polygon2D.hpp"
#include "Shapes/2D/ConvexHull2D.hpp"

#include "Shapes/2D/Intersections2D/Intersections2D.hpp"

// 3D
#include "Shapes/3D/Algorithms3D/Algorithms3D.hpp"

#include "Shapes/3D/Ray3D.hpp"
#include "Shapes/3D/Line3D.hpp"
#include "Shapes/3D/BBox3D.hpp"
#include "Shapes/3D/IShape3D.hpp"
#include "Shapes/3D/Sphere.hpp"
#include "Shapes/3D/Cylinder.hpp"
#include "Shapes/3D/Capsule.hpp"
#include "Shapes/3D/Triangle3D.hpp"

#include "Shapes/3D/Intersection/Intersections3D.hpp"