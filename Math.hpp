/*
 * MIT License
 *
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <math.h>

#include "CommonMath.hpp"

#include "Geometry/Vector2D.hpp"
#include "Geometry/Vector3D.hpp"
#include "Geometry/Vector4D.hpp"

#include "Algebra/EulerAngles.hpp"
#include "Algebra/Fraction.hpp"
#include "Algebra/Quaternion.hpp"
#include "Algebra/Matrices/IMatrix.hpp"
#include "Algebra/Matrices/Matrix.hpp"
#include "Algebra/Matrices/Matrix3x3.hpp"
#include "Algebra/Matrices/Matrix4x4.hpp"
#include "Algebra/Matrices/Operators.hpp"
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
#include "Shapes/2D/ConvexPolygon2D.hpp"
#include "Shapes/2D/TriangleMesh2D.hpp"

#include "Shapes/2D/Intersections2D/Intersections2D.hpp"

#include "Shapes/2D/Algorithms2D/Delaunay/DelaunayTriangulation.hpp"
#include "Shapes/2D/Structures/ShapeStore2D.hpp"
#include "Shapes/2D/Structures/ISpatialIndex2D.hpp"

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