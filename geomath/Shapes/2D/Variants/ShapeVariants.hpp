/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "../Circle2D.hpp"
#include "../Triangle2D.hpp"
#include "../Rectangle2D.hpp"
#include "../Polygon2D.hpp"
#include "../ConvexPolygon2D.hpp"

#include <variant>

namespace Arns::geomath
{


using PolygonalShapeVariant2D = std::variant<
    Triangle2D,
    Rectangle2D,
    Polygon2D,
    ConvexPolygon2D
>;

using FiniteShapeVariant2D = std::variant<
    Circle2D,
    Triangle2D,
    Rectangle2D,
    Polygon2D,
    ConvexPolygon2D
>;

using FiniteShapePtrVariant2D = std::variant<
    std::monostate,
    const Circle2D*,
    const Triangle2D*,
    const Rectangle2D*,
    const Polygon2D*,
    const ConvexPolygon2D*
>;

using ShapeVariant2D = FiniteShapeVariant2D;//Todo: Move this

using GeometryVariant2D = std::variant<
    //Vector2D,
    Segment2D,
    Circle2D,
    Triangle2D,
    Rectangle2D,
    Polygon2D,
    ConvexPolygon2D
>;

}