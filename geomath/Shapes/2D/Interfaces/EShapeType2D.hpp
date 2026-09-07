/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Formatting/Formatting.hpp"
#include <string_view>

namespace Arns::geomath
{


enum ShapeType2D
{
    SHAPE2D_TRIANGLE,
    SHAPE2D_RECTANGLE,
    SHAPE2D_CONVEX_POLYGON,
    SHAPE2D_POLYGON,
    SHAPE2D_CIRCLE
};

constexpr std::string_view to_string(ShapeType2D type)
{
    switch (type)
    {
        case SHAPE2D_TRIANGLE: return "Triangle2D";
        case SHAPE2D_RECTANGLE: return "Rectangle2D";
        case SHAPE2D_CONVEX_POLYGON: return "ConvexPolygon2D";
        case SHAPE2D_POLYGON: return "Polygon2D";
        case SHAPE2D_CIRCLE: return "Circle2D";
        default: return "UnknownShapeType2D";
    }
}

inline std::ostream& operator<<(std::ostream& os, ShapeType2D type)
{
    os << to_string(type);
    return os;
}

constexpr bool isPolygonalShape(ShapeType2D type)
{
    switch (type)
    {
        case SHAPE2D_TRIANGLE:
        case SHAPE2D_RECTANGLE:
        case SHAPE2D_CONVEX_POLYGON:
        case SHAPE2D_POLYGON:
            return true;
        default:
            return false;
    }
}

constexpr bool isConvexPolygonal(ShapeType2D type)
{
    switch (type)
    {
        case SHAPE2D_TRIANGLE:
        case SHAPE2D_RECTANGLE:
        case SHAPE2D_CONVEX_POLYGON:
            return true;
        default:
            return false;
    }
}


}


template <>
struct std::formatter<geomath::ShapeType2D> : std::formatter<std::string_view>
{
    constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const geomath::ShapeType2D& type, FormatContext& ctx) const
    {
        return std::format_to(ctx.out(), "{}", geomath::to_string(type));
    }
};