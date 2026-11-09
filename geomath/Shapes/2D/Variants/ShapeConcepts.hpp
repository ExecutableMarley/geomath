/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Formatting/Formatting.hpp"
#include "Geometry/Vector2D.hpp"
#include "VariantTraits.hpp"
#include "../Interfaces/EShapeType2D.hpp"
#include <concepts>
#include <span>

namespace Arns::geomath
{


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

template <typename T>
concept Shape2D = requires(T& t, const T& ct, const Vector2D& p) {
    //{ ct.contains(p) } -> std::same_as<bool>;
    //{ t.translate(p) } -> std::same_as<T&>;
    { ct.type() } -> std::same_as<ShapeType2D>;
};

template <typename T>
concept FiniteShape2D =
    //Shape2D<T> &&
    requires(const T& t) {
        { t.area() } -> std::same_as<real_t>;
        { t.perimeter() } -> std::same_as<real_t>;
        { t.centroid() } -> std::same_as<Vector2D>;
        { t.boundingBox() } -> std::same_as<BBox2D>;
    };

template <typename T>
concept PolygonalShape = 
    FiniteShape2D<T> && 
    requires(const T& t) {
    { t.vertices() } -> std::same_as<const std::span<const Vector2D>>;
    { t.vertexCount() } -> std::same_as<std::size_t>;
    { t.isConvex() } -> std::same_as<bool>;
};

template <typename T>
concept ValidShape2D = 
    std::is_same_v<T, Circle2D>    ||
    std::is_same_v<T, Rectangle2D>   ||
    std::is_same_v<T, Triangle2D>    ||
    std::is_same_v<T, Polygon2D>     ||
    std::is_same_v<T, ConvexPolygon2D>||
    is_variant_v<T>;


template <typename T>
struct PolygonalShapeFormatter
{
    int precision = 6;
    bool hasPrecision = false;

    constexpr auto parse(std::format_parse_context& ctx)
    {
        return parse_optional_float_format(ctx, precision, hasPrecision);
    }

    template <typename FormatContext>
    auto format(const T& shape, FormatContext& ctx) const
    {
        std::string verticesStr;

        for (size_t i = 0; i < shape.vertexCount(); ++i)
        {
            if (hasPrecision)
                verticesStr += std::format(
                    "[{:.{}f}, {:.{}f}]",
                    shape[i].x,
                    precision,
                    shape[i].y,
                    precision
                );
            else
                verticesStr += std::format("{}", shape[i]);

            if (i + 1 < shape.vertexCount())
                verticesStr += ", ";
        }

        return std::format_to(
            ctx.out(),
            "PolygonalShape2D(type: {}, vertices: [{}])",
            shape.type(),
            verticesStr
        );
    }
};

}