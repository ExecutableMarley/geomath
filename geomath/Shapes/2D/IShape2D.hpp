/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>
#include <span>
#include <memory>

#include "CommonMath.hpp"
#include "Formatting/Formatting.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"

namespace Arns
{

namespace geomath
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


//Possibly infinite
class IBaseShape2D
{
public:
    virtual ShapeType2D type() const = 0;

    virtual IBaseShape2D& translate(const Vector2D &translation) = 0;

    virtual bool contains(const Vector2D &point) const = 0;

    template <class T>
    const T* shape_cast() const
    {
        return (this->type() == T::shapeType) ? dynamic_cast<const T*>(this) : nullptr;
    }

    //virtual std::unique_ptr<IBaseShape2D> clone() const = 0;
};

class IPolygonalShape2D;

//Finite Shape
class IFiniteShape2D : public IBaseShape2D
{
public:
    virtual real_t area() const = 0;

    virtual real_t perimeter() const = 0;

    virtual Vector2D centroid() const = 0;

    virtual BBox2D boundingBox() const = 0;

    virtual std::unique_ptr<IFiniteShape2D> clone() const = 0;

    const IPolygonalShape2D* polygonal() const;
};

class IPolygonalShape2D : public IFiniteShape2D
{
public:
    virtual ~IPolygonalShape2D() = default;

    /// Number of vertices
    virtual size_t vertexCount() const = 0;

    /// Read-only access to all vertices
    virtual const std::span<const Vector2D> vertices() const = 0;

    /// Indexed vertex access
    virtual const Vector2D& operator[](size_t index) const = 0;
};

template <class T>
const T* shape_cast(const IBaseShape2D* shape)
{
    return (shape->type() == T::shapeType) ? dynamic_cast<const T*>(shape) : nullptr;
}

inline const IPolygonalShape2D *IFiniteShape2D::polygonal() const
{
    if (!isPolygonalShape(this->type()))
        return nullptr;

    return dynamic_cast<const IPolygonalShape2D *>(this);
}

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


} // namespace geomath

} // namespace Arns

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