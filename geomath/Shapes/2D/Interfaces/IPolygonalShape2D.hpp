/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "IFiniteShape2D.hpp"

#include <span>

namespace Arns::geomath
{


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


} // namespace Arns::geomath