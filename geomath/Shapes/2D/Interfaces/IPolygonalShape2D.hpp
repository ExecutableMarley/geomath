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

    virtual bool isConvex() const = 0;

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




} // namespace Arns::geomath