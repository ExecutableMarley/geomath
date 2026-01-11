/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>
#include <memory>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"

namespace Arns
{

namespace Math
{

enum ShapeType2D
{
    SHAPE2D_TRIANGLE,
    SHAPE2D_RECTANGLE,
    SHAPE2D_CONVEX_POLYGON,
    SHAPE2D_POLYGON,
    SHAPE2D_CIRCLE
};

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

/*
Consider using std::span for efficient vertices access
*/

class IPolygonalShape2D
{
public:
    virtual ~IPolygonalShape2D() = default;

    /// Number of vertices
    virtual size_t vertexCount() const = 0;

    /// Read-only access to all vertices
    virtual std::vector<Vector2D> getVertices() const = 0;

    /// Indexed vertex access
    virtual const Vector2D& operator[](size_t index) const = 0;
};


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

//Finite Shape
class IFiniteShape2D : public IBaseShape2D
{
public:
    virtual real_t area() const = 0;

    virtual real_t perimeter() const = 0;

    virtual Vector2D centroid() const = 0;

    virtual BBox2D boundingBox() const = 0;

    virtual std::unique_ptr<IFiniteShape2D> clone() const = 0;

    const IPolygonalShape2D* polygonal() const
    {
        if (!isPolygonalShape(this->type()))
            return nullptr;

        return dynamic_cast<const IPolygonalShape2D*>(this);
    }
};

template <class T>
const T* shape_cast(const IBaseShape2D* shape)
{
    return (shape->type() == T::shapeType) ? dynamic_cast<const T*>(shape) : nullptr;
}


} // namespace Math

} // namespace Arns