/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <memory>

#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"

namespace Utility
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

/*class IShape2D
{
public:
    virtual ShapeType2D type() const = 0;

    virtual float area() const = 0;

    virtual float perimeter() const = 0;

    virtual Vector2D centroid() const = 0;

    virtual IShape2D& translate(const Vector2D &translation) = 0;

    virtual bool contains(const Vector2D &point) const = 0;
    
    //Todo: Implement
    //virtual bool intersects(const IShape2D &shape) const = 0;

    virtual BBox2D boundingBox() const = 0;
};*/

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

class IFiniteShape2D : public IBaseShape2D
{
public:
    virtual float area() const = 0;

    virtual float perimeter() const = 0;

    virtual Vector2D centroid() const = 0;

    virtual BBox2D boundingBox() const = 0;

    virtual std::unique_ptr<IFiniteShape2D> clone() const = 0;
};

template <class T>
const T* shape_cast(const IBaseShape2D* shape)
{
    return (shape->type() == T::shapeType) ? dynamic_cast<const T*>(shape) : nullptr;
}

} // namespace Math

} // namespace Utility