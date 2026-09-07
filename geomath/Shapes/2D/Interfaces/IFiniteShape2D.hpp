/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "Formatting/Formatting.hpp"
#include "EShapeType2D.hpp"
#include "Geometry/Vector2D.hpp"
#include "../BBox2D.hpp"

#include <memory>

namespace Arns::geomath
{


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

template <class T>
const T* shape_cast(const IBaseShape2D* shape)
{
    return (shape->type() == T::shapeType) ? dynamic_cast<const T*>(shape) : nullptr;
}


}