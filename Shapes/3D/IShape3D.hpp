/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector3D.hpp"
#include "BBox3D.hpp"
#include "Ray3D.hpp"

namespace Utility
{

namespace Math
{

enum ShapeType3D
{
    SPHERE,
    CYLINDER,
    CAPSULE,
    BOX,
    TRIANGLE,
    MESH
};

/*
class IShape3D
{
public:
    virtual ShapeType3D type() const = 0;

    virtual float volume() const = 0;

    virtual float surfaceArea() const = 0;

    virtual Vector3D centroid() const = 0;

    virtual IShape3D& translate(const Vector3D &translation) = 0;
    
    virtual bool contains(const Vector3D &point) const = 0;

    //Todo: Implement
    //virtual bool intersects(const IShape3D &shape) const = 0;

    virtual BBox3D boundingBox() const = 0;
};*/

//Todo: Implement this
class IBaseShape3D
{
    virtual ShapeType3D type() const = 0;

    virtual IBaseShape3D& translate(const Vector3D &translation) = 0;

    virtual bool contains(const Vector3D &point) const = 0;
};

class IFiniteShape3D : IBaseShape3D
{
    virtual float volume() const = 0;

    virtual float surfaceArea() const = 0;

    virtual Vector3D centroid() const = 0;

    virtual BBox3D boundingBox() const = 0;
};
//

} // namespace Math

} // namespace Utility