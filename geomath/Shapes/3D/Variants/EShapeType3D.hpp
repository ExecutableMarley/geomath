/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

namespace Arns::geomath
{


enum ShapeType3D
{
    SHAPE3D_TRIANGLE,
    SHAPE3D_SPHERE,
    SHAPE3D_CYLINDER,
    SHAPE3D_CAPSULE,
    Count3D
};

constexpr std::string_view to_string(ShapeType3D type)
{
    switch (type)
    {
        case SHAPE3D_TRIANGLE: return "SHAPE3D_TRIANGLE";
        case SHAPE3D_SPHERE: return "SHAPE3D_SPHERE";
        case SHAPE3D_CYLINDER: return "SHAPE3D_CYLINDER";
        case SHAPE3D_CAPSULE: return "SHAPE3D_CAPSULE";
        default: return "UnknownShapeType3D";
    }
}

inline std::ostream& operator<<(std::ostream& os, ShapeType3D type)
{
    os << to_string(type);
    return os;
}


}