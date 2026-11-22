/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Formatting/Formatting.hpp"
#include "Geometry/Vector3D.hpp"
#include "VariantTraits.hpp"
#include "EShapeType3D.hpp"
#include <concepts>
#include <span>

namespace Arns::geomath
{

class Matrix3x3;
class Ray3D;
class Segment3D;
class BBox3D;
class Triangle3D;
class Sphere3D;
class Cylinder3D;
class Capsule3D;

template <typename T>
concept Translatable3D = requires(T& t, const Vector3D& translation) {
    { t.translate(translation) } -> std::same_as<T&>;
};

template <typename T>
concept Rotatable3D = requires(T& t, const Matrix3x3& rotationMatrix, const Vector3D& point) {
    { t.rotate(rotationMatrix) } -> std::same_as<T&>;
    { t.rotateAround(rotationMatrix, point) } -> std::same_as<T&>;
};

template <typename T>
concept UniformlyScalable3D = requires(T& t, real_t scale) {
    { t.scale(scale) } -> std::same_as<T&>;
};

template <typename T>
concept Scalable3D = requires(T& t, real_t factor, const Vector3D& point) {
    { t.scale(factor) } -> std::same_as<T&>;
    { t.scale(factor, point) } -> std::same_as<T&>;
};

template <typename T>
concept Transformable3D = Translatable3D<T> && Rotatable3D<T> && Scalable3D<T>;

template <typename T>
concept PointContainable3D = requires(const T& t, const Vector3D& point) {
    { t.contains(point) } -> std::same_as<bool>;
};

template <typename T>
concept Bounded3D = requires(const T& t) {
    { t.boundingBox() } -> std::convertible_to<BBox3D>;
};

}