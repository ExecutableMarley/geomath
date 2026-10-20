#pragma once

#include "../CommonMath.hpp"

#include <concepts>
#include <cstddef>

namespace Arns::geomath
{

template <typename Vector>
concept VectorLike = std::default_initializable<Vector> &&
    requires(Vector& vector, const Vector& constVector, std::size_t index)
    {
        { Vector::size() } -> std::convertible_to<std::size_t>;
        { vector[index] } -> std::same_as<real_t&>;
        { constVector[index] } -> std::same_as<const real_t&>;
    } && (Vector::size() > 0);

}