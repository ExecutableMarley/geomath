/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once
#include <type_traits>
#include <variant>

namespace Arns::geomath
{


// 1. Detect if a type is any std::variant
template <typename T>
struct is_variant : std::false_type {};

template <typename... Ts>
struct is_variant<std::variant<Ts...>> : std::true_type {};

template <typename T>
inline constexpr bool is_variant_v = is_variant<T>::value;

// 2. Extract or query variant alternative types
template <typename T, typename Variant>
struct is_variant_member : std::false_type {};

template <typename T, typename... Ts>
struct is_variant_member<T, std::variant<Ts...>> 
    : std::disjunction<std::is_same<T, Ts>...> {};

template <typename T, typename Variant>
inline constexpr bool is_variant_member_v = is_variant_member<T, Variant>::value;


}