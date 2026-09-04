#pragma once

#include <concepts>
#include <format>
#include <ostream>
#include <string>
#include <type_traits>

namespace Arns::Math
{


template <typename T>
concept Formattable = requires
{
    typename std::formatter<std::remove_cvref_t<T>, char>;
};

template <Formattable T>
std::string to_string(const T& value)
{
    return std::format("{}", value);
}

/*
template <Formattable T>
std::ostream& operator<<(std::ostream& stream, const T& value)
{
    return stream << to_string(value);
}
*/


template <typename ParseContext>
constexpr auto parse_optional_float_format(
    ParseContext& context,
    int& precision,
    bool& hasPrecision)
{
    auto iterator = context.begin();
    const auto end = context.end();
    precision = 6;
    hasPrecision = false;
    if (iterator != end && *iterator == '.')
    {
        ++iterator;
        hasPrecision = true;
        precision = 0;
        while (iterator != end &&
               *iterator >= '0' &&
               *iterator <= '9')
        {
            precision = precision * 10 + (*iterator - '0');
            ++iterator;
        }
    }
    if (iterator != end && *iterator == 'f')
        ++iterator;
    if (iterator != end && *iterator != '}')
        throw std::format_error("Invalid floating-point format specification");
    return iterator;
}


}