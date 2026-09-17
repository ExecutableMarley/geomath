#pragma once

#include <cstddef>
#include <stdexcept>
#include <cstdint>

namespace Arns::geomath
{

struct TriangleIndices
{
    TriangleIndices(std::uint32_t v0, std::uint32_t v1, std::uint32_t v2) : v0(v0), v1(v1), v2(v2) {}

    std::uint32_t v0, v1, v2;

    uint32_t operator[](size_t index) const
    {
        switch (index)
        {
        case 0: return v0;
        case 1: return v1;
        case 2: return v2;
        default: throw std::out_of_range("Triangle vertex index out of range");
        }
    }
};

}