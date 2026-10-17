/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <cstdint>
#include <limits>

namespace Arns::geomath
{


struct SlotID
{
    uint32_t index;
    uint32_t generation;

    static constexpr SlotID invalid() 
    {
        return SlotID{std::numeric_limits<uint32_t>::max(), 0};
    }

    bool operator==(const SlotID& other) const
    {
        return index == other.index && generation == other.generation;
    }
};

using ShapeID = SlotID;


}