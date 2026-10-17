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
    uint32_t index      = invalid_index;
    uint32_t generation = 0;

    static constexpr uint32_t invalid_index = std::numeric_limits<uint32_t>::max();

    static constexpr SlotID invalid() noexcept
    {
        return SlotID{std::numeric_limits<uint32_t>::max(), 0};
    }

    friend constexpr bool operator==(const SlotID&, const SlotID&) = default;
};

using ShapeID = SlotID;


}