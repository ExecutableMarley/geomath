/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <format>

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

inline std::ostream& operator<<(std::ostream& os, const SlotID& slotID)
{
    return os << "SlotID(index=" << slotID.index << ", generation=" << slotID.generation << ")";
}


} // namespace Arns::geomath

template <>
struct std::formatter<Arns::geomath::SlotID>
{
    constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const Arns::geomath::SlotID& slotID, FormatContext& ctx)
    {
        return std::format_to(ctx.out(), "SlotID(index={}, generation={})", slotID.index, slotID.generation);
    }
};
