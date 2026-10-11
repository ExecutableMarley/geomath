/*
 * Geomath adapter for mapbox/earcut.hpp.
 *
 * The complete upstream implementation is vendored at:
 * mapbox/earcut.hpp
 *
 * Copyright (c) 2015, Mapbox. ISC License.
 */

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <tuple>
#include <vector>

#include "Shapes/2D/TriangleMesh2D.hpp"
#include "Shapes/2D/Polygon2D.hpp"
#include "earcut.hpp"

namespace Arns::geomath::detail
{

using PolygonRings2D = std::vector<std::vector<Vector2D>>;

// Todo: Unecessary copies

inline std::vector<TriangleIndices> earcut(std::span<const Vector2D> simplePolygon)
{
    using EarcutPoint = std::array<real_t, 2>;

    std::vector<std::vector<EarcutPoint>> points(1);
    points.front().reserve(simplePolygon.size());
    for (const auto& point : simplePolygon)
        points.front().push_back({point.x, point.y});

    const auto indices = mapbox::earcut<std::uint32_t>(points);
    std::vector<TriangleIndices> triangles;
    triangles.reserve(indices.size() / 3);
    for (size_t i = 0; i + 2 < indices.size(); i += 3)
    {
        triangles.emplace_back(static_cast<int>(indices[i]),
                               static_cast<int>(indices[i + 1]),
                               static_cast<int>(indices[i + 2]));
    }
    return triangles;
}

inline std::vector<TriangleIndices> earcut(const PolygonRings2D& rings)
{
    using EarcutPoint = std::array<real_t, 2>;

    std::vector<std::vector<EarcutPoint>> points;
    points.reserve(rings.size());
    for (const auto& ring : rings)
    {
        auto& convertedRing = points.emplace_back();
        convertedRing.reserve(ring.size());
        for (const auto& point : ring)
            convertedRing.push_back({point.x, point.y});
    }

    const auto indices = mapbox::earcut<std::uint32_t>(points);
    std::vector<TriangleIndices> triangles;
    triangles.reserve(indices.size() / 3);
    for (size_t i = 0; i + 2 < indices.size(); i += 3)
    {
        triangles.emplace_back(static_cast<int>(indices[i]),
                               static_cast<int>(indices[i + 1]),
                               static_cast<int>(indices[i + 2]));
    }
    return triangles;
}

} // namespace Arns::Math::geomath
