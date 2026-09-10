#include "../../../third_party/doctest.h"
#include "Shapes/2D/Algorithms2D/Triangulation/Detail/EarcutTriangulation.hpp"

using namespace Arns::geomath;

TEST_SUITE("Earcut triangulation")
{
    TEST_CASE("handles degenerate and minimal inputs")
    {
        SUBCASE("empty polygon")
        {
            const PolygonRings2D rings = {};
            CHECK(earcut(rings).empty());
        }

        SUBCASE("fewer than 3 vertices")
        {
            const PolygonRings2D rings = {{{0, 0}, {1, 1}}};
            CHECK(earcut(rings).empty());
        }

        SUBCASE("zero-area / collinear points")
        {
            const PolygonRings2D rings = {{{0, 0}, {2, 0}, {4, 0}}};
            CHECK(earcut(rings).empty());
        }
    }

    TEST_CASE("handles duplicate and collinear vertices in valid polygons")
    {
        // Polygon with a redundant point along an edge
        const PolygonRings2D rings = {{{0, 0}, {2, 0}, {4, 0}, {4, 4}, {0, 4}}};
        const auto triangles = earcut(rings);

        // Expected triangles: N - 2 = 3
        CHECK(triangles.size() == 3);
    }

    TEST_CASE("triangulates a concave polygon")
    {
        const PolygonRings2D rings = {{{0, 0}, {4, 0}, {4, 4}, {2, 2}, {0, 4}}};
        const auto triangles = earcut(rings);

        CHECK(triangles.size() == 3);
        real_t area = 0;
        for (const auto& triangle : triangles)
            area += std::abs((rings[0][triangle.v1] - rings[0][triangle.v0]).cross(
                rings[0][triangle.v2] - rings[0][triangle.v0])) / 2;
        CHECK(area == doctest::Approx(real_t{12}));
    }

    TEST_CASE("triangulates a polygon with a hole")
    {
        const PolygonRings2D rings = {
            {{0, 0}, {10, 0}, {10, 10}, {0, 10}},
            {{3, 3}, {3, 7}, {7, 7}, {7, 3}}
        };
        const auto mesh = earcutMesh(rings);

        CHECK(mesh.validate());
        real_t area = 0;
        for (const auto& triangle : mesh.m_triangles)
            area += std::abs((mesh.m_vertices[triangle.v1] - mesh.m_vertices[triangle.v0]).cross(
                mesh.m_vertices[triangle.v2] - mesh.m_vertices[triangle.v0])) / 2;
        CHECK(area == doctest::Approx(real_t{84}));
    }

    TEST_CASE("preserves concatenated indices for either winding order")
    {
        const PolygonRings2D rings = {
            {{0, 0}, {0, 6}, {6, 6}, {6, 0}},
            {{2, 2}, {4, 2}, {4, 4}, {2, 4}}
        };
        const auto triangles = earcut(rings);

        CHECK(triangles.size() == 8);
        for (const auto& triangle : triangles)
        {
            CHECK(triangle.v0 < 8);
            CHECK(triangle.v1 < 8);
            CHECK(triangle.v2 < 8);
        }
    }

    TEST_CASE("triangulates a larger concave polygon")
    {
        PolygonRings2D rings = {{{0, 0}}};
        for (int i = 1; i <= 64; ++i)
            rings[0].push_back({static_cast<real_t>(i), static_cast<real_t>(i % 2 ? 4 : 8)});
        rings[0].push_back({64, 0});

        const auto triangles = earcut(rings);
        CHECK(triangles.size() == rings[0].size() - 2);
    }

    TEST_CASE("triangulates a polygon with multiple internal holes")
    {
        const PolygonRings2D rings = {
            {{0, 0}, {12, 0}, {12, 12}, {0, 12}}, // 12x12 outer (area 144)
            {{1, 1}, {1, 4}, {4, 4}, {4, 1}},     // 3x3 hole (area 9)
            {{7, 7}, {7, 10}, {10, 10}, {10, 7}}  // 3x3 hole (area 9)
        };
        const auto mesh = earcutMesh(rings);

        CHECK(mesh.validate());

        real_t area = 0;
        for (const auto &triangle : mesh.m_triangles)
            area += std::abs((mesh.m_vertices[triangle.v1] - mesh.m_vertices[triangle.v0]).cross(mesh.m_vertices[triangle.v2] - mesh.m_vertices[triangle.v0])) / 2;

        // Expected area: 144 - 9 - 9 = 126
        CHECK(area == doctest::Approx(real_t{126}));
    }

    TEST_CASE("handles touching hole and outer boundary")
    {
        // Hole touches outer ring at vertex {0, 0}
        const PolygonRings2D rings = {
            {{0, 0}, {10, 0}, {10, 10}, {0, 10}},
            {{0, 0}, {2, 2}, {2, 4}, {0, 4}}};
        const auto mesh = earcutMesh(rings);
        CHECK(mesh.validate());
    }

    TEST_CASE("handles extreme coordinate scales and precision limits")
    {
        SUBCASE("large coordinates")
        {
            const real_t shift = real_t{1e6};
            const PolygonRings2D rings = {{{shift, shift}, {shift + 4, shift}, {shift + 4, shift + 4}, {shift, shift + 4}}};
            const auto triangles = earcut(rings);
            CHECK(triangles.size() == 2);
        }

        SUBCASE("tiny coordinates")
        {
            const real_t scale = AbsEpsilon;
            const PolygonRings2D rings = {{{0, 0}, {4 * scale, 0}, {4 * scale, 4 * scale}, {0, 4 * scale}}};
            const auto triangles = earcut(rings);
            CHECK(triangles.size() == 2);
        }
    }
}