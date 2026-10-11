#include "../../../third_party/doctest.h"

#include "Shapes/2D/Algorithms2D/Triangulation/Triangulation.hpp"

using namespace Arns::geomath;

TEST_SUITE("Earcut triangulation")
{
    TEST_CASE("handles degenerate and minimal inputs")
    {
        SUBCASE("empty polygon")
        {
            const Polygon2D poly;
            CHECK(triangulate(poly).m_triangles.empty());
        }

        SUBCASE("fewer than 3 vertices")
        {
            const Polygon2D poly({{0, 0}, {1, 1}});
            CHECK(triangulate(poly).m_triangles.empty());
        }

        SUBCASE("zero-area / collinear points")
        {
            const Polygon2D poly({{0, 0}, {2, 0}, {4, 0}});
            CHECK(triangulate(poly).m_triangles.empty());
        }
    }

    TEST_CASE("handles basic convex polygons")
    {
        SUBCASE("simple triangle")
        {
            const Polygon2D poly({{0, 0}, {1, 0}, {0, 1}});
            auto mesh = triangulate(poly);

            CHECK(mesh.isValid());
            CHECK(mesh.m_vertices.size() == 3);
            CHECK(mesh.m_triangles.size() == 1);
        }

        SUBCASE("axis-aligned square")
        {
            // Square: (0,0) to (1,1) -> should be split into 2 triangles
            const Polygon2D poly({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
            auto mesh = triangulate(poly);

            CHECK(mesh.isValid());
            CHECK(mesh.m_vertices.size() == 4);
            CHECK(mesh.m_triangles.size() == 2);
        }
    }

    TEST_CASE("handles concave polygons (reflex vertices)")
    {
        SUBCASE("L-shaped polygon")
        {
            // An L-shape has a reflex (interior > 180 degree) vertex
            // Coordinates: (0,0), (2,0), (2,1), (1,1), (1,2), (0,2)
            const Polygon2D poly({
                {0.0, 0.0}, 
                {2.0, 0.0}, 
                {2.0, 1.0}, 
                {1.0, 1.0}, 
                {1.0, 2.0}, 
                {0.0, 2.0}
            });
            
            auto mesh = triangulate(poly);

            CHECK(mesh.isValid());
            // decomposes into 4 triangles (N - 2 rule)
            CHECK(mesh.m_triangles.size() == 4);
        }
    }

    TEST_CASE("handles duplicate and collinear vertices in valid polygons")
    {
        // Polygon with a redundant point along an edge
        const Polygon2D poly = {{{0, 0}, {2, 0}, {4, 0}, {4, 4}, {0, 4}}};
        const auto mesh = triangulate(poly);

        // decomposes into 3 triangle
        CHECK(mesh.m_triangles.size() == 3);
    }

    TEST_CASE("triangulates a larger concave polygon")
    {
        std::vector<Vector2D> vertices{{0,0}};
        vertices.reserve(65);
        for (int i = 1; i <= 64; ++i)
            vertices.push_back({static_cast<real_t>(i), static_cast<real_t>(i % 2 ? 4 : 8)});
        vertices.push_back({64, 0});
        Polygon2D poly(std::move(vertices));

        const auto mesh = triangulate(poly);
        CHECK(mesh.m_triangles.size() == poly.vertexCount() - 2);
    }

    TEST_CASE("winding order independence / robustness")
    {
        SUBCASE("counter-clockwise vs clockwise square")
        {
            const Polygon2D polyCCW({{0, 0}, {1, 0}, {1, 1}, {0, 1}});
            const Polygon2D polyCW({{0, 0}, {0, 1}, {1, 1}, {1, 0}});

            auto meshCCW = triangulate(polyCCW);
            auto meshCW = triangulate(polyCW);

            CHECK(meshCCW.isValid());
            CHECK(meshCW.isValid());
            CHECK(meshCCW.m_triangles.size() == meshCW.m_triangles.size());
        }
    }

    TEST_CASE("handles extreme coordinate scales and precision limits")
    {
        SUBCASE("large coordinates")
        {
            const real_t shift = real_t{1e6};
            const Polygon2D poly = {{{shift, shift}, {shift + 4, shift}, {shift + 4, shift + 4}, {shift, shift + 4}}};
            const auto mesh = triangulate(poly);
            CHECK(mesh.m_triangles.size() == 2);
        }

        SUBCASE("tiny coordinates")
        {
            const real_t scale = AbsEpsilon;
            const Polygon2D poly = {{{0, 0}, {4 * scale, 0}, {4 * scale, 4 * scale}, {0, 4 * scale}}};
            const auto mesh = triangulate(poly);
            CHECK(mesh.m_triangles.size() == 2);
        }
    }
}


/*
TEST_SUITE("Earcut triangulation")
{
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

    
}

*/