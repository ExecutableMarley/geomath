#include "../../third_party/doctest.h"

#include "Shapes/3D/TriangleMesh3D.hpp"

using namespace Arns::geomath;
using namespace Arns::geomath::literals;

TEST_SUITE("TriangleMesh3D")
{
    TEST_CASE("Construction and access")
    {
        TriangleMesh3D mesh(
            {{0, 0, 0}, {4, 0, 0}, {0, 3, 0}},
            {{0, 1, 2}});

        CHECK(mesh.vertices().size() == 3);
        CHECK(mesh.triangles().size() == 1);
        CHECK(mesh.vertices()[1] == Vector3D{4, 0, 0});
        CHECK(mesh[0].v0 == 0);
        CHECK(mesh[0].v1 == 1);
        CHECK(mesh[0].v2 == 2);
        CHECK(mesh.getTriangle(0).m_a == Vector3D{0, 0, 0});
        CHECK(mesh.getTriangle(0).m_b == Vector3D{4, 0, 0});
        CHECK(mesh.getTriangle(0).m_c == Vector3D{0, 3, 0});
        CHECK(mesh.isValid());
    }

    TEST_CASE("Data validation rejects bad indices and degenerate triangles")
    {
        CHECK_FALSE(TriangleMeshData3D({{0, 0, 0}, {1, 0, 0}}, {{0, 1, 2}}).isValid());
        CHECK_FALSE(TriangleMeshData3D({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, {{0, 1, 2}}).isValid());
        CHECK_FALSE(TriangleMeshData3D({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, {{0, 0, 2}}).isValid());
        CHECK(TriangleMeshData3D({{0, 0, 0}, {0, 1, 0}, {1, 0, 0}}, {{0, 1, 2}}).isValid());
    }

    TEST_CASE("Derived geometry")
    {
        TriangleMesh3D mesh(
            {{0, 0, 0}, {4, 0, 0}, {0, 3, 0}},
            {{0, 1, 2}});

        CHECK(mesh.localSurfaceArea() == doctest::Approx(real_t{6}));
        CHECK(mesh.worldSurfaceArea() == doctest::Approx(real_t{6}));
        CHECK(mesh.localCentroid() == Vector3D{real_t{4} / 3, 1, 0});
        CHECK(mesh.worldCentroid() == mesh.localCentroid());

        const BBox3D bounds = mesh.localBounds();
        CHECK(bounds.m_min == Vector3D{0, 0, 0});
        CHECK(bounds.m_max == Vector3D{4, 3, 0});
    }

    TEST_CASE("World geometry follows transforms")
    {
        TriangleMesh3D mesh(
            {{0, 0, 0}, {4, 0, 0}, {0, 3, 0}},
            {{0, 1, 2}});

        mesh.translate(Vector3D{2, -1, 5}).rotate(0_deg, 0_deg, 90_deg);
        const BBox3D bounds = mesh.worldBounds();
        CHECK(bounds.m_min == Vector3D{-1, -1, 5});
        CHECK(bounds.m_max == Vector3D{2, 3, 5});
        CHECK(mesh.worldCentroid() == Vector3D{1, real_t{1} / 3, 5});
        CHECK(mesh.worldSurfaceArea() == doctest::Approx(real_t{6}));

        const Triangle3D triangle = mesh.getTriangleWorld(0);
        CHECK(triangle.m_a == Vector3D{2, -1, 5});
        CHECK(triangle.m_b == Vector3D{2, 3, 5});
        CHECK(triangle.m_c == Vector3D{-1, -1, 5});
    }

    TEST_CASE("Scaling around a pivot and baking transforms")
    {
        TriangleMesh3D mesh(
            {{0, 0, 0}, {4, 0, 0}, {0, 3, 0}},
            {{0, 1, 2}});
        const auto localVertices = mesh.vertices();
        const Vector3D localCentroid = mesh.localCentroid();

        mesh.translate(Vector3D{2, -1, 1});
        mesh.scale(real_t{2}, Vector3D{0, 0, 0});

        CHECK(mesh.vertices() == localVertices);
        CHECK(mesh.localCentroid() == localCentroid);
        CHECK(mesh.worldCentroid() == Vector3D{real_t{20} / 3, 0, 2});
        CHECK(mesh.worldSurfaceArea() == doctest::Approx(real_t{24}));

        const Triangle3D worldTriangle = mesh.getTriangleWorld(0);
        mesh.bakeTransform();
        CHECK(mesh.getTriangle(0).m_a == worldTriangle.m_a);
        CHECK(mesh.getTriangle(0).m_b == worldTriangle.m_b);
        CHECK(mesh.getTriangle(0).m_c == worldTriangle.m_c);
        CHECK(mesh.worldSurfaceArea() == doctest::Approx(mesh.localSurfaceArea()));
        CHECK(mesh.worldCentroid() == mesh.localCentroid());
    }

    TEST_CASE("TriangleMeshBuilder3D builds mesh data")
    {
        TriangleMeshBuilder3D builder;
        builder.reserve(3, 1)
            .addVertex(0, 0, 0)
            .addVertex(4, 0, 0)
            .addVertex(0, 3, 0)
            .addTriangle(0, 1, 2);

        TriangleMesh3D mesh(std::move(builder).build());
        CHECK(mesh.vertices().size() == 3);
        CHECK(mesh.triangles().size() == 1);
        CHECK(mesh.localSurfaceArea() == doctest::Approx(real_t{6}));
    }

    TEST_CASE("TriangleMeshBuilder3D appends and merges mesh data")
    {
        TriangleMeshBuilder3D builder;
        builder.addVertices({{0, 0, 0}, {2, 0, 0}, {0, 2, 0}})
            .addTriangles({{0, 1, 2}});

        const TriangleMeshData3D other(
            {{3, 0, 0}, {5, 0, 0}, {3, 2, 0}},
            {{0, 1, 2}});
        builder.merge(other);

        const TriangleMeshData3D& data = builder.data();
        REQUIRE(data.m_vertices.size() == 6);
        REQUIRE(data.m_triangles.size() == 2);
        CHECK(data.m_triangles[0].v0 == 0);
        CHECK(data.m_triangles[0].v1 == 1);
        CHECK(data.m_triangles[0].v2 == 2);
        CHECK(data.m_triangles[1].v0 == 3);
        CHECK(data.m_triangles[1].v1 == 4);
        CHECK(data.m_triangles[1].v2 == 5);
        CHECK(data.isValid());
    }

    TEST_CASE("TriangleMeshBuilder3D filters triangles using their vertices")
    {
        TriangleMeshBuilder3D builder(TriangleMeshData3D(
            {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {3, 0, 0}, {5, 0, 0}, {3, 2, 0}},
            {{0, 1, 2}, {3, 4, 5}}));

        builder.filterTriangles([](const Vector3D& a, const Vector3D& b, const Vector3D& c) {
            return (a.x + b.x + c.x) / real_t{3} > real_t{2};
        });

        const TriangleMeshData3D& data = builder.data();
        REQUIRE(data.m_triangles.size() == 1);
        CHECK(data.m_triangles[0].v0 == 0);
        CHECK(data.m_triangles[0].v1 == 1);
        CHECK(data.m_triangles[0].v2 == 2);
        CHECK(data.m_vertices.size() == 6);
    }

    TEST_CASE("TriangleMeshBuilder3D transforms every vertex")
    {
        TriangleMeshBuilder3D builder(TriangleMeshData3D(
            {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}},
            {{0, 1, 2}}));

        builder.transformVertices([](Vector3D& vertex) {
            vertex += Vector3D{3, -1, 2};
        });

        const TriangleMeshData3D& data = builder.data();
        CHECK(data.m_vertices[0] == Vector3D{3, -1, 2});
        CHECK(data.m_vertices[1] == Vector3D{5, -1, 2});
        CHECK(data.m_vertices[2] == Vector3D{3, 1, 2});
        CHECK(data.isValid());
    }

    TEST_CASE("TriangleMeshBuilder3D removes orphan vertices and remaps indices")
    {
        TriangleMeshBuilder3D builder(TriangleMeshData3D(
            {{-10, -10, 0}, {0, 0, 0}, {2, 0, 0}, {99, 99, 0}, {0, 2, 0}, {2, 2, 0}},
            {{1, 2, 4}, {2, 5, 4}}));

        builder.removeUnusedVertices();

        const TriangleMeshData3D& data = builder.data();
        REQUIRE(data.m_vertices.size() == 4);
        REQUIRE(data.m_triangles.size() == 2);
        CHECK(data.m_vertices[0] == Vector3D{0, 0, 0});
        CHECK(data.m_vertices[1] == Vector3D{2, 0, 0});
        CHECK(data.m_vertices[2] == Vector3D{0, 2, 0});
        CHECK(data.m_vertices[3] == Vector3D{2, 2, 0});
        CHECK(data.m_triangles[0].v0 == 0);
        CHECK(data.m_triangles[0].v1 == 1);
        CHECK(data.m_triangles[0].v2 == 2);
        CHECK(data.m_triangles[1].v0 == 1);
        CHECK(data.m_triangles[1].v1 == 3);
        CHECK(data.m_triangles[1].v2 == 2);
        CHECK(data.isValid());

        TriangleMeshBuilder3D emptyAfterRemoval(TriangleMeshData3D(
            {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, {}));
        emptyAfterRemoval.removeUnusedVertices();
        CHECK(emptyAfterRemoval.data().m_vertices.empty());
    }
}