#include "test_shape2D_utility.hpp"

#include "Shapes/2D/TriangleMesh2D.hpp"

using namespace Arns::geomath;

TEST_SUITE("TriangleMesh2D")
{

TEST_CASE("Construction and access")
{
    TriangleMesh2D mesh(
        {{0, 0}, {4, 0}, {0, 3}},
        {{0, 1, 2}});

    CHECK(mesh.vertices().size() == 3);
    CHECK(mesh.triangles().size() == 1);
    CHECK(mesh.vertices()[1] == Vector2D{4, 0});
    CHECK(mesh[0].v0 == 0);
    CHECK(mesh[0].v1 == 1);
    CHECK(mesh[0].v2 == 2);
    CHECK(mesh.getTriangle(0).a() == Vector2D{0, 0});
    CHECK(mesh.getTriangle(0).b() == Vector2D{4, 0});
    CHECK(mesh.getTriangle(0).c() == Vector2D{0, 3});
    CHECK(mesh.isValid());
}

TEST_CASE("Derived geometry")
{
    TriangleMesh2D mesh(
        {{0, 0}, {4, 0}, {0, 3}},
        {{0, 1, 2}});

    CHECK(mesh.localArea() == doctest::Approx(real_t{6}));
    CHECK(mesh.worldArea() == doctest::Approx(real_t{6}));
    CHECK(mesh.localCentroid() == Vector2D{real_t{4} / 3, 1});
    CHECK(mesh.worldCentroid() == mesh.localCentroid());

    const BBox2D bounds = mesh.localBounds();
    CHECK(bounds.m_min == Vector2D{0, 0});
    CHECK(bounds.m_max == Vector2D{4, 3});
}

TEST_CASE("World geometry follows transforms")
{
    TriangleMesh2D mesh(
        {{0, 0}, {4, 0}, {0, 3}},
        {{0, 1, 2}});

    Triangle2D triangle = mesh.getTriangleWorld(0);
    CHECK(mesh.worldBounds() == BBox2D({0, 0}, {4, 3}));

    mesh.translate(Vector2D{2, -1});
    triangle.translate(Vector2D{2, -1});
    CHECK(mesh.worldBounds() == BBox2D({2, -1}, {6, 2}));
    CHECK(mesh.worldCentroid() == Vector2D{real_t{10} / 3, 0});

    mesh.scale(real_t{2});
    triangle.scale(real_t{2}, triangle.centroid());
    CHECK(mesh.getTriangleWorld(0) == triangle);
}

TEST_CASE("Scaling around a world pivot preserves relative geometry")
{
    TriangleMesh2D mesh(
        {{0, 0}, {4, 0}, {0, 3}},
        {{0, 1, 2}});

    const auto localVertices = mesh.vertices();
    const Vector2D localCentroid = mesh.localCentroid();
    const real_t localArea = mesh.localArea();

    mesh.translate(Vector2D{2, -1});
    mesh.scale(real_t{2}, Vector2D{0, 0});

    CHECK(mesh.vertices() == localVertices);
    CHECK(mesh.localCentroid() == localCentroid);
    CHECK(mesh.localArea() == localArea);

    const Triangle2D worldTriangle = mesh.getTriangleWorld(0);
    CHECK(worldTriangle.a() == Vector2D{4, -2});
    CHECK(worldTriangle.b() == Vector2D{12, -2});
    CHECK(worldTriangle.c() == Vector2D{4, 4});
}

TEST_CASE("Rejects invalid triangle access")
{
    TriangleMesh2D mesh(
        {{0, 0}, {1, 0}, {0, 1}},
        {{0, 1, 2}});

    CHECK_THROWS_AS(mesh.getTriangle(1), std::out_of_range);
}

TEST_CASE("TriangleMeshBuilder2D builds mesh data")
{
    TriangleMeshBuilder2D builder;
    builder.reserve(3, 1)
        .addVertex(0, 0)
        .addVertex(4, 0)
        .addVertex(0, 3)
        .addTriangle(0, 1, 2);

    const TriangleMeshData2D data = std::move(builder).build();
    TriangleMesh2D mesh(std::move(data));

    CHECK(mesh.vertices().size() == 3);
    CHECK(mesh.triangles().size() == 1);
    CHECK(mesh.localArea() == doctest::Approx(real_t{6}));
}

}