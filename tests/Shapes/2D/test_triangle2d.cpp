#include "test_shape2D_utility.hpp"
#include "Shapes/2D/Triangle2D.hpp"

TEST_CASE("Triangle2D polygonal interface")
{
    Triangle2D tri({
        Vector2D{0, 0},
        Vector2D{1, 0},
        Vector2D{0, 1}
    });

    SUBCASE("Shape type is correct")
    {
        CHECK(tri.type() == SHAPE2D_TRIANGLE);
    }

    check_is_polygonal<Triangle2D>(tri, 3);
}

TEST_CASE("Triangle2D vertex access and geometry")
{
    Triangle2D tri({
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{2, 3}
    });

    SUBCASE("Vertex access via a(), b(), c()")
    {
        CHECK(tri.a() == Vector2D{0, 0});
        CHECK(tri.b() == Vector2D{4, 0});
        CHECK(tri.c() == Vector2D{2, 3});
    }

    SUBCASE("Vertex access via operator[]")
    {
        CHECK(tri[0] == Vector2D{0, 0});
        CHECK(tri[1] == Vector2D{4, 0});
        CHECK(tri[2] == Vector2D{2, 3});
    }

    SUBCASE("Vertex access via vertexAt()")
    {
        CHECK(tri.vertexAt(0) == Vector2D{0, 0});
        CHECK(tri.vertexAt(1) == Vector2D{4, 0});
        CHECK(tri.vertexAt(2) == Vector2D{2, 3});
    }

    SUBCASE("Vertex mutation via operator[]")
    {
        tri[0] = Vector2D{1, 1};
        CHECK(tri.a() == Vector2D{1, 1});
    }

    SUBCASE("Vertex mutation via vertexAt()")
    {
        tri.vertexAt(1) = Vector2D{3, 1};
        CHECK(tri.b() == Vector2D{3, 1});
    }

    SUBCASE("vertices() returns span of all vertices")
    {
        auto verts = tri.vertices();
        CHECK(verts.size() == 3);
        CHECK(verts[0] == Vector2D{0, 0});
        CHECK(verts[1] == Vector2D{4, 0});
        CHECK(verts[2] == Vector2D{2, 3});
    }

    SUBCASE("Area calculation")
    {
        CHECK(tri.area() == doctest::Approx(6.0f));
    }

    SUBCASE("Perimeter calculation")
    {
        real_t expected = 4.0f + 2.0f * std::sqrt(13.0f);
        CHECK(tri.perimeter() == doctest::Approx(expected));
    }

    SUBCASE("Centroid is average of vertices")
    {
        CHECK(tri.centroid() == Vector2D{2, 1});
    }
}

TEST_CASE("Triangle2D transformations")
{
    Triangle2D tri({
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{2, 3}
    });

    SUBCASE("Translate moves all vertices")
    {
        tri.translate(Vector2D{1, -2});

        CHECK(tri.a() == Vector2D{1, -2});
        CHECK(tri.b() == Vector2D{5, -2});
        CHECK(tri.c() == Vector2D{3, 1});
    }

    SUBCASE("Rotate around centroid")
    {
        Vector2D originalCentroid = tri.centroid();
        tri.rotate(real_t{90});

        CHECK(tri.centroid() == originalCentroid);
        CHECK(tri.area() == doctest::Approx(6.0f));
    }

    SUBCASE("Rotate around arbitrary point")
    {
        Vector2D pivot{0, 0};
        tri.rotate(real_t{90}, pivot);

        CHECK(tri.a().x == doctest::Approx(0.0f));
        CHECK(tri.a().y == doctest::Approx(0.0f));
        CHECK(tri.area() == doctest::Approx(6.0f));
    }

    SUBCASE("Scale around centroid")
    {
        Vector2D originalCentroid = tri.centroid();
        tri.scale(real_t{2.0});

        CHECK(tri.centroid() == originalCentroid);
        CHECK(tri.area() == doctest::Approx(24.0f));
    }

    SUBCASE("Scale by 0.5")
    {
        tri.scale(real_t{0.5});

        CHECK(tri.area() == doctest::Approx(1.5f));
    }
}

TEST_CASE("Triangle2D copy and clone")
{
    Triangle2D original({
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{2, 3}
    });

    SUBCASE("Copy creates identical triangle")
    {
        Triangle2D copy = original.copy();

        CHECK(copy.a() == original.a());
        CHECK(copy.b() == original.b());
        CHECK(copy.c() == original.c());
    }

    SUBCASE("Clone returns unique_ptr IFiniteShape2D")
    {
        auto clone = original.clone();
        REQUIRE(clone != nullptr);

        CHECK(clone->type() == original.type());

        const Triangle2D* clonedTri = dynamic_cast<const Triangle2D*>(clone.get());
        REQUIRE(clonedTri != nullptr);

        CHECK(clonedTri->a() == original.a());
        CHECK(clonedTri->b() == original.b());
        CHECK(clonedTri->c() == original.c());
    }
}

TEST_CASE("Triangle2D bounding box")
{
    Triangle2D tri({
        Vector2D{1, 2},
        Vector2D{5, 1},
        Vector2D{3, 6}
    });

    SUBCASE("Bounding box encompasses all vertices")
    {
        BBox2D bbox = tri.boundingBox();

        CHECK(bbox.m_min == Vector2D{1, 1});
        CHECK(bbox.m_max == Vector2D{5, 6});
    }

    SUBCASE("Bounding box contains centroid")
    {
        BBox2D bbox = tri.boundingBox();
        Vector2D centroid = tri.centroid();

        CHECK(bbox.contains(centroid));
    }

    SUBCASE("Bounding box after translation")
    {
        tri.translate(Vector2D{-1, 3});
        BBox2D bbox = tri.boundingBox();

        CHECK(bbox.m_min == Vector2D{0, 4});
        CHECK(bbox.m_max == Vector2D{4, 9});
    }
}