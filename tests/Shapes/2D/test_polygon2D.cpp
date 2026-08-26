#include "test_shape2D_utility.hpp"
#include "Shapes/2D/Polygon2D.hpp"

TEST_CASE("Polygon2D polygonal interface")
{
    Polygon2D poly({
        Vector2D{0, 0},
        Vector2D{2, 0},
        Vector2D{1, 1},
        Vector2D{2, 2},
        Vector2D{0, 2}
    });

    SUBCASE("Shape type is correct")
    {
        CHECK(poly.type() == SHAPE2D_POLYGON);
    }

    check_is_polygonal<Polygon2D>(poly, 5);
}

TEST_CASE("Polygon2D vertex access and basic properties")
{
    Polygon2D poly({
        Vector2D{0, 0},
        Vector2D{2, 0},
        Vector2D{1, 1},
        Vector2D{2, 2},
        Vector2D{0, 2}
    });

    SUBCASE("Vertex count")
    {
        CHECK(poly.vertexCount() == 5);
    }

    SUBCASE("Vertex access via operator[]")
    {
        CHECK(poly[0] == Vector2D{0, 0});
        CHECK(poly[1] == Vector2D{2, 0});
        CHECK(poly[4] == Vector2D{0, 2});
    }

    SUBCASE("Vertex access via vertexAt()")
    {
        CHECK(poly.vertexAt(0) == Vector2D{0, 0});
        CHECK(poly.vertexAt(2) == Vector2D{1, 1});
    }

    SUBCASE("Wrapped vertex access")
    {
        CHECK(poly.wrappedVertexAt(0) == Vector2D{0, 0});
        CHECK(poly.wrappedVertexAt(5) == Vector2D{0, 0});
        CHECK(poly.wrappedVertexAt(-1) == Vector2D{0, 2});
    }

    SUBCASE("vertices() returns span")
    {
        auto verts = poly.vertices();
        CHECK(verts.size() == 5);
        CHECK(verts[0] == Vector2D{0, 0});
    }

    SUBCASE("Vertex mutation")
    {
        poly[0] = Vector2D{1, 1};
        CHECK(poly.vertexAt(0) == Vector2D{1, 1});
    }

    SUBCASE("Iterators")
    {
        size_t count = 0;
        for (auto& v : poly)
        {
            count++;
        }
        CHECK(count == 5);
    }
}

TEST_CASE("Polygon2D convexity and geometry")
{
    SUBCASE("Convex polygon detection - square")
    {
        Polygon2D square({
            Vector2D{0, 0},
            Vector2D{1, 0},
            Vector2D{1, 1},
            Vector2D{0, 1}
        });
        CHECK(square.isConvex());
    }

    SUBCASE("Concave polygon detection")
    {
        Polygon2D concave({
            Vector2D{0, 0},
            Vector2D{2, 0},
            Vector2D{1, 1},
            Vector2D{2, 2},
            Vector2D{0, 2}
        });
        CHECK(!concave.isConvex());
    }

    SUBCASE("Area calculation - triangle")
    {
        Polygon2D tri({
            Vector2D{0, 0},
            Vector2D{2, 0},
            Vector2D{0, 2}
        });
        CHECK(tri.area() == doctest::Approx(2.0f));
    }

    SUBCASE("Area calculation - square")
    {
        Polygon2D square({
            Vector2D{0, 0},
            Vector2D{2, 0},
            Vector2D{2, 2},
            Vector2D{0, 2}
        });
        CHECK(square.area() == doctest::Approx(4.0f));
    }

    SUBCASE("Perimeter calculation - square")
    {
        Polygon2D square({
            Vector2D{0, 0},
            Vector2D{3, 0},
            Vector2D{3, 3},
            Vector2D{0, 3}
        });
        CHECK(square.perimeter() == doctest::Approx(12.0f));
    }

    SUBCASE("Centroid calculation")
    {
        Polygon2D square({
            Vector2D{0, 0},
            Vector2D{2, 0},
            Vector2D{2, 2},
            Vector2D{0, 2}
        });
        CHECK(square.centroid() == Vector2D{1, 1});
    }
}

TEST_CASE("Polygon2D transformations")
{
    Polygon2D poly({
        Vector2D{0, 0},
        Vector2D{2, 0},
        Vector2D{2, 2},
        Vector2D{0, 2}
    });

    SUBCASE("Translate")
    {
        poly.translate(Vector2D{1, -1});
        CHECK(poly[0] == Vector2D{1, -1});
        CHECK(poly[2] == Vector2D{3, 1});
    }

    SUBCASE("Rotate around centroid preserves area")
    {
        real_t originalArea = poly.area();
        poly.rotate(real_t{45});
        CHECK(poly.area() == doctest::Approx(originalArea));
    }

    SUBCASE("Rotate around arbitrary point")
    {
        real_t originalArea = poly.area();
        poly.rotate(real_t{90}, Vector2D{0, 0});
        CHECK(poly.area() == doctest::Approx(originalArea));
    }

    SUBCASE("Scale around centroid")
    {
        Vector2D originalCentroid = poly.centroid();
        poly.scale(real_t{2.0});
        CHECK(poly.centroid() == originalCentroid);
        CHECK(poly.area() == doctest::Approx(16.0f));
    }

    SUBCASE("Scale around arbitrary point")
    {
        Vector2D pivot{0, 0};
        poly.scale(real_t{0.5}, pivot);
        CHECK(poly[0] == Vector2D{0, 0});
        CHECK(poly[1] == Vector2D{1, 0});
        CHECK(poly[2] == Vector2D{1, 1});
        CHECK(poly[3] == Vector2D{0, 1});
        CHECK(poly.area() == doctest::Approx(1.0));
    }
}

TEST_CASE("Polygon2D copy and clone")
{
    Polygon2D original({
        Vector2D{0, 0},
        Vector2D{2, 0},
        Vector2D{2, 2},
        Vector2D{0, 2}
    });

    SUBCASE("Copy creates identical polygon")
    {
        Polygon2D copy = original.copy();
        CHECK(copy.vertexCount() == original.vertexCount());
        CHECK(copy[0] == original[0]);
        CHECK(copy.area() == doctest::Approx(original.area()));
    }

    SUBCASE("Clone returns unique_ptr IFiniteShape2D")
    {
        auto clone = original.clone();
        REQUIRE(clone != nullptr);
        CHECK(clone->type() == original.type());

        const Polygon2D* clonedPoly = dynamic_cast<const Polygon2D*>(clone.get());
        REQUIRE(clonedPoly != nullptr);
        CHECK(clonedPoly->vertexCount() == original.vertexCount());
        CHECK(clonedPoly->area() == doctest::Approx(original.area()));
    }
}

TEST_CASE("Polygon2D bounding box")
{
    Polygon2D poly({
        Vector2D{1, 2},
        Vector2D{4, 1},
        Vector2D{3, 4},
        Vector2D{2, 3}
    });

    SUBCASE("Bounding box encompasses all vertices")
    {
        BBox2D bbox = poly.boundingBox();
        CHECK(bbox.m_min == Vector2D{1, 1});
        CHECK(bbox.m_max == Vector2D{4, 4});
    }

    SUBCASE("Bounding box after translation")
    {
        poly.translate(Vector2D{-1, 2});
        BBox2D bbox = poly.boundingBox();
        CHECK(bbox.m_min == Vector2D{0, 3});
        CHECK(bbox.m_max == Vector2D{3, 6});
    }

    SUBCASE("Empty polygon bounding box")
    {
        Polygon2D emptyPoly;
        BBox2D bbox = emptyPoly.boundingBox();
        CHECK(bbox.m_min == Vector2D{0, 0});
        CHECK(bbox.m_max == Vector2D{0, 0});
    }
}
