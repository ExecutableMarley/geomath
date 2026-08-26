#include "test_shape2D_utility.hpp"
#include "Shapes/2D/Rectangle2D.hpp"

TEST_CASE("Rectangle2D polygonal interface")
{
    Rectangle2D rect({
        Vector2D{0, 0},
        Vector2D{1, 0},
        Vector2D{1, 1},
        Vector2D{0, 1}
    });

    SUBCASE("Shape type is correct")
    {
        CHECK(rect.type() == SHAPE2D_RECTANGLE);
    }

    check_is_polygonal<Rectangle2D>(rect, 4);
}

TEST_CASE("Rectangle2D vertex access and geometry")
{
    Rectangle2D rect({
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 3},
        Vector2D{0, 3}
    });

    SUBCASE("Vertex access via a(), b(), c(), d()")
    {
        CHECK(rect.a() == Vector2D{0, 0});
        CHECK(rect.b() == Vector2D{4, 0});
        CHECK(rect.c() == Vector2D{4, 3});
        CHECK(rect.d() == Vector2D{0, 3});
    }

    SUBCASE("Vertex access via operator[]")
    {
        CHECK(rect[0] == Vector2D{0, 0});
        CHECK(rect[1] == Vector2D{4, 0});
        CHECK(rect[2] == Vector2D{4, 3});
        CHECK(rect[3] == Vector2D{0, 3});
    }

    SUBCASE("Vertex access via vertexAt()")
    {
        CHECK(rect.vertexAt(0) == Vector2D{0, 0});
        CHECK(rect.vertexAt(1) == Vector2D{4, 0});
        CHECK(rect.vertexAt(2) == Vector2D{4, 3});
        CHECK(rect.vertexAt(3) == Vector2D{0, 3});
    }

    SUBCASE("Vertex mutation")
    {
        rect[0] = Vector2D{1, 1};
        CHECK(rect.a() == Vector2D{1, 1});
    }

    SUBCASE("vertices() returns span of all vertices")
    {
        auto verts = rect.vertices();
        CHECK(verts.size() == 4);
    }

    SUBCASE("Width and height")
    {
        CHECK(rect.width() == doctest::Approx(4.0f));
        CHECK(rect.height() == doctest::Approx(3.0f));
    }

    SUBCASE("Area calculation")
    {
        CHECK(rect.area() == doctest::Approx(12.0f));
    }

    SUBCASE("Perimeter calculation")
    {
        CHECK(rect.perimeter() == doctest::Approx(14.0f));
    }

    SUBCASE("Centroid")
    {
        CHECK(rect.centroid() == Vector2D{2, 1.5f});
    }
}

TEST_CASE("Rectangle2D constructors")
{
    SUBCASE("Constructor from 4 vertices")
    {
        Rectangle2D rect({Vector2D{0, 0}, Vector2D{5, 0}, Vector2D{5, 3}, Vector2D{0, 3}});
        CHECK(rect.area() == doctest::Approx(15.0f));
    }

    SUBCASE("Constructor from position, width, height")
    {
        Rectangle2D rect(Vector2D{1, 2}, real_t{4}, real_t{3});
        CHECK(rect.a() == Vector2D{1, 2});
        CHECK(rect.b() == Vector2D{5, 2});
        CHECK(rect.c() == Vector2D{5, 5});
        CHECK(rect.d() == Vector2D{1, 5});
        CHECK(rect.area() == doctest::Approx(12.0f));
    }

    SUBCASE("fromMinMax static factory")
    {
        Rectangle2D rect = Rectangle2D::fromMinMax(Vector2D{1, 2}, Vector2D{5, 5});
        CHECK(rect.a() == Vector2D{1, 2});
        CHECK(rect.width() == doctest::Approx(4.0f));
        CHECK(rect.height() == doctest::Approx(3.0f));
    }
}

TEST_CASE("Rectangle2D transformations")
{
    Rectangle2D rect({
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 3},
        Vector2D{0, 3}
    });

    SUBCASE("Translate moves all vertices")
    {
        rect.translate(Vector2D{2, -1});
        CHECK(rect.a() == Vector2D{2, -1});
        CHECK(rect.area() == doctest::Approx(12.0f));
    }

    SUBCASE("Rotate around centroid")
    {
        Vector2D originalCentroid = rect.centroid();
        rect.rotate(real_t{90});
        CHECK(rect.centroid() == originalCentroid);
        CHECK(rect.area() == doctest::Approx(12.0f));
    }

    SUBCASE("Rotate around arbitrary point")
    {
        rect.rotate(real_t{90}, Vector2D{0, 0});
        CHECK(rect.area() == doctest::Approx(12.0f));
    }

    SUBCASE("Scale around centroid")
    {
        Vector2D originalCentroid = rect.centroid();
        rect.scale(real_t{2.0});
        CHECK(rect.centroid() == originalCentroid);
        CHECK(rect.area() == doctest::Approx(48.0f));
    }

    SUBCASE("Scale around arbitrary point")
    {
        rect.scale(real_t{0.5}, Vector2D{0, 0});
        CHECK(rect.a() == Vector2D{0, 0});
        CHECK(rect.b() == Vector2D{2, 0});
        CHECK(rect.c() == Vector2D{2, 1.5f});
        CHECK(rect.d() == Vector2D{0, 1.5f});
        CHECK(rect.area() == doctest::Approx(3.0f));
    }
}

TEST_CASE("Rectangle2D copy and clone")
{
    Rectangle2D original({
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 3},
        Vector2D{0, 3}
    });

    SUBCASE("Copy creates identical rectangle")
    {
        Rectangle2D copy = original.copy();
        CHECK(copy.a() == original.a());
        CHECK(copy.area() == doctest::Approx(original.area()));
    }

    SUBCASE("Clone returns unique_ptr IFiniteShape2D")
    {
        auto clone = original.clone();
        REQUIRE(clone != nullptr);
        CHECK(clone->type() == original.type());

        const Rectangle2D* clonedRect = dynamic_cast<const Rectangle2D*>(clone.get());
        REQUIRE(clonedRect != nullptr);
        CHECK(clonedRect->area() == doctest::Approx(original.area()));
    }
}

TEST_CASE("Rectangle2D bounding box")
{
    Rectangle2D rect({
        Vector2D{1, 2},
        Vector2D{5, 2},
        Vector2D{5, 5},
        Vector2D{1, 5}
    });

    SUBCASE("Bounding box encompasses all vertices")
    {
        BBox2D bbox = rect.boundingBox();
        CHECK(bbox.m_min == Vector2D{1, 2});
        CHECK(bbox.m_max == Vector2D{5, 5});
    }

    SUBCASE("Bounding box contains centroid")
    {
        BBox2D bbox = rect.boundingBox();
        CHECK(bbox.contains(rect.centroid()));
    }
}