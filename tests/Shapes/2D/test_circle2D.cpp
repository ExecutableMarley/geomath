#include "test_shape2D_utility.hpp"
#include "Shapes/2D/Circle2D.hpp"

TEST_CASE("Circle2D is not polygonal")
{
    Circle2D circle(Vector2D{0, 0}, real_t{1.0});

    SUBCASE("Shape type is correct")
    {
        CHECK(circle.type() == SHAPE2D_CIRCLE);
    }

    check_is_not_polygonal(circle);
}

TEST_CASE("Circle2D geometry and transformations")
{
    Circle2D circle(Vector2D{1, 2}, real_t{3.0});

    SUBCASE("Area, circumference, and perimeter")
    {
        CHECK(circle.area() == doctest::Approx(PI * real_t{9}));
        CHECK(circle.circumference() == doctest::Approx(real_t{2} * PI * real_t{3}));
        CHECK(circle.perimeter() == doctest::Approx(circle.circumference()));
    }

    SUBCASE("Centroid is center")
    {
        CHECK(circle.centroid() == Vector2D{1, 2});
    }

    SUBCASE("Get vertices returns requested resolution")
    {
        auto vertices = circle.getVertices(8);

        CHECK(vertices.size() == 8);
        CHECK(vertices[0].x == doctest::Approx(real_t{4}));
        CHECK(vertices[0].y == doctest::Approx(real_t{2}));
    }

    SUBCASE("Translate moves center without changing radius")
    {
        circle.translate(Vector2D{2, -1});

        CHECK(circle.centroid() == Vector2D{3, 1});
        CHECK(circle.m_radius == doctest::Approx(real_t{3}));
    }

    SUBCASE("Scale changes radius but keeps center")
    {
        circle.scale(real_t{0.5});

        CHECK(circle.centroid() == Vector2D{1, 2});
        CHECK(circle.m_radius == doctest::Approx(real_t{1.5}));
    }

    SUBCASE("Rotate around another point moves the center")
    {
        Circle2D rotatedCircle(circle);
        rotatedCircle.rotate(real_t{90}, Vector2D{0, 2});

        CHECK(rotatedCircle.centroid().x == doctest::Approx(real_t{0}));
        CHECK(rotatedCircle.centroid().y == doctest::Approx(real_t{3}));
        CHECK(rotatedCircle.m_radius == doctest::Approx(real_t{3}));
    }

    SUBCASE("Copy duplicates the circle")
    {
        Circle2D copy = circle.copy();

        CHECK(copy.centroid() == circle.centroid());
        CHECK(copy.m_radius == doctest::Approx(circle.m_radius));
    }

    SUBCASE("Clone returns an equivalent unique_ptr")
    {
        auto clone = circle.clone();
        REQUIRE(clone != nullptr);

        CHECK(clone->type() == circle.type());

        const Circle2D* clonedCircle = dynamic_cast<const Circle2D*>(clone.get());
        REQUIRE(clonedCircle != nullptr);

        CHECK(clonedCircle->centroid() == circle.centroid());
        CHECK(clonedCircle->m_radius == doctest::Approx(circle.m_radius));
    }

    SUBCASE("Bounding box matches circle extents")
    {
        BBox2D bbox = circle.boundingBox();

        CHECK(bbox.m_min == Vector2D{-2, -1});
        CHECK(bbox.m_max == Vector2D{4 ,  5});
        CHECK(bbox.width()  == doctest::Approx(real_t{6}));
        CHECK(bbox.height() == doctest::Approx(real_t{6}));
        CHECK(bbox.centroid() == circle.centroid());
    }
}