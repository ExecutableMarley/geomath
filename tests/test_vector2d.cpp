#include "third_party/doctest.h"

#include "Geometry/Vector2D.hpp"

using namespace Arns::Math;

TEST_CASE("Vector2D addition")
{
    Vector2D a{1, 2};
    Vector2D b{3, 4};

    Vector2D c = a + b;
    CHECK(c.x == doctest::Approx(4));
    CHECK(c.y == doctest::Approx(6));
}

TEST_CASE("Vector2D length")
{
    Vector2D v{3, 4};

    CHECK(v.length() == doctest::Approx(5));
}