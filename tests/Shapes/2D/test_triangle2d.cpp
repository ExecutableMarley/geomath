#include "../third_party/doctest.h"

#include "../../CommonMath.hpp"
#include "../../Shapes/2D/Triangle2D.hpp"
#include "test_shape2D_utility.hpp"

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