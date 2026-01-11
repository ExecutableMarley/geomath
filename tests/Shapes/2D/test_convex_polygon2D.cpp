#include "../third_party/doctest.h"

#include "../../CommonMath.hpp"
#include "../../Shapes/2D/ConvexPolygon2D.hpp"
#include "test_shape2D_utility.hpp"

TEST_CASE("ConvexPolygon2D polygonal interface")
{
    ConvexPolygon2D poly({
        Vector2D{0, 0},
        Vector2D{2, 0},
        Vector2D{3, 1},
        Vector2D{1, 2}
    });

    SUBCASE("Shape type is correct")
    {
        CHECK(poly.type() == SHAPE2D_CONVEX_POLYGON);
    }

    check_is_polygonal<ConvexPolygon2D>(poly, 4);
}