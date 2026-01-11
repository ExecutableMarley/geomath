#include "../third_party/doctest.h"

#include "../../CommonMath.hpp"
#include "../../Shapes/2D/Polygon2D.hpp"
#include "test_shape2D_utility.hpp"

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