#include "../third_party/doctest.h"

#include "../../CommonMath.hpp"
#include "../../Shapes/2D/Circle2D.hpp"
#include "test_shape2D_utility.hpp"

TEST_CASE("Circle2D is not polygonal")
{
    Circle2D circle(Vector2D{0, 0}, real_t{1.0});

    SUBCASE("Shape type is correct")
    {
        CHECK(circle.type() == SHAPE2D_CIRCLE);
    }

    check_is_not_polygonal(circle);
}