#include "../third_party/doctest.h"

#include "../../CommonMath.hpp"
#include "../../Shapes/2D/Rectangle2D.hpp"
#include "test_shape2D_utility.hpp"

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