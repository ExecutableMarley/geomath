#include "../../../third_party/doctest.h"

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Shapes/2D/Algorithms2D/Algorithms2D.hpp"
#include "Shapes/2D/Ray2D.hpp"
#include "Shapes/2D/Line2D.hpp"
#include "Shapes/2D/Triangle2D.hpp"
#include "Shapes/2D/Rectangle2D.hpp"
#include "Shapes/2D/Polygon2D.hpp"
#include "Shapes/2D/ConvexPolygon2D.hpp"
#include "Shapes/2D/Circle2D.hpp"


using namespace Arns::Math;

// ============================================================================
// TRIANGLE2D CONTAINMENT TESTS
// ============================================================================

TEST_CASE("Triangle2D point containment")
{
    Triangle2D triangle(Vector2D{0, 0}, Vector2D{4, 0}, Vector2D{2, 3});

    SUBCASE("Point at centroid should be contained")
    {
        Vector2D centroid = triangle.centroid();
        CHECK(triangle.contains(centroid));
    }

    SUBCASE("Point at vertex should be contained")
    {
        CHECK(triangle.contains(Vector2D{0, 0}));
        CHECK(triangle.contains(Vector2D{4, 0}));
        CHECK(triangle.contains(Vector2D{2, 3}));
    }

    SUBCASE("Point on edge should be contained")
    {
        CHECK(triangle.contains(Vector2D{2, 0}));
        CHECK(triangle.contains(Vector2D{1, 1.5}));
    }

    SUBCASE("Point inside triangle should be contained")
    {
        CHECK(triangle.contains(Vector2D{2, 1}));
        CHECK(triangle.contains(Vector2D{1, 0.5}));
    }

    SUBCASE("Point outside triangle should not be contained")
    {
        CHECK(!triangle.contains(Vector2D{-1, 0}));
        CHECK(!triangle.contains(Vector2D{5, 0}));
        CHECK(!triangle.contains(Vector2D{2, 4}));
        CHECK(!triangle.contains(Vector2D{0, 2}));
    }
}

TEST_CASE("Triangle2D segment containment")
{
    Triangle2D triangle(Vector2D{0, 0}, Vector2D{4, 0}, Vector2D{2, 4});

    SUBCASE("Segment completely inside triangle should be contained")
    {
        Segment2D segment(Vector2D{1, 0.5}, Vector2D{3, 0.5});
        CHECK(triangle.contains(segment));
    }

    SUBCASE("Segment along edge should be contained")
    {
        Segment2D segment(Vector2D{1, 0}, Vector2D{3, 0});
        CHECK(triangle.contains(segment));
    }

    SUBCASE("Segment partially outside should not be contained")
    {
        Segment2D segment(Vector2D{1, 0}, Vector2D{5, 0});
        CHECK(!triangle.contains(segment));
    }

    SUBCASE("Segment completely outside should not be contained")
    {
        Segment2D segment(Vector2D{-2, 0}, Vector2D{-1, 0});
        CHECK(!triangle.contains(segment));
    }
}

TEST_CASE("Triangle2D shape containment")
{
    Triangle2D triangle(Vector2D{0, 0}, Vector2D{6, 0}, Vector2D{3, 6});

    SUBCASE("Smaller triangle inside should be contained")
    {
        Triangle2D smaller(Vector2D{2, 1}, Vector2D{4, 1}, Vector2D{3, 3});
        CHECK(triangle.contains(smaller));
    }

    SUBCASE("Triangle identical to original should be contained")
    {
        Triangle2D identical = triangle;
        CHECK(triangle.contains(identical));
    }

    SUBCASE("Triangle outside should not be contained")
    {
        Triangle2D outside(Vector2D{7, 0}, Vector2D{9, 0}, Vector2D{8, 2});
        CHECK(!triangle.contains(outside));
    }

    SUBCASE("Triangle partially overlapping should not be contained")
    {
        Triangle2D overlapping(Vector2D{4, 0}, Vector2D{8, 0}, Vector2D{6, 2});
        CHECK(!triangle.contains(overlapping));
    }
}

// ============================================================================
// RECTANGLE2D CONTAINMENT TESTS
// ============================================================================

TEST_CASE("Rectangle2D point containment")
{
    Rectangle2D rect(
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 3},
        Vector2D{0, 3}
    );

    SUBCASE("Point at centroid should be contained")
    {
        Vector2D centroid = rect.centroid();
        CHECK(rect.contains(centroid));
    }

    SUBCASE("Point at vertex should be contained")
    {
        CHECK(rect.contains(Vector2D{0, 0}));
        CHECK(rect.contains(Vector2D{4, 0}));
        CHECK(rect.contains(Vector2D{4, 3}));
        CHECK(rect.contains(Vector2D{0, 3}));
    }

    SUBCASE("Point on edge should be contained")
    {
        CHECK(rect.contains(Vector2D{2, 0}));    // on bottom edge
        CHECK(rect.contains(Vector2D{4, 1.5})); // on right edge
        CHECK(rect.contains(Vector2D{2, 3}));    // on top edge
        CHECK(rect.contains(Vector2D{0, 1.5})); // on left edge
    }

    SUBCASE("Point inside rectangle should be contained")
    {
        CHECK(rect.contains(Vector2D{2, 1.5}));
        CHECK(rect.contains(Vector2D{1, 1}));
        CHECK(rect.contains(Vector2D{3, 2}));
    }

    SUBCASE("Point outside rectangle should not be contained")
    {
        CHECK(!rect.contains(Vector2D{-1, 1.5}));
        CHECK(!rect.contains(Vector2D{5, 1.5}));
        CHECK(!rect.contains(Vector2D{2, -1}));
        CHECK(!rect.contains(Vector2D{2, 4}));
    }
}

TEST_CASE("Rectangle2D segment containment")
{
    Rectangle2D rect(
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 3},
        Vector2D{0, 3}
    );

    SUBCASE("Segment completely inside should be contained")
    {
        Segment2D segment(Vector2D{1, 1}, Vector2D{3, 1});
        CHECK(rect.contains(segment));
    }

    SUBCASE("Segment along edge should be contained")
    {
        Segment2D segment(Vector2D{1, 0}, Vector2D{3, 0});
        CHECK(rect.contains(segment));
    }

    SUBCASE("Segment from inside to outside should not be contained")
    {
        Segment2D segment(Vector2D{2, 1.5}, Vector2D{5, 1.5});
        CHECK(!rect.contains(segment));
    }

    SUBCASE("Segment completely outside should not be contained")
    {
        Segment2D segment(Vector2D{5, 1.5}, Vector2D{6, 1.5});
        CHECK(!rect.contains(segment));
    }
}

TEST_CASE("Rectangle2D shape containment")
{
    Rectangle2D rect(
        Vector2D{0, 0},
        Vector2D{6, 0},
        Vector2D{6, 4},
        Vector2D{0, 4}
    );

    SUBCASE("Smaller rectangle inside should be contained")
    {
        Rectangle2D smaller(
            Vector2D{1, 1},
            Vector2D{5, 1},
            Vector2D{5, 3},
            Vector2D{1, 3}
        );
        CHECK(rect.contains(smaller));
    }

    SUBCASE("Rectangle identical should be contained")
    {
        Rectangle2D identical = rect;
        CHECK(rect.contains(identical));
    }

    SUBCASE("Rectangle outside should not be contained")
    {
        Rectangle2D outside(
            Vector2D{7, 0},
            Vector2D{10, 0},
            Vector2D{10, 4},
            Vector2D{7, 4}
        );
        CHECK(!rect.contains(outside));
    }

    SUBCASE("Triangle inside rectangle should be contained")
    {
        Triangle2D triangle(Vector2D{1, 1}, Vector2D{5, 1}, Vector2D{3, 3});
        CHECK(rect.contains(triangle));
    }

    SUBCASE("Circle inside rectangle should be contained")
    {
        Circle2D circle(Vector2D{3, 2}, 1.0);
        CHECK(rect.contains(circle));
    }
}

// ============================================================================
// CIRCLE2D CONTAINMENT TESTS
// ============================================================================

TEST_CASE("Circle2D point containment")
{
    Circle2D circle(Vector2D{2, 2}, 2.0);

    SUBCASE("Point at center should be contained")
    {
        CHECK(circle.contains(Vector2D{2, 2}));
    }

    SUBCASE("Point on circumference should be contained")
    {
        CHECK(circle.contains(Vector2D{4, 2}));  // right point
        CHECK(circle.contains(Vector2D{0, 2}));  // left point
        CHECK(circle.contains(Vector2D{2, 4}));  // top point
        CHECK(circle.contains(Vector2D{2, 0}));  // bottom point
    }

    SUBCASE("Point inside circle should be contained")
    {
        CHECK(circle.contains(Vector2D{2, 2.5}));
        CHECK(circle.contains(Vector2D{3, 2}));
        CHECK(circle.contains(Vector2D{1, 1}));
    }

    SUBCASE("Point outside circle should not be contained")
    {
        CHECK(!circle.contains(Vector2D{5, 2}));
        CHECK(!circle.contains(Vector2D{2, 5}));
        CHECK(!circle.contains(Vector2D{5, 5}));
    }
}

TEST_CASE("Circle2D segment containment")
{
    Circle2D circle(Vector2D{0, 0}, 3.0);

    SUBCASE("Segment inside circle should be contained")
    {
        Segment2D segment(Vector2D{-1, 0}, Vector2D{1, 0});
        CHECK(circle.contains(segment));
    }

    SUBCASE("Segment with endpoints on circumference should be contained")
    {
        Segment2D segment(Vector2D{3, 0}, Vector2D{0, 3});
        CHECK(circle.contains(segment));
    }

    SUBCASE("Segment extending outside should not be contained")
    {
        Segment2D segment(Vector2D{-1, 0}, Vector2D{4, 0});
        CHECK(!circle.contains(segment));
    }

    SUBCASE("Segment completely outside should not be contained")
    {
        Segment2D segment(Vector2D{4, 0}, Vector2D{5, 0});
        CHECK(!circle.contains(segment));
    }
}

TEST_CASE("Circle2D shape containment")
{
    Circle2D circle(Vector2D{0, 0}, 5.0);

    SUBCASE("Smaller circle inside should be contained")
    {
        Circle2D smaller(Vector2D{0, 0}, 2.0);
        CHECK(circle.contains(smaller));
    }

    SUBCASE("Concentric smaller circle should be contained")
    {
        Circle2D concentric(Vector2D{0, 0}, 3.0);
        CHECK(circle.contains(concentric));
    }

    SUBCASE("Identical circle should be contained")
    {
        Circle2D identical(Vector2D{0, 0}, 5.0);
        CHECK(circle.contains(identical));
    }

    SUBCASE("Circle outside should not be contained")
    {
        Circle2D outside(Vector2D{7, 0}, 2.0);
        CHECK(!circle.contains(outside));
    }

    SUBCASE("Overlapping circle should not be contained")
    {
        Circle2D overlapping(Vector2D{4, 0}, 2.0);
        CHECK(!circle.contains(overlapping));
    }

    SUBCASE("Triangle inside circle should be contained")
    {
        Triangle2D triangle(Vector2D{-1, -1}, Vector2D{1, -1}, Vector2D{0, 1});
        CHECK(circle.contains(triangle));
    }

    SUBCASE("Rectangle inside circle should be contained")
    {
        Rectangle2D rect(
            Vector2D{-2, -2},
            Vector2D{2, -2},
            Vector2D{2, 2},
            Vector2D{-2, 2}
        );
        CHECK(circle.contains(rect));
    }
}

// ============================================================================
// POLYGON2D CONTAINMENT TESTS
// ============================================================================

TEST_CASE("Polygon2D point containment")
{
    // Pentagon (regular-ish)
    std::vector<Vector2D> vertices = {
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{5, 2},
        Vector2D{2.5, 4},
        Vector2D{-1, 2}
    };
    Polygon2D polygon(vertices);

    SUBCASE("Point at centroid should be contained")
    {
        Vector2D centroid = polygon.centroid();
        CHECK(polygon.contains(centroid));
    }

    SUBCASE("Point at vertex should be contained")
    {
        for (const auto& vertex : vertices)
        {
            CHECK(polygon.contains(vertex));
        }
    }

    SUBCASE("Point inside polygon should be contained")
    {
        CHECK(polygon.contains(Vector2D{2, 1.5}));
    }

    SUBCASE("Point outside polygon should not be contained")
    {
        CHECK(!polygon.contains(Vector2D{-2, 0}));
        CHECK(!polygon.contains(Vector2D{6, 2}));
        CHECK(!polygon.contains(Vector2D{2, 5}));
    }
}

TEST_CASE("Polygon2D segment containment")
{
    std::vector<Vector2D> vertices = {
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 3},
        Vector2D{0, 3}
    };
    Polygon2D polygon(vertices);

    SUBCASE("Segment inside polygon should be contained")
    {
        Segment2D segment(Vector2D{1, 1}, Vector2D{3, 1});
        CHECK(polygon.contains(segment));
    }

    SUBCASE("Segment partially outside should not be contained")
    {
        Segment2D segment(Vector2D{2, 1}, Vector2D{5, 1});
        CHECK(!polygon.contains(segment));
    }
}

TEST_CASE("Polygon2D shape containment")
{
    std::vector<Vector2D> vertices = {
        Vector2D{0, 0},
        Vector2D{6, 0},
        Vector2D{6, 4},
        Vector2D{0, 4}
    };
    Polygon2D polygon(vertices);

    SUBCASE("Triangle inside polygon should be contained")
    {
        Triangle2D triangle(Vector2D{1, 1}, Vector2D{5, 1}, Vector2D{3, 3});
        CHECK(polygon.contains(triangle));
    }

    SUBCASE("Triangle outside polygon should not be contained")
    {
        Triangle2D triangle(Vector2D{7, 0}, Vector2D{9, 0}, Vector2D{8, 2});
        CHECK(!polygon.contains(triangle));
    }
}

// ============================================================================
// CONVEX_POLYGON2D CONTAINMENT TESTS
// ============================================================================

TEST_CASE("ConvexPolygon2D point containment")
{
    // Square (convex)
    std::vector<Vector2D> vertices = {
        Vector2D{0, 0},
        Vector2D{4, 0},
        Vector2D{4, 4},
        Vector2D{0, 4}
    };
    ConvexPolygon2D convex(vertices);

    SUBCASE("Point at centroid should be contained")
    {
        Vector2D centroid = convex.centroid();
        CHECK(convex.contains(centroid));
    }

    SUBCASE("Point at vertex should be contained")
    {
        for (const auto& vertex : vertices)
        {
            CHECK(convex.contains(vertex));
        }
    }

    SUBCASE("Point on edge should be contained")
    {
        CHECK(convex.contains(Vector2D{2, 0}));
        CHECK(convex.contains(Vector2D{4, 2}));
    }

    SUBCASE("Point inside convex polygon should be contained")
    {
        CHECK(convex.contains(Vector2D{2, 2}));
        CHECK(convex.contains(Vector2D{1, 3}));
    }

    SUBCASE("Point outside convex polygon should not be contained")
    {
        CHECK(!convex.contains(Vector2D{-1, 2}));
        CHECK(!convex.contains(Vector2D{5, 2}));
    }
}

TEST_CASE("ConvexPolygon2D segment containment")
{
    std::vector<Vector2D> vertices = {
        Vector2D{0, 0},
        Vector2D{5, 0},
        Vector2D{5, 5},
        Vector2D{0, 5}
    };
    ConvexPolygon2D convex(vertices);

    SUBCASE("Segment inside convex polygon should be contained")
    {
        Segment2D segment(Vector2D{1, 1}, Vector2D{4, 1});
        CHECK(convex.contains(segment));
    }

    SUBCASE("Segment on edge should be contained")
    {
        Segment2D segment(Vector2D{1, 0}, Vector2D{4, 0});
        CHECK(convex.contains(segment));
    }

    SUBCASE("Segment extending outside should not be contained")
    {
        Segment2D segment(Vector2D{2, 1}, Vector2D{6, 1});
        CHECK(!convex.contains(segment));
    }
}

TEST_CASE("ConvexPolygon2D shape containment")
{
    std::vector<Vector2D> vertices = {
        Vector2D{0, 0},
        Vector2D{8, 0},
        Vector2D{8, 6},
        Vector2D{0, 6}
    };
    ConvexPolygon2D convex(vertices);

    SUBCASE("Smaller convex polygon inside should be contained")
    {
        std::vector<Vector2D> smallerVerts = {
            Vector2D{2, 1},
            Vector2D{6, 1},
            Vector2D{6, 5},
            Vector2D{2, 5}
        };
        ConvexPolygon2D smaller(smallerVerts);
        CHECK(convex.contains(smaller));
    }

    SUBCASE("Triangle inside should be contained")
    {
        Triangle2D triangle(Vector2D{1, 1}, Vector2D{7, 1}, Vector2D{4, 5});
        CHECK(convex.contains(triangle));
    }

    SUBCASE("Circle inside should be contained")
    {
        Circle2D circle(Vector2D{4, 3}, 1.5);
        CHECK(convex.contains(circle));
    }

    SUBCASE("Shape extending outside should not be contained")
    {
        Triangle2D triangle(Vector2D{7, 1}, Vector2D{10, 1}, Vector2D{8.5, 3});
        CHECK(!convex.contains(triangle));
    }
}