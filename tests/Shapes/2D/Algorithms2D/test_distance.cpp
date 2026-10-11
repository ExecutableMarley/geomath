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
#include "Shapes/2D/Algorithms2D/ShapeVariantAlgorithms.hpp"

using namespace Arns::geomath;

// ============================================================================
// 
// ============================================================================

TEST_CASE("distance point-segment")
{
    Segment2D seg{Vector2D(1, 0), Vector2D(1, 1)};
    Vector2D closest;

    SUBCASE("point projects to segment endpoint")
    {
        real_t d1 = distance(Vector2D(0, 0), seg, nullptr, &closest);
        CHECK(doctest::Approx(d1) == real_t{1.0});
        CHECK(closest == Vector2D(1, 0));
    }

    SUBCASE("point on the segment interior")
    {
        real_t dInside = distance(Vector2D(1, 0.5), seg, nullptr, &closest);
        CHECK(doctest::Approx(dInside) == real_t{0.0});
        CHECK(closest == Vector2D(1, 0.5));
    }

    SUBCASE("point projects orthogonally onto interior")
    {
        real_t dPerp = distance(Vector2D(3, 0.5), seg, nullptr, &closest);
        CHECK(doctest::Approx(dPerp) == real_t{2.0});
        CHECK(closest == Vector2D(1, 0.5));
    }

    SUBCASE("point beyond top endpoint")
    {
        real_t dFarTop = distance(Vector2D(1, 3), seg, nullptr, &closest);
        CHECK(doctest::Approx(dFarTop) == real_t{2.0});
        CHECK(closest == Vector2D(1, 1));
    }
}

TEST_CASE("distance point-bbox")
{
    BBox2D bbox(Vector2D(0, 0), Vector2D(1, 1));
    Vector2D closest;

    SUBCASE("point outside corner (top-right)")
    {
        real_t d1 = distance(Vector2D(2, 3), bbox, nullptr, &closest);
        CHECK(doctest::Approx(d1) == std::sqrt(real_t{5.0}));
        CHECK(closest == Vector2D(1, 1));
    }

    SUBCASE("point strictly inside bbox")
    {
        real_t dInside = distance(Vector2D(0.5, 0.5), bbox, nullptr, &closest);
        CHECK(doctest::Approx(dInside) == real_t{0.0});
        CHECK(closest == Vector2D(0.5, 0.5));
    }

    SUBCASE("point outside edge (left)")
    {
        real_t dLeft = distance(Vector2D(-2, 0.5), bbox, NULL, &closest);
        CHECK(doctest::Approx(dLeft) == real_t{2.0});
        CHECK(closest == Vector2D(0, 0.5));
    }
}

TEST_CASE("distance point-polygon and point-circle")
{
    Vector2D closest;

    SUBCASE("triangle queries")
    {
        Triangle2D tri(Vector2D(0, 0), Vector2D(2, 0), Vector2D(0, 2));

        SUBCASE("point inside triangle")
        {
            real_t dInside = distance(Vector2D(0.5, 0.5), tri, nullptr, &closest);
            CHECK(doctest::Approx(dInside) == real_t{0.0});
            CHECK(closest == Vector2D(0.5, 0.5));
        }

        SUBCASE("point on hypotenuse boundary")
        {
            real_t dOnEdge = distance(Vector2D(1, 1), tri, nullptr, &closest);
            CHECK(doctest::Approx(dOnEdge) == real_t{0.0});
            CHECK(closest == Vector2D(1, 1));
        }

        SUBCASE("point outside triangle hypotenuse")
        {
            real_t dOutside = distance(Vector2D(3, 3), tri, nullptr, &closest);
            CHECK(doctest::Approx(dOutside) == std::sqrt(real_t{8.0}));
            CHECK(closest == Vector2D(1, 1)); // Closest point on hypotenuse x + y = 2
        }
    }

    SUBCASE("circle queries")
    {
        Circle2D circle(Vector2D(2, 0), real_t{1.0});

        SUBCASE("point outside circle")
        {
            real_t dc1 = distance(Vector2D(0, 0), circle, nullptr, &closest);
            CHECK(doctest::Approx(dc1) == real_t{1.0});
            CHECK(closest == Vector2D(1, 0));
        }

        SUBCASE("point at circle center")
        {
            real_t dcCenter = distance(Vector2D(2, 0), circle, nullptr, &closest);
            CHECK(doctest::Approx(dcCenter) == real_t{0.0});
            CHECK(closest == Vector2D(2, 0));
        }

        SUBCASE("point inside circle but non-center")
        {
            real_t dcInside = distance(Vector2D(2.5, 0), circle, nullptr, &closest);
            CHECK(doctest::Approx(dcInside) == real_t{0.0});
            CHECK(closest == Vector2D(2.5, 0));
        }
    }
}

TEST_CASE("distance segment-segment")
{
    Vector2D c1, c2;

    SUBCASE("parallel segments")
    {
        Segment2D a(Vector2D(0, 0), Vector2D(1, 0));
        Segment2D b(Vector2D(0, 1), Vector2D(1, 1));

        real_t d = distance(a, b, &c1, &c2);
        CHECK(doctest::Approx(d) == real_t{1.0});
        CHECK(c1 == Vector2D(0, 0));
        CHECK(c2 == Vector2D(0, 1));
    }

    SUBCASE("intersecting segments")
    {
        Segment2D a(Vector2D(0, 0), Vector2D(1, 0));
        Segment2D c(Vector2D(0.5, -1), Vector2D(0.5, 1));

        real_t d0 = distance(a, c, &c1, &c2);
        CHECK(doctest::Approx(d0) == real_t{0.0});
    }

    SUBCASE("degenerate point segment")
    {
        Segment2D point(Vector2D(0, 0), Vector2D(0, 0));
        Segment2D other(Vector2D(-1, 1), Vector2D(1, 1));

        real_t dDegenerate = distance(point, other, &c1, &c2);
        CHECK(doctest::Approx(dDegenerate) == real_t{1.0});
        CHECK(c1 == Vector2D(0, 0));
        CHECK(c2 == Vector2D(0, 1));
    }

    SUBCASE("skewed line intersection outside segment bounds")
    {
        Segment2D a(Vector2D(0, 1), Vector2D(2, 3));
        Segment2D b(Vector2D(0, 0), Vector2D(1, -1));

        real_t d = distance(a, b, &c1, &c2);

        // Expected closest point on A is (0, 1) and on B is (0, 0)
        CHECK(c1 == Vector2D(0, 1));
        CHECK(c2 == Vector2D(0, 0));
        CHECK(doctest::Approx(d) == real_t{1.0});
    }
}

TEST_CASE("distance segment-bbox")
{
    Vector2D c1, c2;

    SUBCASE("separation and partial intersection")
    {
        BBox2D bbox(Vector2D(0, 0), Vector2D(1, 1));

        Segment2D seg(Vector2D(2, 2), Vector2D(3, 2));
        real_t d = distance(seg, bbox, &c1, &c2);
        CHECK(doctest::Approx(d) == std::sqrt(real_t{2}));
        CHECK(c1 == Vector2D(2, 2));
        CHECK(c2 == Vector2D(1, 1));

        Segment2D seg2(Vector2D(0.5, 0.5), Vector2D(0.5, -1.0));
        real_t d0 = distance(seg2, bbox, &c1, &c2);
        CHECK(doctest::Approx(d0) == real_t{0});
    }

    SUBCASE("containment and boundary intersection")
    {
        const BBox2D box(Vector2D(0, 0), Vector2D(2, 2));

        SUBCASE("segment fully inside box")
        {
            const Segment2D segment(Vector2D(0.5, 0.5), Vector2D(1.5, 1.5));
            CHECK(distance(segment, box, &c1, &c2) == doctest::Approx(0));
        }

        SUBCASE("segment crosses box boundary")
        {
            const Segment2D segment(Vector2D(-1, 1), Vector2D(1, 1));
            CHECK(distance(segment, box, &c1, &c2) == doctest::Approx(0));
        }
    }
}

TEST_CASE("distance segment-polygon and segment-circle")
{
    Vector2D c1, c2;

    SUBCASE("segment-polygon")
    {
        Triangle2D tri(Vector2D(0, 0), Vector2D(2, 0), Vector2D(0, 2));

        SUBCASE("separation")
        {
            Segment2D seg(Vector2D(3, 3), Vector2D(4, 3));
            real_t d = distance(seg, tri, &c1, &c2);
            CHECK(doctest::Approx(d) == std::sqrt(real_t{8}));
            CHECK(c1 == Vector2D(3, 3));
            CHECK(c2 == Vector2D(1, 1));
        }

        SUBCASE("segment fully inside polygon")
        {
            const Segment2D segment(Vector2D(0.25, 0.25), Vector2D(0.75, 0.75));
            CHECK(distance(segment, tri, &c1, &c2) == doctest::Approx(0));
        }

        SUBCASE("segment crosses polygon boundary")
        {
            const Segment2D segment(Vector2D(-1, 0.5), Vector2D(1, 0.5));
            CHECK(distance(segment, tri, &c1, &c2) == doctest::Approx(0));
        }
    }

    SUBCASE("segment-circle")
    {
        Circle2D circle(Vector2D(0, 0), real_t{1});
        Segment2D seg2(Vector2D(2, 0), Vector2D(3, 0));
        real_t d2 = distance(seg2, circle, &c1, &c2);
        CHECK(doctest::Approx(d2) == real_t{1});
        CHECK(c2 == Vector2D(1, 0));
    }
}

// ============================================================================
// 
// ============================================================================

TEST_CASE("Polygon vs Polygon (Triangle vs Rectangle) - Corner to Edge Collision")
{
    // Right triangle pointing its tip toward the right
    Triangle2D tri(
        Vector2D(0, -1),
        Vector2D(1, 0), // Tip at (1,0)
        Vector2D(0, 1));

    // Vertical wall rectangle at X = 3
    Rectangle2D wall(
        Vector2D(3, -2),
        Vector2D(4, -2),
        Vector2D(4, 2),
        Vector2D(3, 2));

    Vector2D velTri(2.0, 0.0);
    Vector2D velWall(0.0, 0.0);

    Vector2D bestA, bestB;
    
    real_t d = distance(tri, wall, &bestA, &bestB);
    CHECK(d == doctest::Approx(2.0));
    CHECK(bestA == Vector2D(1, 0));
    CHECK(bestB == Vector2D(3, 0));
}

TEST_CASE("distance bbox-polygon separated, containment, and boundary intersection")
{
    Vector2D c1, c2;

    SUBCASE("separated shapes")
    {
        const BBox2D box(Vector2D(0, 0), Vector2D(1, 1));
        const Triangle2D triangle(Vector2D(3, 0), Vector2D(4, 0), Vector2D(3, 1));
        CHECK(distance(box, triangle, &c1, &c2) == doctest::Approx(2));
        CHECK(c1 == Vector2D(1, 0));
        CHECK(c2 == Vector2D(3, 0));
    }

    SUBCASE("box fully inside polygon")
    {
        const BBox2D box(Vector2D(0, 0), Vector2D(1, 1));
        const Triangle2D triangle(Vector2D(-2, -2), Vector2D(4, -2), Vector2D(-2, 4));
        CHECK(distance(box, triangle, &c1, &c2) == doctest::Approx(0));
    }

    SUBCASE("polygon fully inside box")
    {
        const BBox2D box(Vector2D(0, 0), Vector2D(3, 3));
        const Triangle2D triangle(Vector2D(1, 1), Vector2D(2, 1), Vector2D(1, 2));
        CHECK(distance(box, triangle, &c1, &c2) == doctest::Approx(0));
    }

    SUBCASE("boundaries intersect")
    {
        const BBox2D box(Vector2D(0, 0), Vector2D(1, 1));
        const Triangle2D triangle(Vector2D(0.5, 0.5), Vector2D(2, 0.5), Vector2D(0.5, 2));
        CHECK(distance(box, triangle, &c1, &c2) == doctest::Approx(0));
    }
}

TEST_CASE("distance polygon-circle separated, containment, and boundary intersection")
{
    Vector2D c1, c2;

    SUBCASE("separated shapes")
    {
        const Triangle2D triangle(Vector2D(0, 0), Vector2D(4, 0), Vector2D(0, 4));
        const Circle2D circle(Vector2D(5, 1), real_t{0.5});
        CHECK(distance(triangle, circle, &c1, &c2) == doctest::Approx(std::sqrt(real_t{2}) - real_t{0.5}));
        CHECK(c1 == Vector2D(4, 0));
        const real_t radiusOffset = real_t{0.5} / std::sqrt(real_t{2});
        CHECK(c2 == Vector2D(real_t{5} - radiusOffset, real_t{1} - radiusOffset));
    }

    SUBCASE("circle fully inside polygon")
    {
        const Triangle2D triangle(Vector2D(0, 0), Vector2D(4, 0), Vector2D(0, 4));
        const Circle2D circle(Vector2D(0.5, 0.5), real_t{0.25});
        CHECK(distance(triangle, circle, &c1, &c2) == doctest::Approx(0));
        CHECK(c1 == circle.m_center);
        CHECK(c2 == circle.m_center);
    }

    SUBCASE("polygon fully inside circle")
    {
        const Triangle2D triangle(Vector2D(0, 0), Vector2D(1, 0), Vector2D(0, 1));
        const Circle2D circle(Vector2D(0.25, 0.25), real_t{2});
        CHECK(distance(triangle, circle, &c1, &c2) == doctest::Approx(0));
    }

    SUBCASE("boundaries intersect")
    {
        const Triangle2D triangle(Vector2D(0, 0), Vector2D(4, 0), Vector2D(0, 4));
        const Circle2D circle(Vector2D(1, 1), real_t{1});
        CHECK(distance(triangle, circle, &c1, &c2) == doctest::Approx(0));
    }
}

TEST_CASE("distance polygon-polygon containment and boundary intersection")
{
    Vector2D c1, c2;

    SUBCASE("one polygon fully inside the other")
    {
        const Triangle2D outer(Vector2D(0, 0), Vector2D(4, 0), Vector2D(0, 4));
        const Triangle2D inner(Vector2D(0.5, 0.5), Vector2D(1, 0.5), Vector2D(0.5, 1));
        CHECK(distance(outer, inner, &c1, &c2) == doctest::Approx(0));
    }

    SUBCASE("boundaries intersect")
    {
        const Rectangle2D first(Vector2D(0, 0), real_t{2}, real_t{2});
        const Rectangle2D second(Vector2D(1, 1), real_t{2}, real_t{2});
        CHECK(distance(first, second, &c1, &c2) == doctest::Approx(0));
    }
}

TEST_CASE("distance variant dispatch handles reversed shape order")
{
    const Triangle2D triangle(Vector2D(0, 0), Vector2D(4, 0), Vector2D(0, 4));
    const Circle2D circle(Vector2D(5, 1), real_t{0.5});
    const FiniteShapeVariant2D triangleVariant{triangle};
    const FiniteShapeVariant2D circleVariant{circle};
    const real_t expectedDistance = std::sqrt(real_t{2}) - real_t{0.5};

    Vector2D circlePoint;
    Vector2D trianglePoint;
    const real_t d = distance(circleVariant, triangleVariant, &circlePoint, &trianglePoint);

    CHECK(d == doctest::Approx(expectedDistance));
    CHECK(circlePoint == Vector2D(real_t{5} - real_t{0.5} / std::sqrt(real_t{2}),
                                  real_t{1} - real_t{0.5} / std::sqrt(real_t{2})));
    CHECK(trianglePoint == Vector2D(4, 0));
    CHECK(distance(circleVariant, triangleVariant, nullptr, nullptr) == doctest::Approx(expectedDistance));
}