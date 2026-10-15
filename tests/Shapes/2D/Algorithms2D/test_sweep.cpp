#include "../../../third_party/doctest.h"

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Shapes/2D/Algorithms2D/Algorithms2D.hpp"
#include "Shapes/2D/BBox2D.hpp"
#include "Shapes/2D/Ray2D.hpp"
#include "Shapes/2D/Line2D.hpp"
#include "Shapes/2D/Triangle2D.hpp"
#include "Shapes/2D/Rectangle2D.hpp"
#include "Shapes/2D/Polygon2D.hpp"
#include "Shapes/2D/ConvexPolygon2D.hpp"
#include "Shapes/2D/Circle2D.hpp"
#include "Shapes/2D/Algorithms2D/ShapeVariantAlgorithms.hpp"

using namespace Arns::geomath;

TEST_SUITE("sweep")
{
    TEST_CASE("Circle vs Circle")
    {
        SUBCASE("Head-on Collision")
        {
            Circle2D a(Vector2D(0, 0), 1.0);
            Circle2D b(Vector2D(5, 0), 1.0);

            Vector2D velA(2.0, 0.0);
            Vector2D velB(-2.0, 0.0);

            // Distance between centers is 5. Sum of radii is 2.0
            // Initial gap = 3.0 units.
            // Relative velocity = 4.0 u/s towards each other.
            // Expected TOI = 3.0 / 4.0 = 0.75
            SweepResult2D result = sweep(a, velA, b, velB, 1.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(0.75));
            CHECK(result.normal == Vector2D(1.0, 0.0));
        }

        SUBCASE("Near Miss (Parallel Motion)")
        {
            Circle2D a(Vector2D(0, 0), 1.0);
            Circle2D b(Vector2D(5, 2.5), 1.0);

            Vector2D velA(4.0, 0.0);
            Vector2D velB(0.0, 0.0);

            // Closest approach along Y is 2.5, but combined radii = 2.0
            // The circles will bypass each other without overlapping.
            SweepResult2D result = sweep(a, velA, b, velB, 1.0);

            CHECK_FALSE(result.hit);
        }

        SUBCASE("Initial Penetration / Touching at t = 0")
        {
            Circle2D a(Vector2D(0, 0), 1.0);
            Circle2D b(Vector2D(2, 0), 1.0); // Exact touch at (1,0)

            Vector2D velA(1.0, 0.0);
            Vector2D velB(-1.0, 0.0);

            SweepResult2D result = sweep(a, velA, b, velB, 1.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(0.0));
        }

        SUBCASE("Touching while Separating")
        {
            Circle2D a(Vector2D(0, 0), 1.0);
            Circle2D b(Vector2D(2, 0), 1.0);

            SweepResult2D result = sweep(a, Vector2D(-1, 0), b, Vector2D(0, 0), 1.0);

            CHECK_FALSE(result.hit);
        }

        SUBCASE("Penetrating while Separating")
        {
            Circle2D a(Vector2D(0, 0), 1.0);
            Circle2D b(Vector2D(1.5, 0), 1.0);

            SweepResult2D result = sweep(a, Vector2D(-1, 0), b, Vector2D(0, 0), 1.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(0.0));
        }

        SUBCASE("Collision Outside maxTime Horizon")
        {
            Circle2D a(Vector2D(0, 0), 1.0);
            Circle2D b(Vector2D(10, 0), 1.0);

            Vector2D velA(2.0, 0.0);
            Vector2D velB(0.0, 0.0);

            // Gap = 8.0, Speed = 2.0 -> TOI = 4.0.
            // Limit maxTime to 1.0 second.
            SweepResult2D result = sweep(a, velA, b, velB, 1.0);

            CHECK_FALSE(result.hit);
        }
    }

    TEST_CASE("Circle vs Polygon (Rectangle)")
    {
        SUBCASE("Dynamic moving into Static")
        {
            Circle2D circle(Vector2D(0, 0), 0.5);

            // A 2x2 box centered around (4, 0) -> X spans [3, 5], Y spans [-1, 1]
            Rectangle2D box(
                Vector2D(3, -1),
                Vector2D(5, -1),
                Vector2D(5, 1),
                Vector2D(3, 1));

            Vector2D velCircle(4.0, 0.0);
            Vector2D velBox(0.0, 0.0);

            // Circle edge is at X = 0.5. Box left edge is at X = 3.0
            // Initial gap = 2.5. Velocity = 4.0.
            // Expected TOI = 2.5 / 4.0 = 0.625
            SweepResult2D result = sweep(circle, velCircle, box, velBox, 1.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(0.625));
            CHECK(result.normal == Vector2D(1.0, 0.0));
        }
    }

    TEST_CASE("BBox vs BBox")
    {
        SUBCASE("Touching while Separating")
        {
            BBox2D a(Vector2D(-1, -1), Vector2D(1, 1));
            BBox2D b(Vector2D(1, -1), Vector2D(3, 1));

            SweepResult2D result = detail::SolveBBoxToBBoxCCD(
                a, Vector2D(-1, 0), b, Vector2D(0, 0), 1.0);

            CHECK_FALSE(result.hit);
        }

        SUBCASE("Touching while Closing")
        {
            BBox2D a(Vector2D(-1, -1), Vector2D(1, 1));
            BBox2D b(Vector2D(1, -1), Vector2D(3, 1));

            SweepResult2D result = detail::SolveBBoxToBBoxCCD(
                a, Vector2D(1, 0), b, Vector2D(0, 0), 1.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(0.0));
        }

        SUBCASE("Penetrating while Separating")
        {
            BBox2D a(Vector2D(-1, -1), Vector2D(1, 1));
            BBox2D b(Vector2D(0.5, -1), Vector2D(2.5, 1));

            SweepResult2D result = detail::SolveBBoxToBBoxCCD(
                a, Vector2D(-1, 0), b, Vector2D(0, 0), 1.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(0.0));
        }
    }

    TEST_CASE("Polygon vs Polygon")
    {
        SUBCASE("Corner to Edge Collision (Triangle vs Rectangle)")
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

            // Tri tip at X=1, Wall edge at X=3. Gap = 2.0. Speed = 2.0.
            // Expected TOI = 1.0.
            SweepResult2D result = sweep(tri, velTri, wall, velWall, 2.0);

            CHECK(result.hit);
            CHECK(result.toi == doctest::Approx(1.0));
        }

        SUBCASE("Arbitrary Polygon vs Polygon")
        {
            // Concave/convex generic 5-point polygon
            Polygon2D polyA({Vector2D(-1, -1),
                             Vector2D(1, -1),
                             Vector2D(1, 1),
                             Vector2D(0, 1.5),
                             Vector2D(-1, 1)});

            Polygon2D polyB({Vector2D(4, 0),
                             Vector2D(6, -1),
                             Vector2D(6, 1)});

            Vector2D velA(3.0, 0.0);
            Vector2D velB(-1.0, 0.0);

            // Closing speed = 4.0.
            SweepResult2D result = sweep(polyA, velA, polyB, velB, 2.0);

            CHECK(result.hit);
            CHECK(result.toi > 0.0);
            CHECK(result.toi < 2.0);
        }
    }
}