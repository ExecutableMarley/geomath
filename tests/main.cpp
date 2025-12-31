#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "third_party/doctest.h"

#include "../Math.hpp"

using namespace Arns::Math;

TEST_CASE("ReadMe example")
{
    // 1. Geometric Intersections
    Triangle2D triangle({1, 0}, {0, 0}, {0, 1});
    Circle2D circle(Vector2D(-1, 0), 1.0f);

    SUBCASE("Triangle-Circle Intersection") {
        bool collides = intersect(triangle, circle);
        CHECK(collides); // expects true
    }

    // 2. Method Chaining 
    SUBCASE("Triangle scaling and perimeter") {
        Triangle2D scaled = triangle.copy().scale(2.0f);
        real_t p = scaled.perimeter();

        // Expected perimeter = original perimeter * 2
        real_t expected = triangle.perimeter() * 2;
        CHECK(doctest::Approx(p) == expected);
    }

    // 3. Algebraic Types (Fraction)
    SUBCASE("Fraction addition") {
        Fraction f1(1, 3);
        Fraction f2(1, 6);
        Fraction result = f1 + f2;

         // expects 1/2
        CHECK(result == Fraction(1, 2));
    }
    
    return;
}