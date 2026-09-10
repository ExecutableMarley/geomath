#include "../third_party/doctest.h"
#include "Geometry/Vector2D.hpp"
#include <sstream>

using namespace Arns::geomath;

TEST_SUITE("Vector2D")
{
    TEST_CASE("Constructors")
    {
        const Vector2D defaultVector;
        CHECK(defaultVector.x == doctest::Approx(real_t{0}));
        CHECK(defaultVector.y == doctest::Approx(real_t{0}));
        CHECK(defaultVector.isZero());

        const Vector2D value(3, 4);
        CHECK(value.x == doctest::Approx(real_t{3}));
        CHECK(value.y == doctest::Approx(real_t{4}));
    }

    TEST_CASE("arithmetic operators produce expected results")
    {
        const Vector2D a(3, 4);
        const Vector2D b(1, 2);

        CHECK((a + b) == Vector2D(4, 6));
        CHECK((a - b) == Vector2D(2, 2));
        CHECK((a * real_t{2}) == Vector2D(6, 8));
        CHECK((a / real_t{2}) == Vector2D(real_t{1.5}, 2));

        Vector2D c = a;
        c += b;
        CHECK(c == Vector2D(4, 6));

        c = a;
        c -= b;
        CHECK(c == Vector2D(2, 2));

        c = a;
        c *= real_t{2};
        CHECK(c == Vector2D(6, 8));

        c = a;
        c /= real_t{2};
        CHECK(c == Vector2D(real_t{1.5}, 2));

        c = -a;
        CHECK(c == Vector2D(-3, -4));
    }

    TEST_CASE("Geometric Properties")
    {
        const Vector2D value(3, 4);

        CHECK(value.length() == doctest::Approx(real_t{5}));
        CHECK(value.lengthSquared() == doctest::Approx(real_t{25}));
        CHECK(value.distance(Vector2D(1, 1)) == doctest::Approx(real_t{3.605551}));
        CHECK(value.distanceSquared(Vector2D(1, 1)) == doctest::Approx(real_t{13}));
        CHECK(value.dot(Vector2D(1, 1)) == doctest::Approx(real_t{7}));
        CHECK(value.cross(Vector2D(1, 1)) == doctest::Approx(real_t{-1}));

        const Vector2D zeroDistance = Vector2D(1, 1);
        CHECK(zeroDistance.distance(zeroDistance) == doctest::Approx(real_t{0}));
        CHECK(zeroDistance.distanceSquared(zeroDistance) == doctest::Approx(real_t{0}));
    }

    TEST_CASE("State Queries")
    {
        const Vector2D zero;
        const Vector2D horizontal(4, 0);
        const Vector2D parallel(8, 0);
        const Vector2D orthogonal(0, 3);

        CHECK(zero.isZero());
        CHECK(horizontal.isParallel(parallel));
        CHECK(horizontal.isOrthogonal(orthogonal));
        CHECK(horizontal.isPerpendicular(orthogonal));
        CHECK_FALSE(horizontal.isParallel(orthogonal));
    }

    TEST_CASE("Transform / Modification (Normalize & resize)")
    {
        Vector2D value(3, 4);
        value.normalize();
        CHECK(value == Vector2D(real_t{0.6}, real_t{0.8}));
        CHECK(value.isNormalized());

        value = Vector2D(3, 4);
        value.resize(real_t{10});
        CHECK(value == Vector2D(6, 8));

        Vector2D zeroVector;
        zeroVector.normalize();
        CHECK(zeroVector == Vector2D::zero());
        CHECK(zeroVector.length() == doctest::Approx(real_t{0}));
    }

    TEST_CASE("Transform / Modification (Clamping)")
    {
        Vector2D value(5, 2);
        value.clamp(Vector2D(0, 0), Vector2D(4, 4));
        CHECK(value == Vector2D(4, 2));

        const Vector2D perpendicular = Vector2D(3, 4).createPerpendicular();
        CHECK(perpendicular == Vector2D(-4, 3));

        const Vector2D unitPerpendicular = Vector2D(3, 4).createUnitPerpendicular();
        CHECK(unitPerpendicular.length() == doctest::Approx(real_t{1}));
        CHECK(unitPerpendicular == Vector2D(real_t{-0.8}, real_t{0.6}));
    }

    TEST_CASE("Transform / Modification (rotation)")
    {
        Vector2D value(1, 0);
        value.rotate(real_t{90});
        CHECK(value == Vector2D(0, 1));

        Vector2D point(2, 0);
        point.rotateAround(real_t{90}, Vector2D(1, 0));
        CHECK(point == Vector2D(1, 1));
    }

    TEST_CASE("Constants & Helpers")
    {
        CHECK(Vector2D::zero() == Vector2D(0, 0));
        CHECK(Vector2D::unitX() == Vector2D(1, 0));
        CHECK(Vector2D::unitY() == Vector2D(0, 1));

        CHECK(Vector2D::min(Vector2D(4, 2), Vector2D(1, 5)) == Vector2D(1, 2));
        CHECK(Vector2D::max(Vector2D(4, 2), Vector2D(1, 5)) == Vector2D(4, 5));

        CHECK(geomath::dot(Vector2D(1, 2), Vector2D(3, 4)) == doctest::Approx(real_t{11}));
        CHECK(geomath::cross(Vector2D(1, 2), Vector2D(3, 4)) == doctest::Approx(real_t{-2}));
        CHECK(isCCW(Vector2D(0, 0), Vector2D(1, 0), Vector2D(0, 1)));
        CHECK(isCW(Vector2D(0, 0), Vector2D(1, 0), Vector2D(0, -1)));
        CHECK(isColinear(Vector2D(0, 0), Vector2D(1, 0), Vector2D(2, 0)));
    }

    TEST_CASE("Vector2D stream and format output")
    {
        Vector2D vector(1, 2);
        SUBCASE("Stream output is formatted correctly")
        {
            std::ostringstream oss;
            oss << vector;

            CHECK(oss.str() == "[1, 2]");
        }

        SUBCASE("std::format output is formatted correctly")
        {
            std::string formatted = std::format("{}", vector);

            CHECK(formatted == "[1, 2]");
        }

        SUBCASE("std::format with precision outputs correctly")
        {
            std::string formatted = std::format("{:.2f}", vector);

            CHECK(formatted == "[1.00, 2.00]");
        }

        SUBCASE("to_string outputs correctly")
        {
            std::string str = to_string(vector);

            CHECK(str == "[1, 2]");
        }
    }
}