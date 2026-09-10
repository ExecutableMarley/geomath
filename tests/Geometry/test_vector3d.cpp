#include "../third_party/doctest.h"
#include "Geometry/Vector3D.hpp"
#include <sstream>

using namespace Arns::geomath;

TEST_SUITE("Vector3D")
{
    TEST_CASE("Constructors")
    {
        const Vector3D defaultVector;
        CHECK(defaultVector.x == doctest::Approx(real_t{0}));
        CHECK(defaultVector.y == doctest::Approx(real_t{0}));
        CHECK(defaultVector.z == doctest::Approx(real_t{0}));
        CHECK(defaultVector.isZero());

        const Vector3D value(1, 2, 3);
        CHECK(value.x == doctest::Approx(real_t{1}));
        CHECK(value.y == doctest::Approx(real_t{2}));
        CHECK(value.z == doctest::Approx(real_t{3}));
    }

    TEST_CASE("arithmetic operators produce expected results")
    {
        const Vector3D a(1, 2, 3);
        const Vector3D b(4, 5, 6);

        CHECK((a + b) == Vector3D(5, 7, 9));
        CHECK((a - b) == Vector3D(-3, -3, -3));
        CHECK((a * real_t{2}) == Vector3D(2, 4, 6));
        CHECK((a / real_t{2}) == Vector3D(real_t{0.5}, 1, real_t{1.5}));

        Vector3D c = a;
        c += b;
        CHECK(c == Vector3D(5, 7, 9));

        c = a;
        c -= b;
        CHECK(c == Vector3D(-3, -3, -3));

        c = a;
        c *= real_t{2};
        CHECK(c == Vector3D(2, 4, 6));

        c = a;
        c /= real_t{2};
        CHECK(c == Vector3D(real_t{0.5}, 1, real_t{1.5}));

        Vector3D value(2, -3, 4);
        const Vector3D negated = -value;
        CHECK(negated == Vector3D(-2, 3, -4));
        CHECK(negated.lengthSquared() == doctest::Approx(value.lengthSquared()));
    }

    TEST_CASE("Geometric Properties")
    {
        const Vector3D value(1, 2, 2);

        CHECK(value.length() == doctest::Approx(real_t{3}));
        CHECK(value.lengthSquared() == doctest::Approx(real_t{9}));
        CHECK(value.distance(Vector3D(4, 6, 2)) == doctest::Approx(real_t{5}));
        CHECK(value.distanceSquared(Vector3D(4, 6, 2)) == doctest::Approx(real_t{25}));
        CHECK(value.dot(Vector3D(2, 0, 1)) == doctest::Approx(real_t{4}));
        CHECK(value.cross(Vector3D(2, 0, 1)) == Vector3D(2, 3, -4));
    }

    TEST_CASE("State Queries")
    {
        const Vector3D zero;
        const Vector3D axis(3, 0, 0);
        const Vector3D parallel(6, 0, 0);
        const Vector3D orthogonal(0, 1, 0);

        CHECK(zero.isZero());
        CHECK(axis.isParallel(parallel));
        CHECK(axis.isOrthogonal(orthogonal));
        CHECK_FALSE(axis.isParallel(orthogonal));
    }

    TEST_CASE("Transform / Modification (Normalize & resize)")
    {
        Vector3D value(3, 4, 0);
        value.normalize();
        CHECK(value == Vector3D(real_t{0.6}, real_t{0.8}, real_t{0.0}));
        CHECK(value.isNormalized());

        value = Vector3D(3, 4, 0);
        value.resize(real_t{10});
        CHECK(value == Vector3D(6, 8, 0));

        Vector3D zeroVector;
        zeroVector.normalize();
        CHECK(zeroVector == Vector3D::zero());
        CHECK(zeroVector.length() == doctest::Approx(real_t{0}));
    }

    TEST_CASE("Transform / Modification (rotation)")
    {
        Vector3D value(0, 1, 0);
        value.rotateAroundX(PI / real_t{2});
        CHECK(value == Vector3D(0, 0, 1));

        value = Vector3D(0, 0, 1);
        value.rotateAroundY(PI / real_t{2});
        CHECK(value == Vector3D(1, 0, 0));

        value = Vector3D(1, 0, 0);
        value.rotateAroundZ(PI / real_t{2});
        CHECK(value == Vector3D(0, 1, 0));
    }

    TEST_CASE("Transform / Modification (Clamping)")
    {
        Vector3D value(5, 2, 7);
        value.clamp(Vector3D(0, 0, 0), Vector3D(4, 4, 6));
        CHECK(value == Vector3D(4, 2, 6));

        const Vector3D copied = value.copy();
        CHECK(copied == value);
    }

    TEST_CASE("Constants & Helpers")
    {
        CHECK(Vector3D::zero() == Vector3D(0, 0, 0));
        CHECK(Vector3D::unitX() == Vector3D(1, 0, 0));
        CHECK(Vector3D::unitY() == Vector3D(0, 1, 0));
        CHECK(Vector3D::unitZ() == Vector3D(0, 0, 1));

        CHECK(Vector3D::min(Vector3D(4, 2, 7), Vector3D(1, 5, 3)) == Vector3D(1, 2, 3));
        CHECK(Vector3D::max(Vector3D(4, 2, 7), Vector3D(1, 5, 3)) == Vector3D(4, 5, 7));

        CHECK(geomath::dot(Vector3D(1, 2, 3), Vector3D(4, 5, 6)) == doctest::Approx(real_t{32}));
        CHECK(geomath::cross(Vector3D(1, 0, 0), Vector3D(0, 1, 0)) == Vector3D(0, 0, 1));
    }

    TEST_CASE("Vector3D stream and format output")
    {
        Vector3D vector(1, 2, 3);
        SUBCASE("Stream output is formatted correctly")
        {
            std::ostringstream oss;
            oss << vector;

            CHECK(oss.str() == "[1, 2, 3]");
        }

        SUBCASE("std::format output is formatted correctly")
        {
            std::string formatted = std::format("{}", vector);

            CHECK(formatted == "[1, 2, 3]");
        }

        SUBCASE("std::format with precision outputs correctly")
        {
            std::string formatted = std::format("{:.2f}", vector);

            CHECK(formatted == "[1.00, 2.00, 3.00]");
        }

        SUBCASE("to_string outputs correctly")
        {
            std::string str = to_string(vector);

            CHECK(str == "[1, 2, 3]");
        }
    }
}