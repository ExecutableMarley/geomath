#include "../third_party/doctest.h"
#include "Geometry/Vector3D.hpp"

using namespace Arns::Math;

TEST_SUITE("Vector3D")
{
    TEST_CASE("Constructors")
    {
        const Vector3D defaultVector;
        CHECK(defaultVector.x == doctest::Approx(0.0f));
        CHECK(defaultVector.y == doctest::Approx(0.0f));
        CHECK(defaultVector.z == doctest::Approx(0.0f));
        CHECK(defaultVector.isZero());

        const Vector3D value(1.0f, 2.0f, 3.0f);
        CHECK(value.x == doctest::Approx(1.0f));
        CHECK(value.y == doctest::Approx(2.0f));
        CHECK(value.z == doctest::Approx(3.0f));
    }

    TEST_CASE("arithmetic operators produce expected results")
    {
        const Vector3D a(1.0f, 2.0f, 3.0f);
        const Vector3D b(4.0f, 5.0f, 6.0f);

        CHECK((a + b) == Vector3D(5.0f, 7.0f, 9.0f));
        CHECK((a - b) == Vector3D(-3.0f, -3.0f, -3.0f));
        CHECK((a * 2.0f) == Vector3D(2.0f, 4.0f, 6.0f));
        CHECK((a / 2.0f) == Vector3D(0.5f, 1.0f, 1.5f));

        Vector3D c = a;
        c += b;
        CHECK(c == Vector3D(5.0f, 7.0f, 9.0f));

        c = a;
        c -= b;
        CHECK(c == Vector3D(-3.0f, -3.0f, -3.0f));

        c = a;
        c *= 2.0f;
        CHECK(c == Vector3D(2.0f, 4.0f, 6.0f));

        c = a;
        c /= 2.0f;
        CHECK(c == Vector3D(0.5f, 1.0f, 1.5f));

        Vector3D value(2.0f, -3.0f, 4.0f);
        const Vector3D negated = -value;
        CHECK(negated == Vector3D(-2.0f, 3.0f, -4.0f));
        CHECK(negated.lengthSquared() == doctest::Approx(value.lengthSquared()));
    }

    TEST_CASE("Geometric Properties")
    {
        const Vector3D value(1.0f, 2.0f, 2.0f);

        CHECK(value.length() == doctest::Approx(3.0f));
        CHECK(value.lengthSquared() == doctest::Approx(9.0f));
        CHECK(value.distance(Vector3D(4.0f, 6.0f, 2.0f)) == doctest::Approx(5.0f));
        CHECK(value.distanceSquared(Vector3D(4.0f, 6.0f, 2.0f)) == doctest::Approx(25.0f));
        CHECK(value.dot(Vector3D(2.0f, 0.0f, 1.0f)) == doctest::Approx(4.0f));
        CHECK(value.cross(Vector3D(2.0f, 0.0f, 1.0f)) == Vector3D(2.0f, 3.0f, -4.0f));
    }

    TEST_CASE("State Queries")
    {
        const Vector3D zero;
        const Vector3D axis(3.0f, 0.0f, 0.0f);
        const Vector3D parallel(6.0f, 0.0f, 0.0f);
        const Vector3D orthogonal(0.0f, 1.0f, 0.0f);

        CHECK(zero.isZero());
        CHECK(axis.isParallel(parallel));
        CHECK(axis.isOrthogonal(orthogonal));
        CHECK_FALSE(axis.isParallel(orthogonal));
    }

    TEST_CASE("Transform / Modification (Normalize & resize)")
    {
        Vector3D value(3.0f, 4.0f, 0.0f);
        value.normalize();
        CHECK(value == Vector3D(0.6f, 0.8f, 0.0f));
        CHECK(value.isNormalized());

        value = Vector3D(3.0f, 4.0f, 0.0f);
        value.resize(10.0f);
        CHECK(value == Vector3D(6.0f, 8.0f, 0.0f));

        Vector3D zeroVector;
        zeroVector.normalize();
        CHECK(zeroVector == Vector3D::zero());
        CHECK(zeroVector.length() == doctest::Approx(0.0f));
    }

    TEST_CASE("Transform / Modification (rotation)")
    {
        Vector3D value(0.0f, 1.0f, 0.0f);
        value.rotateAroundX(PI / 2.0f);
        CHECK(value == Vector3D(0.0f, 0.0f, 1.0f));

        value = Vector3D(0.0f, 0.0f, 1.0f);
        value.rotateAroundY(PI / 2.0f);
        CHECK(value == Vector3D(1.0f, 0.0f, 0.0f));

        value = Vector3D(1.0f, 0.0f, 0.0f);
        value.rotateAroundZ(PI / 2.0f);
        CHECK(value == Vector3D(0.0f, 1.0f, 0.0f));
    }

    TEST_CASE("Transform / Modification (Clamping)")
    {
        Vector3D value(5.0f, 2.0f, 7.0f);
        value.clamp(Vector3D(0.0f, 0.0f, 0.0f), Vector3D(4.0f, 4.0f, 6.0f));
        CHECK(value == Vector3D(4.0f, 2.0f, 6.0f));

        const Vector3D copied = value.copy();
        CHECK(copied == value);
    }

    TEST_CASE("Constants & Helpers")
    {
        CHECK(Vector3D::zero() == Vector3D(0.0f, 0.0f, 0.0f));
        CHECK(Vector3D::unitX() == Vector3D(1.0f, 0.0f, 0.0f));
        CHECK(Vector3D::unitY() == Vector3D(0.0f, 1.0f, 0.0f));
        CHECK(Vector3D::unitZ() == Vector3D(0.0f, 0.0f, 1.0f));

        CHECK(Vector3D::min(Vector3D(4.0f, 2.0f, 7.0f), Vector3D(1.0f, 5.0f, 3.0f)) == Vector3D(1.0f, 2.0f, 3.0f));
        CHECK(Vector3D::max(Vector3D(4.0f, 2.0f, 7.0f), Vector3D(1.0f, 5.0f, 3.0f)) == Vector3D(4.0f, 5.0f, 7.0f));

        CHECK(Arns::Math::dot(Vector3D(1.0f, 2.0f, 3.0f), Vector3D(4.0f, 5.0f, 6.0f)) == doctest::Approx(32.0f));
        CHECK(Arns::Math::cross(Vector3D(1.0f, 0.0f, 0.0f), Vector3D(0.0f, 1.0f, 0.0f)) == Vector3D(0.0f, 0.0f, 1.0f));
    }
}