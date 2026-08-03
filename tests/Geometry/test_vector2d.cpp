#include "../third_party/doctest.h"
#include "Geometry/Vector2D.hpp"

using namespace Arns::Math;

TEST_SUITE("Vector2D")
{
    TEST_CASE("Constructors")
    {
        const Vector2D defaultVector;
        CHECK(defaultVector.x == doctest::Approx(0.0f));
        CHECK(defaultVector.y == doctest::Approx(0.0f));
        CHECK(defaultVector.isZero());

        const Vector2D value(3.0f, 4.0f);
        CHECK(value.x == doctest::Approx(3.0f));
        CHECK(value.y == doctest::Approx(4.0f));
    }

    TEST_CASE("arithmetic operators produce expected results")
    {
        const Vector2D a(3.0f, 4.0f);
        const Vector2D b(1.0f, 2.0f);

        CHECK((a + b) == Vector2D(4.0f, 6.0f));
        CHECK((a - b) == Vector2D(2.0f, 2.0f));
        CHECK((a * 2.0f) == Vector2D(6.0f, 8.0f));
        CHECK((a / 2.0f) == Vector2D(1.5f, 2.0f));

        Vector2D c = a;
        c += b;
        CHECK(c == Vector2D(4.0f, 6.0f));

        c = a;
        c -= b;
        CHECK(c == Vector2D(2.0f, 2.0f));

        c = a;
        c *= 2.0f;
        CHECK(c == Vector2D(6.0f, 8.0f));

        c = a;
        c /= 2.0f;
        CHECK(c == Vector2D(1.5f, 2.0f));

        c = -a;
        CHECK(c == Vector2D(-3.0f, -4.0f));
    }

    TEST_CASE("Geometric Properties")
    {
        const Vector2D value(3.0f, 4.0f);

        CHECK(value.length() == doctest::Approx(5.0f));
        CHECK(value.lengthSquared() == doctest::Approx(25.0f));
        CHECK(value.distance(Vector2D(1.0f, 1.0f)) == doctest::Approx(3.605551f));
        CHECK(value.distanceSquared(Vector2D(1.0f, 1.0f)) == doctest::Approx(13.0f));
        CHECK(value.dot(Vector2D(1.0f, 1.0f)) == doctest::Approx(7.0f));
        CHECK(value.cross(Vector2D(1.0f, 1.0f)) == doctest::Approx(-1.0f));

        const Vector2D zeroDistance = Vector2D(1.0f, 1.0f);
        CHECK(zeroDistance.distance(zeroDistance) == doctest::Approx(0.0f));
        CHECK(zeroDistance.distanceSquared(zeroDistance) == doctest::Approx(0.0f));
    }

    TEST_CASE("State Queries")
    {
        const Vector2D zero;
        const Vector2D horizontal(4.0f, 0.0f);
        const Vector2D parallel(8.0f, 0.0f);
        const Vector2D orthogonal(0.0f, 3.0f);

        CHECK(zero.isZero());
        CHECK(horizontal.isParallel(parallel));
        CHECK(horizontal.isOrthogonal(orthogonal));
        CHECK(horizontal.isPerpendicular(orthogonal));
        CHECK_FALSE(horizontal.isParallel(orthogonal));
    }

    TEST_CASE("Transform / Modification (Normalize & resize)")
    {
        Vector2D value(3.0f, 4.0f);
        value.normalize();
        CHECK(value == Vector2D(0.6f, 0.8f));
        CHECK(value.isNormalized());

        value = Vector2D(3.0f, 4.0f);
        value.resize(10.0f);
        CHECK(value == Vector2D(6.0f, 8.0f));

        Vector2D zeroVector;
        zeroVector.normalize();
        CHECK(zeroVector == Vector2D::zero());
        CHECK(zeroVector.length() == doctest::Approx(0.0f));
    }

    TEST_CASE("Transform / Modification (Clamping)")
    {
        Vector2D value(5.0f, 2.0f);
        value.clamp(Vector2D(0.0f, 0.0f), Vector2D(4.0f, 4.0f));
        CHECK(value == Vector2D(4.0f, 2.0f));

        const Vector2D perpendicular = Vector2D(3.0f, 4.0f).createPerpendicular();
        CHECK(perpendicular == Vector2D(-4.0f, 3.0f));

        const Vector2D unitPerpendicular = Vector2D(3.0f, 4.0f).createUnitPerpendicular();
        CHECK(unitPerpendicular.length() == doctest::Approx(1.0f));
        CHECK(unitPerpendicular == Vector2D(-0.8f, 0.6f));
    }

    TEST_CASE("Transform / Modification (rotation)")
    {
        Vector2D value(1.0f, 0.0f);
        value.rotate(90.0f);
        CHECK(value == Vector2D(0.0f, 1.0f));

        Vector2D point(2.0f, 0.0f);
        point.rotateAround(90.0f, Vector2D(1.0f, 0.0f));
        CHECK(point == Vector2D(1.0f, 1.0f));
    }

    TEST_CASE("Constants & Helpers")
    {
        CHECK(Vector2D::zero() == Vector2D(0.0f, 0.0f));
        CHECK(Vector2D::unitX() == Vector2D(1.0f, 0.0f));
        CHECK(Vector2D::unitY() == Vector2D(0.0f, 1.0f));

        CHECK(Vector2D::min(Vector2D(4.0f, 2.0f), Vector2D(1.0f, 5.0f)) == Vector2D(1.0f, 2.0f));
        CHECK(Vector2D::max(Vector2D(4.0f, 2.0f), Vector2D(1.0f, 5.0f)) == Vector2D(4.0f, 5.0f));

        CHECK(Arns::Math::dot(Vector2D(1.0f, 2.0f), Vector2D(3.0f, 4.0f)) == doctest::Approx(11.0f));
        CHECK(Arns::Math::cross(Vector2D(1.0f, 2.0f), Vector2D(3.0f, 4.0f)) == doctest::Approx(-2.0f));
        CHECK(isCCW(Vector2D(0.0f, 0.0f), Vector2D(1.0f, 0.0f), Vector2D(0.0f, 1.0f)));
        CHECK(isCW(Vector2D(0.0f, 0.0f), Vector2D(1.0f, 0.0f), Vector2D(0.0f, -1.0f)));
        CHECK(isColinear(Vector2D(0.0f, 0.0f), Vector2D(1.0f, 0.0f), Vector2D(2.0f, 0.0f)));
    }
}