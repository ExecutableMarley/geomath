#include "../third_party/doctest.h"
#include "Algebra/Transform2D.hpp"

using namespace Arns::geomath;
using namespace Arns::geomath::literals;

TEST_SUITE("Transform2D")
{
    TEST_CASE("point, vector, and inverse transformations")
    {
        const Transform2D transform(Vector2D(3, -2), Radians(PI / real_t{2}), real_t{2});

        CHECK(transform.transformPoint(Vector2D(1, 0)) == Vector2D(3, 0));
        CHECK(transform.transformVector(Vector2D(1, 0)) == Vector2D(0, 2));
        CHECK(transform.inverseTransformPoint(transform.transformPoint(Vector2D(1, 2))) == Vector2D(1, 2));
        CHECK(transform.inverseTransformVector(transform.transformVector(Vector2D(1, 2))) == Vector2D(1, 2));
    }

    TEST_CASE("setters and modifiers update cached matrices")
    {
        Transform2D transform;
        transform.translate(Vector2D(1, 2)).rotate(Radians(PI / 2)).scale(real_t{2});
        CHECK(transform.transformPoint(Vector2D(1, 0)) == Vector2D(1, 4));

        transform.setPosition(Vector2D(0, 0));
        transform.setRotation(0_rad);
        transform.setScale(real_t{3});
        CHECK(transform.transformPoint(Vector2D(1, 2)) == Vector2D(3, 6));

        transform.reset();
        CHECK(transform.transformPoint(Vector2D(1, 2)) == Vector2D(1, 2));
        CHECK(transform.getPosition() == Vector2D::zero());
        CHECK(transform.getRotation() == Radians(0));
        CHECK(transform.getScale() == real_t{1});
    }

    TEST_CASE("composition and compound composition")
    {
        const Transform2D parent(Vector2D(2, 3), Radians::halfPi(), real_t{2});
        const Transform2D child(Vector2D(1, 0), Radians::halfPi(), real_t{3});
        const Transform2D composed = parent * child;

        CHECK(composed.getPosition() == Vector2D(2, 5));
        CHECK(composed.transformPoint(Vector2D(1, 0)) == Vector2D(-4, 5));

        Transform2D compound = parent;
        compound *= child;
        CHECK(compound.transformPoint(Vector2D(1, 0)) == composed.transformPoint(Vector2D(1, 0)));
    }

    TEST_CASE("identity and zero scale")
    {
        const Transform2D identity = Transform2D::identity();
        CHECK(identity.transformPoint(Vector2D(2, 3)) == Vector2D(2, 3));
        CHECK(identity.transformVector(Vector2D(2, 3)) == Vector2D(2, 3));

        const Transform2D zeroScale(Vector2D(2, 3), Radians(0), real_t{0});
        const Vector2D inverse = zeroScale.inverseTransformPoint(Vector2D(4, 6));
        CHECK(inverse == Vector2D(0, 0));
    }
}