#include "../third_party/doctest.h"
#include "Algebra/Transform3D.hpp"

using namespace Arns::geomath;
using namespace Arns::geomath::literals;

TEST_SUITE("Transform3D")
{
    TEST_CASE("point, vector, and inverse transformations")
    {
        const Transform3D transform(
            Vector3D(3, -2, 5), EulerAngles(Degrees(0), Degrees(0), Degrees(90)), real_t{2});

        CHECK(transform.transformPoint(Vector3D(1, 0, 0)) == Vector3D(3, 0, 5));
        CHECK(transform.transformVector(Vector3D(1, 0, 0)) == Vector3D(0, 2, 0));
        CHECK(transform.inverseTransformPoint(transform.transformPoint(Vector3D(1, 2, 3))) == Vector3D(1, 2, 3));
        CHECK(transform.inverseTransformVector(transform.transformVector(Vector3D(1, 2, 3))) == Vector3D(1, 2, 3));
    }

    TEST_CASE("setters and modifiers update cached matrices")
    {
        Transform3D transform;
        transform.translate(Vector3D(1, 2, 3)).rotate(0_deg, 0_deg, 90_deg).scale(2);
        CHECK(transform.transformPoint(Vector3D(1, 0, 0)) == Vector3D(1, 4, 3));

        transform.setPosition(Vector3D(0, 0, 0));
        transform.setRotation(EulerAngles());
        transform.setScale(real_t{2});
        CHECK(transform.transformPoint(Vector3D(1, 1, 1)) == Vector3D(2, 2, 2));

        transform.reset();
        CHECK(transform.transformPoint(Vector3D(1, 2, 3)) == Vector3D(1, 2, 3));
    }

    TEST_CASE("Matrix4x4 transforms points using homogeneous coordinates")
    {
        const Matrix4x4 matrix = Matrix4x4::createTranslation3D(4, 5, 6);
        CHECK(matrix.transformPoint(Vector3D(1, 2, 3)) == Vector3D(5, 7, 9));
        CHECK(matrix.transformVector(Vector3D(1, 2, 3)) == Vector3D(1, 2, 3));
    }

    TEST_CASE("composition and compound composition")
    {
        const Transform3D parent(Vector3D(2, 3, 4), EulerAngles(Degrees(0), Degrees(0), Degrees(90)), real_t{2});
        const Transform3D child(Vector3D(1, 0, 0), EulerAngles(Degrees(0), Degrees(0), Degrees(90)), real_t{3});
        const Transform3D composed = parent * child;

        CHECK(composed.getPosition() == Vector3D(2, 5, 4));
        CHECK(composed.transformPoint(Vector3D(1, 0, 0)) == Vector3D(-4, 5, 4));

        Transform3D compound = parent;
        compound *= child;
        CHECK(compound.transformPoint(Vector3D(1, 0, 0)) == composed.transformPoint(Vector3D(1, 0, 0)));
    }

    TEST_CASE("identity and zero scale")
    {
        const Transform3D identity = Transform3D::identity();
        CHECK(identity.transformPoint(Vector3D(2, 3, 0)) == Vector3D(2, 3, 0));
        CHECK(identity.transformVector(Vector3D(2, 3, 0)) == Vector3D(2, 3, 0));

        const Transform3D zeroScale(Vector3D(2, 3, 0), EulerAngles(), real_t{0});
        const Vector3D inverse = zeroScale.inverseTransformPoint(Vector3D(4, 6, 0));
        CHECK(inverse == Vector3D(0, 0, 0));
    }
}