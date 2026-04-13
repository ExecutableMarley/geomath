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

//----

TEST_CASE("distance point-segment")
{
    Segment2D seg{Vector2D(1,0), Vector2D(1,1)};
    Vector2D closest;

    real_t d1 = distance(Vector2D(0, 0), seg, &closest);
    CHECK(doctest::Approx(d1) == real_t{1.0});
    CHECK(closest == Vector2D(1, 0));

    real_t dInside = distance(Vector2D(1, 0.5), seg, &closest);
    CHECK(doctest::Approx(dInside) == real_t{0.0});
    CHECK(closest == Vector2D(1, 0.5));
}

TEST_CASE("distance point-bbox")
{
    BBox2D bbox(Vector2D(0, 0), Vector2D(1, 1));
    Vector2D closest;

    real_t d1 = distance(Vector2D(2, 3), bbox, &closest);
    CHECK(doctest::Approx(d1) == std::sqrt(real_t{5.0}));
    CHECK(closest == Vector2D(1, 1));

    real_t dInside = distance(Vector2D(0.5, 0.5), bbox, &closest);
    CHECK(doctest::Approx(dInside) == real_t{0.0});
    CHECK(closest == Vector2D(0.5, 0.5));
}

TEST_CASE("distance point-polygon and point-circle")
{
    Triangle2D tri(Vector2D(0, 0), Vector2D(2, 0), Vector2D(0, 2));
    Vector2D closest;

    real_t dInside = distance(Vector2D(1, 1), tri, &closest);
    CHECK(doctest::Approx(dInside) == real_t{0.0});
    CHECK(closest == Vector2D(1, 1));

    real_t dOutside = distance(Vector2D(3, 3), tri, &closest);
    CHECK(doctest::Approx(dOutside) == std::sqrt(real_t{8.0}));

    Circle2D circle(Vector2D(2, 0), 1.0);
    real_t dc1 = distance(Vector2D(0,0), circle, &closest);
    CHECK(doctest::Approx(dc1) == real_t{1.0});
    CHECK(closest == Vector2D(1, 0));

    real_t dcInside = distance(Vector2D(2, 0), circle, &closest);
    CHECK(doctest::Approx(dcInside) == real_t{0.0});
    CHECK(closest == Vector2D(2, 0));
}

TEST_CASE("distance segment-segment")
{
    Segment2D a(Vector2D(0, 0), Vector2D(1, 0));
    Segment2D b(Vector2D(0, 1), Vector2D(1, 1));
    Vector2D c1, c2;

    real_t d = distance(a, b, &c1, &c2);
    CHECK(doctest::Approx(d) == real_t{1.0});
    CHECK(c1 == Vector2D(0, 0));
    CHECK(c2 == Vector2D(0, 1));

    Segment2D c(Vector2D(0.5, -1), Vector2D(0.5, 1));
    real_t d0 = distance(a, c, &c1, &c2);
    CHECK(doctest::Approx(d0) == real_t{0.0});
}

TEST_CASE("distance segment-bbox")
{
    Segment2D seg(Vector2D(2, 2), Vector2D(3, 2));
    BBox2D bbox(Vector2D(0, 0), Vector2D(1, 1));
    Vector2D c1, c2;

    real_t d = distance(seg, bbox, &c1, &c2);
    CHECK(doctest::Approx(d) == std::sqrt(real_t{2.0}));
    CHECK(c1 == Vector2D(2, 2));
    CHECK(c2 == Vector2D(1, 1));

    Segment2D seg2(Vector2D(0.5, 0.5), Vector2D(0.5, -1.0));
    real_t d0 = distance(seg2, bbox, &c1, &c2);
    CHECK(doctest::Approx(d0) == 0.0f);
}

TEST_CASE("distance segment-polygon and segment-circle")
{
    Triangle2D tri(Vector2D(0, 0), Vector2D(2, 0), Vector2D(0, 2));
    Segment2D seg(Vector2D(3, 3), Vector2D(4, 3));
    Vector2D c1, c2;

    real_t d = distance(seg, tri, &c1, &c2);
    CHECK(doctest::Approx(d) == std::sqrt(real_t{10.0}));

    Circle2D circle(Vector2D(0, 0), 1.0);
    Segment2D seg2(Vector2D(2, 0), Vector2D(3, 0));
    real_t d2 = distance(seg2, circle, &c1, &c2);
    CHECK(doctest::Approx(d2) == real_t{1.0});
    CHECK(c2 == Vector2D(1, 0));
}
