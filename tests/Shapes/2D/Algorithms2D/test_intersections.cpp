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
#include "Shapes/2D/Variants/ShapeVariants.hpp"
#include "Shapes/2D/Algorithms2D/ShapeVariantAlgorithms.hpp"

using namespace Arns::geomath;


struct ExpectedHit
{
    bool hit;

    // Only checked if hit == true
    Vector2D position;
    real_t t;

    // for later
    Vector2D normal; //
};

template <typename RayOrSegment, typename Shape>
void checkHit(const RayOrSegment& r, const Shape& s, const ExpectedHit& expected,
    real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max())
{
    HitInfo2D info;
    bool result = intersect(r, s, t_min, t_max, &info);

    CHECK(result == expected.hit);

    if (!expected.hit)
        return;

    //Expected t
    CHECK(doctest::Approx(info.t) == expected.t);
    //Expected position
    CHECK(info.intersectionPoint == expected.position);
    //Ray equation
    CHECK(r.origin() + r.direction() * info.t == info.intersectionPoint);
}

TEST_CASE("ray hits primitive")
{
    Ray2D ray({0,0}, {1,0});

    SUBCASE("Axis-Aligned Bounding Box")
    {
        BBox2D box({7, -1}, {9, 1});
        checkHit(ray, box, {true, {7, 0}, 7});
    }

    SUBCASE("Triangle")
    {
        Triangle2D tri({15, -1}, {18, 1}, {18, -1});
        checkHit(ray, tri, {true, {real_t{16.5}, 0}, real_t{16.5}});
    }

    SUBCASE("Rectangle")
    {
        auto rect = Rectangle2D::fromMinMax({10, -1}, {14, 1});
        checkHit(ray, rect, {true, {10, 0}, 10});
    }

    SUBCASE("Circle")
    {
        Circle2D circle({5, 0}, 1);
        checkHit(ray, circle, {true, {4, 0}, 4});
    }
}

/*
TEST_CASE("ray misses all primitives") 
{
    std::vector<std::unique_ptr<IBaseShape2D>> primitives;
    //primitives.push_back(std::make_unique<BBox2D>(BBox2D({1, -1}, {3, 1})));
    primitives.push_back(std::make_unique<Triangle2D>(Triangle2D({1, 1}, {1,-1}, {3, 1})));
    primitives.push_back(std::make_unique<Rectangle2D>(Rectangle2D({1, 1}, {1,-1}, {3, -1}, {3, 1})));
    primitives.push_back(std::make_unique<Circle2D>(Circle2D({1,0}, 2)));

    SUBCASE("Ray passes beside")
    {
        Ray2D ray({0,5}, {1,0});
        for(const auto& p : primitives)
        {
            CAPTURE(p->type());
            CHECK_FALSE(intersect(ray, *p));
        }
    }

    SUBCASE("Ray starts past shape")
    {
        Ray2D ray({5,0}, {1,0});
        for(const auto& p : primitives)
        {
            CAPTURE(p->type());
            CHECK_FALSE(intersect(ray, *p));
        }
    }
}
*/

TEST_CASE("ray starts inside primitive")
{
    Ray2D ray({0,0}, {1,0});

    SUBCASE("BBox containing ray origin")
    {
        BBox2D box({-2, -2}, {5, 2});
        checkHit(ray, box, {true, {5, 0}, 5});
    }

    SUBCASE("Triangle containing ray origin")
    {
        auto tri = Triangle2D({-2, -2}, {2, -2}, {2, 3});
        checkHit(ray, tri, {true, {2, 0}, 2});
    }

    SUBCASE("Rectangle containing ray origin")
    {
        auto rect = Rectangle2D::fromMinMax({-1, -1}, {4, 1});
        checkHit(ray, rect, {true, {4, 0}, 4});
    }

    SUBCASE("Circle containing ray origin")
    {
        Circle2D circle({0, 0}, 3);
        checkHit(ray, circle, {true, {3, 0}, 3});
    }
}

TEST_CASE("ray tangent to primitive")
{
    Ray2D ray({0,1}, {1,0});
    
    SUBCASE("BBox tangent to ray")
    {
        auto rect = BBox2D({8, 0}, {12, 1});
        checkHit(ray, rect, {true, {8, 1}, 8});
    }

    SUBCASE("Triangle tangent to ray")
    {
        Triangle2D tri({20, 1}, {23, 3}, {23, 1});
        checkHit(ray, tri, {true, {20, 1}, 20});
    }

    SUBCASE("Rectangle tangent to ray")
    {
        auto rect = Rectangle2D::fromMinMax({8, 0}, {12, 1});
        checkHit(ray, rect, {true, {8, 1}, 8});
    }

    SUBCASE("Circle tangent to ray")
    {
        Circle2D circle({5, 0}, 1);
        checkHit(ray, circle, {true, {5, 1}, 5});
    }
}

TEST_CASE("ray picks closest hit")
{
    Ray2D ray({0,0}, {1,0});

    BBox2D box({6,-1}, {8,1});
    checkHit(ray, box, {
        true, 
        {6, 0}, 
        6
    });

    Triangle2D tri({3, 1}, {3, -1}, {6, -1});
    checkHit(ray, tri, {
    true, 
    {3, 0}, 
    3
    });

    Rectangle2D rect = Rectangle2D::fromMinMax({5, -1}, {10, 1}); 
    checkHit(ray, rect, {
        true, 
        {5, 0}, 
        5
    });

    Circle2D circ({11, 0}, real_t{1});
    checkHit(ray, circ, {
        true, 
        {10, 0}, 
        10
    });
}

// ============================================================================
// Shape-Shape Intersections
// ============================================================================

TEST_CASE("BBox2D intersection")
{
    BBox2D b1(Vector2D(0, 0), Vector2D(2, 2));
    BBox2D b2(Vector2D(1, 1), Vector2D(3, 3));
    CHECK(intersect(b1, b2) == true);

    BBox2D b3(Vector2D(3, 3), Vector2D(4, 4));
    CHECK(intersect(b1, b3) == false);
}

TEST_CASE("BBox2D touching edge")
{
    BBox2D b1(Vector2D(0, 0), Vector2D(2, 2));
    BBox2D b2(Vector2D(2, 0), Vector2D(4, 2)); // touches at x = 2

    CHECK(intersect(b1, b2) == true);
}

TEST_CASE("BBox2D and IPolygonalShape2D intersection")
{
    BBox2D bbox(Vector2D(0, 0), Vector2D(3, 3));
    Triangle2D triangle(Vector2D(1, 1), Vector2D(2, 1), Vector2D(real_t{1.5}, 2));
    CHECK(intersect(bbox, triangle) == true);

    Triangle2D triangle2(Vector2D(4, 4), Vector2D(5, 4), Vector2D(real_t{4.5}, 5));
    CHECK(intersect(bbox, triangle2) == false);
}

TEST_CASE("BBox2D fully contains triangle")
{
    BBox2D bbox(Vector2D(0, 0), Vector2D(5, 5));
    Triangle2D tri(Vector2D(1, 1), Vector2D(2, 1), Vector2D(real_t{1.5}, 2));

    CHECK(intersect(bbox, tri) == true);
}

TEST_CASE("Triangle fully contains bbox corner")
{
    Triangle2D tri(Vector2D(0, 0), Vector2D(5, 0), Vector2D(0, 5));
    BBox2D bbox(Vector2D(1, 1), Vector2D(2, 2));

    CHECK(intersect(bbox, tri) == true);
}

TEST_CASE("BBox2D and Circle2D intersection")
{
    BBox2D bbox(Vector2D(0, 0), Vector2D(2, 2));
    Circle2D circle(Vector2D(1, 1), real_t{0.5});
    CHECK(intersect(bbox, circle) == true);

    Circle2D circle2(Vector2D(5, 5), 1);
    CHECK(intersect(bbox, circle2) == false);
}

TEST_CASE("Circle touching bbox edge")
{
    BBox2D bbox(Vector2D(0, 0), Vector2D(2, 2));
    Circle2D circle(Vector2D(3, 1), 1); // tangent at x=2

    CHECK(intersect(bbox, circle) == true);
}

TEST_CASE("Zero-size bbox")
{
    BBox2D b(Vector2D(1, 1), Vector2D(1, 1)); // point
    Circle2D c(Vector2D(1, 1), real_t{0.1});

    CHECK(intersect(b, c) == true);
}

TEST_CASE("Circle inside triangle")
{
    Triangle2D tri(Vector2D(0, 0), Vector2D(4, 0), Vector2D(2, 4));
    Circle2D c(Vector2D(2, 1), real_t{0.5});

    CHECK(intersect(tri, c) == true);
}

TEST_CASE("IPolygonalShape2D intersection")
{
    Triangle2D t1(Vector2D(0, 0), Vector2D(2, 0), Vector2D(1, 2));
    Triangle2D t2(Vector2D(1, 1), Vector2D(3, 1), Vector2D(2, 3));
    CHECK(intersect(t1, t2) == true);
    CHECK(intersect(t2, t1) == true);

    Triangle2D t3(Vector2D(3, 3), Vector2D(4, 3), Vector2D(real_t{3.5}, 4));
    CHECK(intersect(t1, t3) == false);
    CHECK(intersect(t3, t1) == false);
}

TEST_CASE("Concave polygon intersection respects its notch")
{
    Polygon2D uShape({
        Vector2D(0, 0), Vector2D(4, 0), Vector2D(4, 4), Vector2D(3, 4),
        Vector2D(3, 1), Vector2D(1, 1), Vector2D(1, 4), Vector2D(0, 4)
    });
    Polygon2D polygonInNotch({
        Vector2D(real_t{1.5}, 2), Vector2D(real_t{2.5}, 2),
        Vector2D(real_t{2.5}, 3), Vector2D(real_t{1.5}, 3)
    });

    CHECK(intersect(uShape, polygonInNotch) == false);
}

TEST_CASE("Degenerate triangle (collinear)")
{
    Triangle2D t(Vector2D(0, 0), Vector2D(1, 1), Vector2D(2, 2)); // line
    Triangle2D t2(Vector2D(0, 1), Vector2D(1, 2), Vector2D(2, 3));

    CHECK(intersect(t, t2) == false);
}

TEST_CASE("IPolygonalShape2D and Circle2D intersection")
{
    Triangle2D triangle(Vector2D(0, 0), Vector2D(2, 0), Vector2D(1, 2));
    Circle2D circle(Vector2D(1, 1), real_t{0.5});
    CHECK(intersect(triangle, circle) == true);

    Circle2D circle2(Vector2D(5, 5), 1);
    CHECK(intersect(triangle, circle2) == false);
}

TEST_CASE("Circle2D intersection")
{
    Circle2D c1(Vector2D(0, 0), 1);
    Circle2D c2(Vector2D(1, 0), 1);
    CHECK(intersect(c1, c2) == true);

    Circle2D c3(Vector2D(3, 0), 1);
    CHECK(intersect(c1, c3) == false);
}

/*
TEST_CASE("IFiniteShape2D intersection")
{
    Triangle2D t(Vector2D(0, 0), Vector2D(2, 0), Vector2D(1, 2));
    Circle2D c(Vector2D(1, 1), real_t{0.5});
    CHECK(intersect(static_cast<const IFiniteShape2D&>(t), static_cast<const IFiniteShape2D&>(c)) == true);

    Circle2D c2(Vector2D(5, 5), 1);
    CHECK(intersect(static_cast<const IFiniteShape2D&>(t), static_cast<const IFiniteShape2D&>(c2)) == false);
}
*/

//Todo: Testing alternative dispatching designs
TEST_CASE("FiniteShapeVariant2D intersection")
{
    Triangle2D triangle(Vector2D(0, 0), Vector2D(2, 0), Vector2D(1, 2));
    Rectangle2D rectangle({1, 1}, {1, -1}, {3, -1}, {3, 1});
    Circle2D intersectingCircle(Vector2D(1, 1), real_t{0.5});
    Circle2D overlappingCircle(Vector2D(1.5, 1), real_t{0.5});
    Circle2D separateCircle(Vector2D(5, 5), 1);

    const FiniteShapeVariant2D triangleShape{triangle};
    const FiniteShapeVariant2D rectangleShape{rectangle};
    const FiniteShapeVariant2D intersectingCircleShape{intersectingCircle};
    const FiniteShapeVariant2D overlappingCircleShape{overlappingCircle};
    const FiniteShapeVariant2D separateCircleShape{separateCircle};

    CHECK(intersect(triangleShape, rectangleShape));
    CHECK(intersect(rectangleShape, triangleShape));
    CHECK(intersect(triangleShape, intersectingCircleShape));
    CHECK(intersect(intersectingCircleShape, triangleShape));
    CHECK(intersect(intersectingCircleShape, overlappingCircleShape));
    CHECK_FALSE(intersect(triangleShape, separateCircleShape));
    CHECK_FALSE(intersect(intersectingCircleShape, separateCircleShape));
}

TEST_CASE("Ray2D and Segment2D intersection with FiniteShapeVariant2D")
{
    Triangle2D triangle(Vector2D(0, 0), Vector2D(2, 0), Vector2D(1, 2));
    const FiniteShapeVariant2D triangleShape{triangle};
    const FiniteShapeVariant2D circleShape{Circle2D(Vector2D(1, 1), real_t{0.5})};
    const FiniteShapeVariant2D distantCircleShape{Circle2D(Vector2D(5, 5), 1)};
    const Ray2D ray(Vector2D(-1, 1), Vector2D(1, 0));
    const Segment2D segment(Vector2D(-1, 1), Vector2D(3, 1));

    CHECK(intersect(ray, triangleShape));
    CHECK(intersect(ray, circleShape));
    CHECK_FALSE(intersect(ray, distantCircleShape));
    CHECK(intersect(segment, triangleShape));
    CHECK(intersect(segment, circleShape));
    CHECK_FALSE(intersect(segment, distantCircleShape));
}