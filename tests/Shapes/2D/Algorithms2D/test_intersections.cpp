#include "../../../third_party/doctest.h"

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "Shapes/2D/Ray2D.hpp"
#include "Shapes/2D/Line2D.hpp"
#include "Shapes/2D/Triangle2D.hpp"
#include "Shapes/2D/Rectangle2D.hpp"
#include "Shapes/2D/Circle2D.hpp"
#include "Shapes/2D/Algorithms2D/Algorithms2D.hpp"

using namespace Arns::Math;


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
        checkHit(ray, box, {true, {7.0f, 0.0f}, 7.0f});
    }

    SUBCASE("Triangle")
    {
        Triangle2D tri({15, -1}, {18, 1}, {18, -1});
        checkHit(ray, tri, {true, {16.5f, 0.0f}, 16.5f});
    }

    SUBCASE("Rectangle")
    {
        auto rect = Rectangle2D::fromMinMax({10.0f, -1.0f}, {14.0f, 1.0f});
        checkHit(ray, rect, {true, {10.0f, 0.0f}, 10.0f});
    }

    SUBCASE("Circle")
    {
        Circle2D circle({5, 0}, 1);
        checkHit(ray, circle, {true, {4, 0}, 4});
    }
}

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
        auto tri = Triangle2D({-2.0f, -2.0f}, {2.0f, -2.0f}, {2.0f, 3.0f});
        checkHit(ray, tri, {true, {2, 0}, 2});
    }

    SUBCASE("Rectangle containing ray origin")
    {
        auto rect = Rectangle2D::fromMinMax({-1.0f, -1.0f}, {4.0f, 1.0f});
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
        auto rect = BBox2D({8.0f, 0.0f}, {12.0f, 1.0f});
        checkHit(ray, rect, {true, {8, 1}, 8});
    }

    SUBCASE("Rectangle tangent to ray")
    {
        auto rect = Rectangle2D::fromMinMax({8.0f, 0.0f}, {12.0f, 1.0f});
        checkHit(ray, rect, {true, {8, 1}, 8});
    }

    SUBCASE("Triangle tangent to ray")
    {
        Triangle2D tri({20, 1}, {23, 3}, {23, 1});
        checkHit(ray, tri, {true, {20, 1}, 20});
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
        {6.0f, 0.0f}, 
        6.0f
    });

    Triangle2D tri({3, 1}, {3, -1}, {6, -1});
    checkHit(ray, tri, {
    true, 
    {3.0f, 0.0f}, 
    3.0f
    });

    Circle2D circ({11.0f, 0.0f}, 1.0f);
    checkHit(ray, circ, {
        true, 
        {10.0f, 0.0f}, 
        10.0f
    });

    Rectangle2D rect = Rectangle2D::fromMinMax({5.0f, -1.0f}, {10.0f, 1.0f}); 
    checkHit(ray, rect, {
        true, 
        {5.0f, 0.0f}, 
        5.0f
    });
}