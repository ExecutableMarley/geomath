/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms2D.hpp"
#include "CommonMath.hpp"
#include "../Ray2D.hpp"
#include "../Line2D.hpp"
#include "../BBox2D.hpp"
#include "../Circle2D.hpp"
#include "../Triangle2D.hpp"
#include "../Rectangle2D.hpp"
#include "../Polygon2D.hpp"
#include "../ConvexPolygon2D.hpp"
#include "Shapes/2D/Variants/ShapeConcepts.hpp"
#include "Shapes/2D/Variants/ShapeVariants.hpp"
#include "ShapeVariantAlgorithms.hpp"

#include <cassert>
#include <span>
#include <type_traits>

namespace Arns::geomath
{

//Post impact velocity response?

//Todo: Refactor the specialized versions at some point to fit them in


SweepResult2D SolveCircleToCircleCCD(const Circle2D& a, const Vector2D& velA,
    const Circle2D& b, const Vector2D& velB, real_t maxTime) 
{
    SweepResult2D result;

    Vector2D p0 = a.m_center - b.m_center;
    Vector2D v = velA - velB;
    real_t R = a.m_radius + b.m_radius;

    real_t c = p0.dot(p0) - R * R;

    // Initial overlap check (t = 0)
    if (!approximatelyGreaterAbs(c, real_t{0})) 
    {
        result.hit = true;
        result.toi = 0;
        real_t dist = p0.length();
        result.normal = approximatelyZero(dist) ? Vector2D{1, 0} : -p0 / dist;
        result.point = b.m_center - result.normal * b.m_radius;
        return result;
    }

    real_t a_coeff = dot(v, v);
    if (approximatelyZero(a_coeff)) // Relative velocity is zero or negligible
    {
        return result;
    }

    real_t b_coeff = 2 * dot(p0, v);
    if (approximatelyGreater(b_coeff, real_t{0})) // Moving away from each other
    {
        return result;
    }

    real_t discriminant = b_coeff * b_coeff - 4 * a_coeff * c;
    if (discriminant < real_t{0})
    {
        if (!approximatelyZero(discriminant))
            return result;

        discriminant = real_t{0};
    }

    // Smallest root gives the entry time
    real_t t = (-b_coeff - std::sqrt(discriminant)) / (2 * a_coeff);

    if (approximatelyBetween(t, real_t{0}, maxTime))
    {
        result.hit = true;
        result.toi = t;
        Vector2D hitPos = p0 + v * t;
        result.normal = -hitPos.normalize();
        result.point = b.m_center + velB * t - result.normal * b.m_radius;
    }

    return result;
}

SweepResult2D SolveBBoxToBBoxCCD(const BBox2D& a, const Vector2D& velA,
    const BBox2D& b, const Vector2D& velB, real_t maxTime) 
{
    SweepResult2D result;

    Vector2D p0 = a.centroid() - b.centroid();
    Vector2D E = a.halfExtents() + b.halfExtents();
    Vector2D v = velA - velB;

    real_t tMin = real_t{0};
    real_t tMax = maxTime;
    Vector2D hitNormal{0, 0};

    // Slab test for X and Y axes
    for (int i = 0; i < 2; ++i)
    {
        real_t pos = p0[i];
        real_t extent = E[i];
        real_t vel = v[i];

        if (approximatelyZero(vel)) 
        {
            // Velocity along axis is zero; check if origin is outside slab
            if (pos < -extent || pos > extent) {
                return result; // Parallel and outside
            }
        } 
        else 
        {
            real_t invV = real_t{1} / vel;
            real_t t1 = (-extent - pos) * invV;
            real_t t2 = ( extent - pos) * invV;

            Vector2D axisNormal{0, 0};
            if (vel > real_t{0.0}) {
                axisNormal[i] = real_t{-1.0}; // Hitting left side of minkowski box
            } else {
                axisNormal[i] = real_t{1.0};  // Hitting right side of minkowski box
            }

            if (approximatelyGreater(t1, t2)) {
                std::swap(t1, t2);
                axisNormal[i] = -axisNormal[i];
            }

            if (approximatelyGreater(t1, tMin)) {
                tMin = t1;
                hitNormal = axisNormal;
            }

            tMax = std::min(tMax, t2);

            if (approximatelyGreater(tMin, tMax)) {
                return result; // Misses the box entirely
            }
        }
    }

    // Check if shapes were already overlapping at t = 0
    if (approximatelyBetween(p0.x, -E.x, E.x) && approximatelyBetween(p0.y, -E.y, E.y))
    {
        result.hit = true;
        result.toi = 0;
        result.normal = hitNormal; 
        result.point = b.closestPoint(a.centroid());
        return result;
    }

    if (approximatelyBetween(tMin, real_t{0}, maxTime))
    {
        result.hit = true;
        result.toi = tMin;
        result.normal = hitNormal;
        const Vector2D centerA = a.centroid() + velA * tMin;
        const BBox2D boxBAtImpact = b.copy().translate(velB * tMin);
        if (!approximatelyZero(hitNormal.x))
        {
            result.point = {
                hitNormal.x > real_t{0} ? boxBAtImpact.m_max.x : boxBAtImpact.m_min.x,
                std::clamp(centerA.y, boxBAtImpact.m_min.y, boxBAtImpact.m_max.y)};
        }
        else
        {
            result.point = {
                std::clamp(centerA.x, boxBAtImpact.m_min.x, boxBAtImpact.m_max.x),
                hitNormal.y > real_t{0} ? boxBAtImpact.m_max.y : boxBAtImpact.m_min.y};
        }
    }

    return result;
}

template <typename A, typename B>
SweepResult2D sweepConservativeAdvancementGeneric(
    const A& a, Vector2D velocityA,
    const B& b, Vector2D velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4))
{
    const Vector2D relativeVelocity = velocityB - velocityA;
    const real_t speed = relativeVelocity.length();

    real_t t = 0;

    for (int i = 0; i < maxIterations; ++i)
    {
        if (t > maxTime)
            break;

        Vector2D cloA, cloB;
        real_t d;

        B currentB = b;
        currentB.translate(relativeVelocity * t);
        d = distance(a, currentB, &cloA, &cloB);

        if (d <= epsilon)
        {
            Vector2D normalDelta = cloA - cloB;
            const real_t normalLengthSquared = normalDelta.lengthSquared();
            const Vector2D normal = normalLengthSquared > epsilon * epsilon
                ? normalDelta.normalize()
                : speed > epsilon ? (-relativeVelocity).normalize() : Vector2D{1, 0};

            return {
                .hit = true,
                .toi = t,
                .normal = normal,
                .point = cloB
            };
        }

        if (speed <= epsilon)
            break;

        t += d / speed;

        if (t > maxTime)
            break;
    }

    return {};
}

template <typename A, typename B, typename Enable = void>
struct SweepTraits {
    static SweepResult2D sweep(
        const A& a, Vector2D va,
        const B& b, Vector2D vb,
        real_t maxTime, int maxIterations, real_t epsilon) 
    {
        return sweepConservativeAdvancementGeneric(a, va, b, vb, maxTime, maxIterations, epsilon);
    }
};

template <>
struct SweepTraits<Circle2D, Circle2D> {
    static SweepResult2D sweep(
        const Circle2D& a, Vector2D va,
        const Circle2D& b, Vector2D vb,
        real_t maxTime, int, real_t)
    {
        return SolveCircleToCircleCCD(a, va, b, vb, maxTime);
    }
};

// Helper function
template <typename A, typename B>
SweepResult2D sweep(
    const A& a, Vector2D va,
    const B& b, Vector2D vb,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4))
{
    return SweepTraits<A, B>::sweep(a, va, b, vb, maxTime, maxIterations, epsilon);
}

/*
SweepResult2D sweep(
    const IFiniteShape2D& shapeA, Vector2D va,
    const IFiniteShape2D& shapeB, Vector2D vb,
    real_t maxTime, int maxIterations, real_t epsilon)
{
    const Circle2D* circle1 = shape_cast<Circle2D>(&shapeA);
    const Circle2D* circle2 = shape_cast<Circle2D>(&shapeB);
    const IPolygonalShape2D* polygon1 = shapeA.polygonal();
    const IPolygonalShape2D* polygon2 = shapeB.polygonal();

    if (circle1 && circle2)
        return SweepTraits<Circle2D, Circle2D>::sweep(*circle1, va, *circle2, vb, maxTime, maxIterations, epsilon);

    if (polygon1 && polygon2)
        return SweepTraits<IPolygonalShape2D, IPolygonalShape2D>::sweep(*polygon1, va, *polygon2, vb, maxTime, maxIterations, epsilon);

    if (polygon1 && circle2)
        return SweepTraits<IPolygonalShape2D, Circle2D>::sweep(*polygon1, va, *circle2, vb, maxTime, maxIterations, epsilon);

    if (circle1 && polygon2)
        return SweepTraits<Circle2D, IPolygonalShape2D>::sweep(*circle1, va, *polygon2, vb, maxTime, maxIterations, epsilon);

    return {};
}
*/

SweepResult2D sweep(
    const FiniteShapeVariant2D& shapeA, Vector2D va, 
    const FiniteShapeVariant2D& shapeB, Vector2D vb,
    real_t maxTime, int maxIterations, real_t epsilon)
{
    return std::visit([&](const auto& a, const auto& b) {
        return sweep(a, va, b, vb, maxTime, maxIterations, epsilon);
    }, shapeA, shapeB);
}

}