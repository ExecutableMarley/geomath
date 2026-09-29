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

#include <cassert>
#include <span>
#include <type_traits>

namespace Arns::geomath
{

//Post impact velocity response?

//Todo: Refactor the specialized versions at some point to fit them in

/*
struct CCDResult 
{
    bool hit{false};
    real_t toi{0.0}; // Time of Impact in range [0, 1]
    Vector2D normal{0.0, 0.0}; // Contact normal pointing from B to A
};

CCDResult SolveCircleToCircleCCD(const Circle2D& a, const Vector2D& velA,
    const Circle2D& b, const Vector2D& velB, real_t maxTime) 
{
    CCDResult result;

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
        result.normal = approximatelyZero(dist) ? Vector2D{1, 0} : p0 / dist;
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
        result.normal = hitPos.normalize(); // Points away from center B
    }

    return result;
}

CCDResult SolveBBoxToBBoxCCD(const BBox2D& a, const Vector2D& velA,
    const BBox2D& b, const Vector2D& velB, real_t maxTime) 
{
    CCDResult result;

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
    //if (p0.x >= -E.x && p0.x <= E.x && p0.y >= -E.y && p0.y <= E.y) 
    if (approximatelyBetween(p0.x, -E.x, E.x) && approximatelyBetween(p0.y, -E.y, E.y))
    {
        result.hit = true;
        result.toi = 0;
        result.normal = hitNormal; 
        return result;
    }

    //if (tMin >= 0 && tMin <= 1)
    if (approximatelyBetween(tMin, real_t{0}, real_t{1}))
    {
        result.hit = true;
        result.toi = tMin;
        result.normal = hitNormal;
    }

    return result;
}

CCDResult SolveBBoxToCircleCCD(
    const BBox2D& box, const Vector2D& velBox,
    const Circle2D& circle, const Vector2D& velCircle,
    real_t maxTime, int maxIterations, real_t epsilon) 
{
    CCDResult result;

    // Compute relative velocity (Circle moving relative to a stationary AABB)
    Vector2D relVel = velCircle - velBox;
    real_t speed = relVel.length();

    // Static case: circle is stationary relative to the box
    if (speed < real_t{1e-8}) {
        real_t dist = box.minDistance(circle.m_center) - circle.m_radius;
        if (dist <= real_t{0}) {
            result.hit = true;
            result.toi = real_t{0};
            Vector2D cp = box.closestPoint(circle.m_center);
            Vector2D normalDir = box.centroid() - circle.m_center; // fallback vector
            Vector2D diff = cp - circle.m_center;
            result.normal = (diff.lengthSquared() > real_t{1e-8}) ? diff.normalize() : normalDir.normalize();
        }
        return result;
    }

    real_t t = real_t{0};

    // Conservative Advancement Iterative Loop
    for (int i = 0; i < maxIterations; ++i) {
        Vector2D currentCircleCenter = circle.m_center + relVel * t;

        // Separation distance between circle and AABB
        real_t boxDist = box.minDistance(currentCircleCenter);
        real_t separationDist = boxDist - circle.m_radius;

        // Collision threshold reached
        if (separationDist <= epsilon) {
            result.hit = true;
            result.toi = t;

            // Calculate contact normal at TOI
            Vector2D cp = box.closestPoint(currentCircleCenter);
            Vector2D diff = cp - currentCircleCenter;

            // Normal points from Circle to Box
            if (diff.lengthSquared() > real_t{1e-8})
            {
                result.normal = diff.normalize();
            }
            else
            {
                // Center is inside the box at contact (fallback normal opposite to velocity)
                result.normal = (relVel * real_t{-1}).normalize();
            }
            return result;
        }

        // Advance time conservatively along relative velocity vector
        real_t dt = separationDist / speed;
        t += dt;

        // Exceeded max time
        if (t > maxTime)
        {
            return result;
        }
    }

    return result;
}
*/

template <typename A, typename B>
SweepResult sweepConservativeAdvancementGeneric(
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

        if constexpr (std::is_base_of_v<IFiniteShape2D, A> &&
                      std::is_base_of_v<IFiniteShape2D, B>)
        {
            std::unique_ptr<IFiniteShape2D> currentB = b.clone();
            currentB->translate(relativeVelocity * t);
            d = distance(a, *currentB, &cloA, &cloB);
        }
        else
        {
            B currentB = b;
            currentB.translate(relativeVelocity * t);
            d = distance(a, currentB, &cloA, &cloB);
        }

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

SweepResult sweepConservativeAdvancement(
    const IFiniteShape2D& a, Vector2D velocityA,
    const IFiniteShape2D& b, Vector2D velocityB,
    real_t maxTime, int maxIterations, real_t epsilon)
{
    //const Circle2D* circleA = shape_cast<Circle2D>(&a);
    //const Circle2D* circleB = shape_cast<Circle2D>(&b);

    return sweepConservativeAdvancementGeneric(a, velocityA, b, velocityB, maxTime, maxIterations, epsilon);
}


// Todo: Consider the correct place for this

BBox2D sweptBounds(const BBox2D& bounds,
                   const Vector2D& displacement)
{
    return {
        Vector2D::min(bounds.m_min, bounds.m_min + displacement),
        Vector2D::max(bounds.m_max, bounds.m_max + displacement)};
}

BBox2D sweptBounds(const IFiniteShape2D& shape,
                   Vector2D velocity,
                   real_t dt)
{
    return sweptBounds(shape.boundingBox(), velocity * dt);
}


}