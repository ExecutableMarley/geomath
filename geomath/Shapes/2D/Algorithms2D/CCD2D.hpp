/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Geometry/Vector2D.hpp"
#include "Distance2D.hpp"

namespace Arns::geomath
{
    
struct SweepResult2D
{
    bool hit = false;
    real_t toi = 1.0;
    Vector2D normal{}; // a to b
    Vector2D point{};
};

namespace detail
{
    SweepResult2D SolveCircleToCircleCCD(const Circle2D& a, const Vector2D& velA,
        const Circle2D& b, const Vector2D& velB, real_t maxTime);

    SweepResult2D SolveBBoxToBBoxCCD(const BBox2D& a, const Vector2D& velA,
        const BBox2D& b, const Vector2D& velB, real_t maxTime);
}

template <typename A, typename Object>
SweepResult2D sweepConservativeAdvancementGeneric(
    const A& a, Vector2D velocityA,
    const Object& b, Vector2D velocityB,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4))
{
    const Vector2D relativeVelocity = velocityB - velocityA;
    const real_t speed = relativeVelocity.length();

    // Note: This does not detect deep penetration. We might ignore things that should be handled 
    // As collision.
    // Todo: Proper DCD with MTV to reuse here

    Vector2D cloA, cloB;
    real_t d = distance(a, b, &cloA, &cloB);
    if (d <= epsilon)
    {
        Vector2D normalDelta = cloA - cloB;
        const real_t normalLengthSquared = normalDelta.lengthSquared();
        const Vector2D normal = normalLengthSquared > epsilon * epsilon
            ? normalDelta.normalize()
            : speed > epsilon ? (-relativeVelocity).normalize() : Vector2D{1, 0};

        // Check if they are separating at t = 0
        real_t separationSpeed = normal.dot(relativeVelocity);

        if (separationSpeed >= -epsilon)
        {
            // Moving apart -> No Collision
            return {};
        }

        // Approaching
        return {
            .hit = true,
            .toi = real_t{0},
            .normal = normal,
            .point = cloB
        };
    }

    if (speed <= epsilon)
    {
        return {};
    }

    real_t t = 0;
    for (int i = 0; i < maxIterations; ++i)
    {
        t += d / speed;
        if (t > maxTime)
            break;

        Vector2D cloA, cloB;
        real_t d;

        Object currentB = b;
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

        //if (t > maxTime)
          //  break;
    }

    return {};
}


template <typename A, typename Object, typename Enable = void>
    requires detail::HasDistance<A, Object>
struct SweepTraits {
    static SweepResult2D sweep(
        const A& a, Vector2D va,
        const Object& b, Vector2D vb,
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
        return detail::SolveCircleToCircleCCD(a, va, b, vb, maxTime);
    }
};

//================================
// Public CCD dispatch function 
//================================

// Helper function
template <typename A, typename Object>
    requires detail::HasDistance<A, Object>
SweepResult2D sweep(
    const A& a, Vector2D va,
    const Object& b, Vector2D vb,
    real_t maxTime = real_t(1), int maxIterations = 15, real_t epsilon = real_t(1e-4))
{
    return SweepTraits<A, Object>::sweep(a, va, b, vb, maxTime, maxIterations, epsilon);
}


} // namespace Arns::geomath
