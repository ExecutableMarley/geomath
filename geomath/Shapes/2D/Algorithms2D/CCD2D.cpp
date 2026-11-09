/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms2D.hpp"
#include "CommonMath.hpp"
#include "../Ray2D.hpp"
#include "../Segment2D.hpp"
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

// Todo: This still occasionally detects separating and minimal overlapping circles occasionally
SweepResult2D detail::SolveCircleToCircleCCD(const Circle2D& a, const Vector2D& velA,
    const Circle2D& b, const Vector2D& velB, real_t maxTime) 
{
    SweepResult2D result;

    Vector2D p0 = a.m_center - b.m_center;
    Vector2D v = velA - velB;
    real_t R = a.m_radius + b.m_radius;

    real_t c = p0.dot(p0) - R * R;
    real_t approach = p0.dot(v); // Half of b_coeff
    bool shouldReportImmediateHit = false;

    if (approximatelyLess(c, real_t{0})) // Geometric overlap on TOI = 0
    {
        shouldReportImmediateHit = true;
    }
    else if (!approximatelyGreaterAbs(c, real_t{0})) // Touching at Toi = 0
    {
        // Check if approaching
        if (approximatelyLess(approach, real_t{0}))
        {
            shouldReportImmediateHit = true;
        }
    }

    if (shouldReportImmediateHit)
    {
        result.hit = true;
        result.toi = real_t{0};
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

SweepResult2D detail::SolveBBoxToBBoxCCD(const BBox2D& a, const Vector2D& velA,
    const BBox2D& b, const Vector2D& velB, real_t maxTime) 
{
    SweepResult2D result;

    Vector2D p0 = a.centroid() - b.centroid();
    Vector2D E = a.halfExtents() + b.halfExtents();
    Vector2D v = velA - velB;

    bool isInsideX = approximatelyBetween(p0.x, -E.x, E.x);
    bool isInsideY = approximatelyBetween(p0.y, -E.y, E.y);

    if (isInsideX && isInsideY)
    {
        bool isStrictlyInside = (p0.x > -E.x && p0.x < E.x) &&
                                (p0.y > -E.y && p0.y < E.y);

        if (!isStrictlyInside)
        {
            bool separating = false;

            if (approximatelyZero(std::abs(std::abs(p0.x) - E.x)))
            {
                if ((p0.x > real_t{0} && v.x > real_t{0}) || (p0.x < real_t{0} && v.x < real_t{0}))
                    separating = true;
            }
            if (approximatelyZero(std::abs(std::abs(p0.y) - E.y)))
            {
                if ((p0.y > real_t{0} && v.y > real_t{0}) || (p0.y < real_t{0} && v.y < real_t{0}))
                    separating = true;
            }

            if (separating)
            {
                return result; // Moving apart
            }
        }

        result.hit = true;
        result.toi = real_t{0};
        Vector2D pen = E - p0;
        if (pen.x < pen.y)
        {
            result.normal = (p0.x < real_t{0}) ? Vector2D{-1, 0} : Vector2D{1, 0};
        }
        else
        {
            result.normal = (p0.y < real_t{0}) ? Vector2D{0, -1} : Vector2D{0, 1};
        }
        result.point = b.closestPoint(a.centroid());
        return result;
    }

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

    /*
    // Check if shapes were already overlapping at t = 0
    if (approximatelyBetween(p0.x, -E.x, E.x) && approximatelyBetween(p0.y, -E.y, E.y))
    {
        result.hit = true;
        result.toi = 0;
        result.normal = hitNormal; 
        result.point = b.closestPoint(a.centroid());
        return result;
    }
    */

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