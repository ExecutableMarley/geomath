/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Distance3D.hpp"
#include "../Line3D.hpp"

namespace Arns::geomath
{

// Based on David Eberly's robust segment-segment distance algorithm.
// References: 
// - Geometric Tools (https://www.geometrictools.com)
// - "3D Game Engine Design" by David Eberly

real_t detail::distanceSegmentToSegment(
    const Vector3D& s1, const Vector3D& s2, 
    const Vector3D& k1, const Vector3D& k2,
    Vector3D* closestPoint1, Vector3D* closestPoint2)
{
    const Vector3D d1 = s2 - s1; // Direction vector of segment 1
    const Vector3D d2 = k2 - k1; // Direction vector of segment 2
    const Vector3D r  = s1 - k1;

    const real_t a = d1.lengthSquared(); // Squared length of segment 1
    const real_t e = d2.lengthSquared(); // Squared length of segment 2
    const real_t f = d2.dot(r);

    const bool deg1 = approximatelyZero(a);
    const bool deg2 = approximatelyZero(e);

    real_t s = real_t{0};
    real_t t = real_t{0};

    // Check if either or both segments degenerate into points
    if (deg1 && deg2) {
        s = real_t{0};
        t = real_t{0};
    } 
    else if (deg1) {
        s = real_t{0};
        t = clamp(f / e, real_t{0}, real_t{1});
    } 
    else {
        const real_t c = d1.dot(r);
        if (deg2) {
            t = real_t{0};
            s = clamp(-c / a, real_t{0}, real_t{1});
        } 
        else {
            // General non-degenerate case
            const real_t b = d1.dot(d2);
            const real_t denom = a * e - b * b;

            // Compute unconstrained s if segments are not parallel
            if (!approximatelyZero(denom)) {
                s = clamp((b * f - c * e) / denom, real_t{0}, real_t{1});
            } else {
                s = real_t{0}; // Segments are parallel
            }

            // Compute optimal t for the chosen s
            t = (b * s + f) / e;

            // If t lies outside [0, 1], clamp t and recompute s on that boundary
            if (approximatelyLess(t, real_t{0})) {
                t = real_t{0};
                s = clamp(-c / a, real_t{0}, real_t{1});
            } else if (approximatelyGreater(t, real_t{1})) {
                t = real_t{1};
                s = clamp((b - c) / a, real_t{0}, real_t{1});
            }
        }
    }

    const Vector3D cp1 = s1 + d1 * s;
    const Vector3D cp2 = k1 + d2 * t;

    if (closestPoint1) *closestPoint1 = cp1;
    if (closestPoint2) *closestPoint2 = cp2;

    return (cp1 - cp2).length();
}

real_t detail::distancePointToSegment(const Vector3D& point, const Vector3D& lineStart, const Vector3D& lineEnd, Vector3D* closestPoint)
{
    const Vector3D direction = lineEnd - lineStart;
    const real_t directionLengthSquared = direction.lengthSquared();
    const real_t parameter = directionLengthSquared > real_t{0}
                                 ? clamp((point - lineStart).dot(direction) / directionLengthSquared, real_t{0}, real_t{1})
                                 : real_t{0};
    const Vector3D closest = lineStart + direction * parameter;

    if (closestPoint)
        *closestPoint = closest;

    return (point - closest).length();
}

real_t detail::distancePointToSegment(const Vector3D& point, const Segment3D& line, Vector3D* closestPoint)
{
    return distancePointToSegment(point, line.m_start, line.m_end, closestPoint);
}

real_t detail::distanceLineToLine(const Vector3D& s1, const Vector3D& s2, const Vector3D& k1, const Vector3D& k2,
                          Vector3D* closestPoint1, Vector3D* closestPoint2)
{
    return distanceSegmentToSegment(s1, s2, k1, k2, closestPoint1, closestPoint2);
}

real_t detail::distanceLineToLine(const Segment3D& line1, const Segment3D& line2,
                          Vector3D* closestPoint1, Vector3D* closestPoint2)
{
    return distanceLineToLine(line1.m_start, line1.m_end, line2.m_start, line2.m_end, closestPoint1, closestPoint2);
}

}