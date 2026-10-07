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
#include "../Variants/ShapeConcepts.hpp"
#include "../Variants/ShapeVariants.hpp"

#include <cassert>
#include <span>

namespace Arns::geomath
{



namespace detail
{

real_t distancePointToSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd, Vector2D* closestPoint)
{
    const Vector2D r = segmentEnd - segmentStart;
    const Vector2D s = point - segmentStart;
    const real_t rDotR = r.dot(r);
    if (rDotR == 0.0f)
    {
        if (closestPoint) *closestPoint = segmentStart;
        return (point - segmentStart).length();
    }
    real_t t = std::clamp(s.dot(r) / rDotR, real_t{0}, real_t{1});
    const Vector2D projection = segmentStart + r * t;
    if (closestPoint)
        *closestPoint = projection;
    return (point - projection).length();
}

// Based on David Eberly's robust segment-segment distance algorithm.
// References: 
// - Geometric Tools (https://www.geometrictools.com)
// - "3D Game Engine Design" by David Eberly

real_t distanceSegmentToSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, 
    Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    // Delta Vectors
    const Vector2D d1 = s2 - s1;
    const Vector2D d2 = k2 - k1;
    const Vector2D r  = s1 - k1;

    // Delta lengths
    const real_t l1 = d1.lengthSquared();
    const real_t l2 = d2.lengthSquared();
    const real_t f = d2.dot(r);

    const bool degen1 = approximatelyZero(l1);
    const bool degen2 = approximatelyZero(l2);

    real_t s = real_t{0};
    real_t t = real_t{0};

    if (degen1 && degen2)
    {
        s = real_t{0};
        t = real_t{0};
    }
    else if (degen1)
    {
        s = real_t{0};
        t = clamp(f / l2, real_t{0}, real_t{1});
    }
    else
    {
        const real_t c = d1.dot(r);
        if (degen2) 
        {
            t = real_t{0};
            s = clamp(-c / l1, real_t{0}, real_t{1});
        } 
        else 
        {
            // General non-degenerate case
            const real_t b = d1.dot(d2);
            const real_t denom = l1 * l2 - b * b;

            // Compute unconstrained s
            if (!approximatelyZero(denom)) 
            {
                s = clamp((b * f - c * l2) / denom, real_t{0}, real_t{1});
            } 
            else
            {
                s = real_t{0}; // parallel
            }

            // Compute optimal t for the chosen s
            t = (b * s + f) / l2;

            // If t lies outside [0, 1], clamp t and recompute s on that boundary
            if (approximatelyLess(t, real_t{0})) 
            {
                t = real_t{0};
                s = clamp(-c / l1, real_t{0}, real_t{1});
            } 
            else if (approximatelyGreater(t, real_t{1})) 
            {
                t = real_t{1};
                s = clamp((b - c) / l1, real_t{0}, real_t{1});
            }
        }
    }

    const Vector2D best1 = s1 + d1 * s;
    const Vector2D best2 = k1 + d2 * t;

    if (closestPoint1) *closestPoint1 = best1;
    if (closestPoint2) *closestPoint2 = best2;

    return (best1 - best2).length();
}

real_t distancePointToPolygonEdges(const Vector2D& point, const std::span<const Vector2D> polygon, Vector2D* closestPoint)
{
    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D best;

    for (size_t i = 0; i < polygon.size(); ++i)
    {
        const Vector2D& a = polygon[i];
        const Vector2D& b = polygon[(i + 1) % polygon.size()];

        Vector2D q;
        real_t dist = distancePointToSegment(point, a, b, &q);

        if (dist < minDist)
        {
            minDist = dist;
            best = q;
        }
    }

    if (closestPoint) *closestPoint = best;
    return minDist;
}

real_t distanceSegmentToEdgeList(const Segment2D& seg, std::span<const Vector2D> vertices, bool closed, Vector2D* closestSeg, Vector2D* closestEdge)
{
    const size_t n = vertices.size();
    if (n < 2)
        return std::numeric_limits<real_t>::max();

    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D bestSeg;
    Vector2D bestEdge;

    const size_t edgeCount = closed ? n : n - 1;

    for (size_t i = 0; i < edgeCount; ++i)
    {
        const Vector2D& a = vertices[i];
        const Vector2D& b = vertices[(i + 1) % n];

        Vector2D cs, ce;
        real_t d = distanceSegmentToSegment(seg.m_start, seg.m_end, a, b, &cs, &ce);

        if (d < minDist)
        {
            minDist  = d;
            bestSeg  = cs;
            bestEdge = ce;

            // => Intersection
            if (minDist == real_t(0))
                break;
        }
    }

    if (closestSeg)  *closestSeg  = bestSeg;
    if (closestEdge) *closestEdge = bestEdge;

    return minDist;
}

}


// Distance calculation algorithms

real_t distance(const Vector2D& point, const Segment2D& segment, Vector2D* closestPoint)
{
    return detail::distancePointToSegment(point, segment.m_start, segment.m_end, closestPoint);
}

real_t distance(const Segment2D& segment1, const Segment2D& segment2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return detail::distanceSegmentToSegment(segment1.m_start, segment1.m_end, segment2.m_start, segment2.m_end, closestPoint1, closestPoint2);
}

//Todo: Distance to shapes

real_t distance(const Vector2D& point, const BBox2D& bbox, Vector2D* closestPoint)
{
    Vector2D c = bbox.closestPoint(point);
    if (closestPoint) *closestPoint = c;
    return point.distance(c);
}

//Todo: Consider renaming slightly



real_t distance(const Vector2D& point, const IPolygonalShape2D& polygon, Vector2D* closestPoint)
{
    //Todo: refactor this
    if (isConvexPolygonal(polygon.type()))
    {
        if (detail::isPointInsideConvexPolygon(point, polygon.vertices()))
        {
            if (closestPoint)
                *closestPoint = point;
            return real_t{0};
        }
        return detail::distancePointToPolygonEdges(point, polygon.vertices(), closestPoint);
    }
    else
    {
        if (detail::isPointInsidePolygon(point, polygon.vertices()))
        {
            if (closestPoint)
                *closestPoint = point;
            return real_t{0};
        }
        return detail::distancePointToPolygonEdges(point, polygon.vertices(), closestPoint);
    }
}

real_t distance(const Vector2D& point, const Circle2D& circle, Vector2D* closestPoint)
{
    Vector2D delta = point - circle.m_center;
    real_t dist = delta.length();

    if (dist <= circle.m_radius)
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }

    //real_t dist = std::sqrt(distSq);
    Vector2D dir = delta / dist;

    if (closestPoint)
        *closestPoint = circle.m_center + dir * circle.m_radius;

    return std::abs(dist - circle.m_radius);
}

real_t boundaryDistance(const Segment2D& segment, const BBox2D& bbox, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return detail::distanceSegmentToEdgeList(segment, bbox.corners(), true, closestPoint1, closestPoint2);
}

real_t distance(const Segment2D& segment, const BBox2D& bbox, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    if (isPointInsideBBox(segment.m_start, bbox))
    {
        if (closestPoint1) *closestPoint1 = segment.m_start;
        if (closestPoint2) *closestPoint2 = segment.m_start;
        return real_t{0};
    }

    return detail::distanceSegmentToEdgeList(segment, bbox.corners(), true, closestPoint1, closestPoint2);
}

real_t boundaryDistance(const Segment2D& segment, const IPolygonalShape2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return detail::distanceSegmentToEdgeList(segment, polygon.vertices(), true, closestPoint1, closestPoint2);
}

real_t distance(const Segment2D& segment, const IPolygonalShape2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    if (isPointInsidePolygon(segment.m_start, polygon))
    {
        if (closestPoint1) *closestPoint1 = segment.m_start;
        if (closestPoint2) *closestPoint2 = segment.m_start;
        return real_t{0};
    }

    return detail::distanceSegmentToEdgeList(segment, polygon.vertices(), true, closestPoint1, closestPoint2);
}

real_t distance(const Segment2D& segment, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    Vector2D p;
    real_t dist = detail::distancePointToSegment(circle.m_center, segment.m_start, segment.m_end, &p);

    //Trivial case
    if (dist <= circle.m_radius)
    {
        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = p;
        return real_t(0);
    }

    //Not intersecting case
    Vector2D dir = p - circle.m_center;
    real_t len = dir.length();

    if (len > real_t{0})
        dir /= len;

    Vector2D q = circle.m_center + dir * circle.m_radius;

    if (closestPoint1) *closestPoint1 = p;
    if (closestPoint2) *closestPoint2 = q;

    return len - circle.m_radius;
}


real_t distance(const BBox2D& a, const BBox2D& b, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    real_t dx = 0;
    real_t dy = 0;

    Vector2D pa;
    Vector2D pb;

    // --- X axis ---
    if (a.m_max.x < b.m_min.x)
    {
        dx = b.m_min.x - a.m_max.x;
        pa.x = a.m_max.x;
        pb.x = b.m_min.x;
    }
    else if (b.m_max.x < a.m_min.x)
    {
        dx = a.m_min.x - b.m_max.x;
        pa.x = a.m_min.x;
        pb.x = b.m_max.x;
    }
    else
    {
        // overlap
        real_t x = std::max(a.m_min.x, b.m_min.x);
        pa.x = pb.x = x;
    }

    // --- Y axis ---
    if (a.m_max.y < b.m_min.y)
    {
        dy = b.m_min.y - a.m_max.y;
        pa.y = a.m_max.y;
        pb.y = b.m_min.y;
    }
    else if (b.m_max.y < a.m_min.y)
    {
        dy = a.m_min.y - b.m_max.y;
        pa.y = a.m_min.y;
        pb.y = b.m_max.y;
    }
    else
    {
        // overlap
        real_t y = std::max(a.m_min.y, b.m_min.y);
        pa.y = pb.y = y;
    }

    if (closestPoint1) *closestPoint1 = pa;
    if (closestPoint2) *closestPoint2 = pb;

    return std::sqrt(dx * dx + dy * dy);
}

//real_t distanceCircleToCircle(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t boundaryDistance(const BBox2D& bbox1, const IPolygonalShape2D& poly2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    //Todo: This can (and should) be optimized

    const auto v = poly2.vertices();
    const size_t n = v.size();

    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D bestP1, bestP2;
    
    auto b = bbox1.corners();

    for (int i = 0; i < 4; ++i)
    {
        Vector2D a1 = b[i];
        Vector2D b1 = b[(i + 1) % 4];

        for (int j = 0; j < n; ++j)
        {
            Vector2D a2 = v[j];
            Vector2D b2 = v[(j + 1) % n];

            Vector2D p1, p2;
            real_t d = detail::distanceSegmentToSegment(a1, b1, a2, b2, &p1, &p2);

            if (d < minDist)
            {
                minDist = d;
                bestP1 = p1;
                bestP2 = p2;

                if (d == real_t(0))
                {
                    if (closestPoint1) *closestPoint1 = bestP1;
                    if (closestPoint2) *closestPoint2 = bestP2;
                    return real_t(0);
                }
            }
        }
    }

    if (closestPoint1) *closestPoint1 = bestP1;
    if (closestPoint2) *closestPoint2 = bestP2;

    return minDist;
}

real_t distance(const BBox2D& bbox1, const IPolygonalShape2D& poly2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const real_t bDist = boundaryDistance(bbox1, poly2, closestPoint1, closestPoint2);
    if (bDist == real_t{0})
        return real_t{0};

    const auto vertices = poly2.vertices();
    for (const Vector2D& vertex : vertices)
    {
        if (isPointInsideBBox(vertex, bbox1))
        {
            if (closestPoint1) *closestPoint1 = vertex;
            if (closestPoint2) *closestPoint2 = vertex;
            return real_t{0};
        }
    }

    for (const Vector2D& corner : bbox1.corners())
    {
        if (poly2.contains(corner))
        {
            if (closestPoint1) *closestPoint1 = corner;
            if (closestPoint2) *closestPoint2 = corner;
            return real_t{0};
        }
    }

    return bDist;
}

real_t distance(const BBox2D& bbox, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const Vector2D& c = circle.m_center;
    const real_t r = circle.m_radius;

    Vector2D p;
    p.x = std::clamp(c.x, bbox.m_min.x, bbox.m_max.x);
    p.y = std::clamp(c.y, bbox.m_min.y, bbox.m_max.y);

    Vector2D delta = c - p;
    real_t distSq = delta.lengthSquared();

    if (distSq == real_t(0))
    {
        // Find closest point on box boundary instead
        real_t left   = c.x - bbox.m_min.x;
        real_t right  = bbox.m_max.x - c.x;
        real_t bottom = c.y - bbox.m_min.y;
        real_t top    = bbox.m_max.y - c.y;

        // Pick smallest penetration direction
        real_t minDist = left;
        p = { bbox.m_min.x, c.y };

        if (right < minDist) { minDist = right;  p = { bbox.m_max.x, c.y }; }
        if (bottom < minDist){ minDist = bottom; p = { c.x, bbox.m_min.y }; }
        if (top < minDist)   { minDist = top;    p = { c.x, bbox.m_max.y }; }

        delta = c - p;
        real_t dist = delta.length();

        Vector2D dir = (dist > real_t(0)) ? (delta / dist) : Vector2D{1, 0};

        if (minDist <= r)
        {
            Vector2D cp = c - dir * r;

            if (closestPoint1) *closestPoint1 = cp;
            if (closestPoint2) *closestPoint2 = cp;
            return real_t(0);
        }

        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = c - dir * r;

        return minDist - r;
    }

    real_t dist = std::sqrt(distSq);
    Vector2D dir = delta / dist;

    // Overlapping or touching
    if (dist <= r)
    {
        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = p;
        return real_t(0);
    }

    // Separated
    if (closestPoint1) *closestPoint1 = p;
    if (closestPoint2) *closestPoint2 = c - dir * r;

    return dist - r;
}

template <PolygonalShape Polygon>
real_t boundaryDistance(const Polygon& poly1, const Polygon& poly2, 
    Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const auto v1 = poly1.vertices();
    const auto v2 = poly2.vertices();

    const int n1 = (int)v1.size();
    const int n2 = (int)v2.size();

    real_t minDist = std::numeric_limits<real_t>::max();

    Vector2D bestP1, bestP2;

    for (int i = 0; i < n1; ++i)
    {
        Vector2D a1 = v1[i];
        Vector2D b1 = v1[(i + 1) % n1];

        for (int j = 0; j < n2; ++j)
        {
            Vector2D a2 = v2[j];
            Vector2D b2 = v2[(j + 1) % n2];

            Vector2D p1, p2;
            real_t d = detail::distanceSegmentToSegment(a1, b1, a2, b2, &p1, &p2);

            if (d < minDist)
            {
                minDist = d;
                bestP1 = p1;
                bestP2 = p2;

                if (minDist == real_t(0))
                {
                    if (closestPoint1) *closestPoint1 = bestP1;
                    if (closestPoint2) *closestPoint2 = bestP2;
                    return real_t(0);
                }
            }
        }
    }

    if (closestPoint1) *closestPoint1 = bestP1;
    if (closestPoint2) *closestPoint2 = bestP2;

    return minDist;
}

real_t distance(const IPolygonalShape2D& poly1, const IPolygonalShape2D& poly2, 
    Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    real_t bDist = boundaryDistance(poly1, poly2, closestPoint1, closestPoint2);
    
    // Boundary Collision
    if (bDist == real_t(0)) return real_t(0);

    const auto v1 = poly1.vertices();
    const auto v2 = poly2.vertices();

    // Test for containment
    if ((!v1.empty() && poly2.contains(v1[0])) || 
        (!v2.empty() && poly1.contains(v2[0])))
    {
        // When contained, the closest witness points coincide inside the region
        Vector2D overlapPt = !v1.empty() ? v1[0] : v2[0];
        if (closestPoint1) *closestPoint1 = overlapPt;
        if (closestPoint2) *closestPoint2 = overlapPt;
        return real_t(0);
    }

    return bDist;
}

template <PolygonalShape Polygon>
real_t signedDistance(const Polygon& poly1, const Polygon& poly2, 
    Vector2D* closestPoint1 = nullptr, Vector2D* closestPoint2 = nullptr)
{
    real_t bDist = boundaryDistance(poly1, poly2, closestPoint1, closestPoint2);
    
    // Boundary Collision
    if (bDist == real_t(0)) return real_t(0);

    const auto v1 = poly1.vertices();
    const auto v2 = poly2.vertices();

    // Test for containment
    if ((!v1.empty() && poly2.contains(v1[0])) || 
        (!v2.empty() && poly1.contains(v2[0])))
    {
        return -bDist; // Negate to signal penetration
    }

    return bDist;
}

real_t boundaryDistance(const IPolygonalShape2D& poly1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const auto v = poly1.vertices();
    const int n = (int)v.size();

    const Vector2D& center = circle2.m_center;
    const real_t radius = circle2.m_radius;

    real_t minDist = std::numeric_limits<real_t>::max();

    Vector2D bestP1, bestP2;

    // --- Edge vs Circle ---
    for (int i = 0; i < n; ++i)
    {
        Vector2D a = v[i];
        Vector2D b = v[(i + 1) % n];

        Vector2D p1, p2;
        real_t d = distance(Segment2D(a, b), circle2, &p1, &p2);

        if (d < minDist)
        {
            minDist = d;
            bestP1 = p1;
            bestP2 = p2;

            if (minDist == real_t(0))
            {
                if (closestPoint1) *closestPoint1 = bestP1;
                if (closestPoint2) *closestPoint2 = bestP2;
                return real_t(0);
            }
        }
    }

    if (closestPoint1) *closestPoint1 = bestP1;
    if (closestPoint2) *closestPoint2 = bestP2;
    return minDist;
}

real_t distance(const IPolygonalShape2D& poly1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    if (poly1.contains(circle2.m_center))
    {
        if (closestPoint1) *closestPoint1 = circle2.m_center;
        if (closestPoint2) *closestPoint2 = circle2.m_center;
        return real_t{0};
    }

    return boundaryDistance(poly1, circle2, closestPoint1, closestPoint2);
}

real_t distance(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    Vector2D delta = circle2.m_center - circle1.m_center;
    real_t dist = delta.length();

    Vector2D dir = (dist > real_t{0}) ? (delta / dist) : Vector2D{1, 0};

    if (dist <= circle1.m_radius + circle2.m_radius)
    {
        Vector2D p = circle1.m_center + dir * circle1.m_radius;
        if (closestPoint1) *closestPoint1 = p;
        if (closestPoint2) *closestPoint2 = p;

        return real_t{0};
    }

    if (closestPoint1) *closestPoint1 = circle1.m_center + dir * circle1.m_radius;
    if (closestPoint2) *closestPoint2 = circle2.m_center - dir * circle2.m_radius;

    return dist - circle1.m_radius - circle2.m_radius;
}

real_t distance(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const Circle2D* circle1 = shape_cast<Circle2D>(&shape1);
    const Circle2D* circle2 = shape_cast<Circle2D>(&shape2);
    const IPolygonalShape2D* polygon1 = shape1.polygonal();
    const IPolygonalShape2D* polygon2 = shape2.polygonal();

    if (circle1 && circle2)
        return distance(*circle1, *circle2, closestPoint1, closestPoint2);

    if (polygon1 && polygon2)
        return distance(*polygon1, *polygon2, closestPoint1, closestPoint2);

    if (polygon1 && circle2)
        return distance(*polygon1, *circle2, closestPoint1, closestPoint2);

    if (circle1 && polygon2)
    {
        return distance(*polygon2, *circle1, closestPoint2, closestPoint1);
    }

    //Todo: Consider behaviour
    return std::numeric_limits<real_t>::max();
}

real_t distance(const FiniteShapeVariant2D& shape1, const FiniteShapeVariant2D& shape2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return std::visit([&](const auto& shape1, const auto& shape2)
    {
        return distance(shape1, shape2, closestPoint1, closestPoint2);
    }, shape1, shape2);
}

    
}