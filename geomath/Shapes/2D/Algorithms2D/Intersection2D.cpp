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

namespace Arns::geomath
{


// Intersection calculation algorithms

//    --- Rays ---


namespace detail
{

// Internal Helper
bool intersectParams(const Vector2D& p, const Vector2D& r, const Vector2D& q, const Vector2D& s, 
    real_t& t, real_t& u)
{
    const real_t rxs = r.cross(s);
    const Vector2D qp = q - p;

    if (approximatelyZero(rxs))
        return false; // parallel

    t = qp.cross(s) / rxs;
    u = qp.cross(r) / rxs;
    return true;
}

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    real_t t, u;
    if (!intersectParams(ray.m_origin, ray.m_direction, p1, p2 - p1, t, u))
        return false;

    if (t < t_min || t > t_max || u < 0.f || u > 1.f)
        return false;

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = ray.m_origin + ray.m_direction * t;
        hitInfo->normal = (p2 - p1).perpendicular().normalize();
    
        if (hitInfo->normal.dot(ray.m_direction) > 0) hitInfo->normal = -hitInfo->normal;
    }

    return true;
}

bool rayBoundsInterval(const Ray2D& ray, const BBox2D& bbox, real_t t_min, real_t t_max, real_t& t_enter, real_t& t_exit)
{
    t_enter = t_min;
    t_exit = t_max;

    for (int i = 0; i < 2; i++)
    {
        if (approximatelyZero(ray.m_direction[i]))
        {
            if (ray.m_origin[i] < bbox.m_min[i] || ray.m_origin[i] > bbox.m_max[i])
                return false;
            continue;
        }

        real_t invD = real_t{1} / ray.m_direction[i];
        real_t t0 = (bbox.m_min[i] - ray.m_origin[i]) * invD;
        real_t t1 = (bbox.m_max[i] - ray.m_origin[i]) * invD;

        if (invD < real_t{0})
            std::swap(t0, t1);
        
        t_enter = t0 > t_enter ? t0 : t_enter;
        t_exit = t1 < t_exit ? t1 : t_exit;

        if (t_exit <= t_enter)
            return false;
    }

    return true;
}

// --- Misc Helper ---

inline void projectPointSetOntoAxis(std::span<const Vector2D> pts, const Vector2D& axis,
                                    real_t& outMin, real_t& outMax)
{
    outMin = std::numeric_limits<real_t>::max();
    outMax = std::numeric_limits<real_t>::lowest();

    for (const auto& p : pts)
    {
        real_t v = p.dot(axis);
        if (v < outMin) outMin = v;
        if (v > outMax) outMax = v;
    }
}

inline void projectBBoxOntoAxis(const BBox2D& bbox, const Vector2D& axis,
                                real_t& outMin, real_t& outMax)
{
    Vector2D center  = (bbox.m_min + bbox.m_max) * 0.5f;
    Vector2D extents = (bbox.m_max - bbox.m_min) * 0.5f;

    // The "radius" of the box projection
    real_t radius = std::abs(extents.x * axis.x) + std::abs(extents.y * axis.y);
    real_t centerProj = center.dot(axis);

    outMin = centerProj - radius;
    outMax = centerProj + radius;
    //return { centerProj - radius, centerProj + radius }
}

// Generic SAT-based intersection for a convex polygon with a BBox
inline bool intersectBBoxWithConvexShape(const BBox2D& bbox, std::span<const Vector2D> verts)
{
    if (verts.size() < 3) return false; // degenerate

    real_t vMinX = verts[0].x, vMaxX = verts[0].x;
    real_t vMinY = verts[0].y, vMaxY = verts[0].y;

    for (size_t i = 1; i < verts.size(); ++i)
    {
        if (verts[i].x < vMinX) vMinX = verts[i].x;
        if (verts[i].x > vMaxX) vMaxX = verts[i].x;
        if (verts[i].y < vMinY) vMinY = verts[i].y;
        if (verts[i].y > vMaxY) vMaxY = verts[i].y;
    }

    if (vMaxX < bbox.m_min.x || vMinX > bbox.m_max.x) return false;
    if (vMaxY < bbox.m_min.y || vMinY > bbox.m_max.y) return false;

    real_t bmin, bmax, tmin, tmax;

    for (size_t i = 0; i < verts.size(); ++i)
    {
        // Get the edge between current and next vertex (wrapping around)
        const Vector2D& p1 = verts[i];
        const Vector2D& p2 = verts[(i + 1) % verts.size()];
        
        Vector2D edge = p2 - p1;
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges (nearly identical vertices)
        if (axis.isZero()) continue;

        projectBBoxOntoAxis    (bbox , axis, bmin, bmax);
        projectPointSetOntoAxis(verts, axis, tmin, tmax);

        if (!intervalsOverlap(bmin, bmax, tmin, tmax))
            return false;
    }

    return true; // No separating axis found
}

}

bool intersect(const Ray2D& ray, const Segment2D& segment, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    return detail::intersectRayWithSegment(ray, segment.m_start, segment.m_end, t_min, t_max, hitInfo);
}

bool intersect(const Ray2D& ray, const BBox2D& bbox, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    real_t t_enter;
    real_t t_exit;
    if (!detail::rayBoundsInterval(ray, bbox, t_min, t_max, t_enter, t_exit))
        return false;

    if (hitInfo)
    {
        real_t hit_t = (t_enter > t_min) ? t_enter : t_exit;
        hitInfo->t = hit_t;
        hitInfo->intersectionPoint = ray.m_origin + ray.m_direction * hit_t;
        hitInfo->normal = bbox.normalAt(hitInfo->intersectionPoint);
    }

    return true;
}

bool intersectRayWithPolygonOptimized(const Ray2D& ray, const IPolygonalShape2D& polygon, real_t t_min, real_t t_max)
{
    for (int i = 0; i < polygon.vertexCount(); i++)
    {
        const Vector2D& p1 = polygon[i];
        const Vector2D& p2 = polygon[(i + 1) % polygon.vertexCount()];

        if (detail::intersectRayWithSegment(ray, p1, p2, t_min, t_max, NULL))
        {
            return true;
        }
    }
    return false;
}

bool intersect(const Ray2D& ray, const IPolygonalShape2D& polygon, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    if (hitInfo == nullptr)
        return intersectRayWithPolygonOptimized(ray, polygon, t_min, t_max);

    bool hit = false;
    real_t t = t_max;

    for (int i = 0; i < polygon.vertexCount(); i++)
    {
        const Vector2D& p1 = polygon[i];
        const Vector2D& p2 = polygon[(i + 1) % polygon.vertexCount()];

        if (detail::intersectRayWithSegment(ray, p1, p2, t_min, t, hitInfo))
        {
            hit = true;
            t = hitInfo->t;
        }
    }
    return hit;
}

bool intersect(const Ray2D& ray, const Circle2D& circle, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    const Vector2D oc = ray.m_origin - circle.m_center;
    const real_t a = ray.m_direction.dot(ray.m_direction);
    const real_t b = 2.0f * oc.dot(ray.m_direction);
    const real_t c = oc.dot(oc) - circle.m_radius * circle.m_radius;
    const real_t discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
        return false;

    const real_t sqrtDiscriminant = sqrt(discriminant);
    real_t t = (-b - sqrtDiscriminant) / (2.0f * a);
    if (t < t_min || t > t_max)
    {
        t = (-b + sqrtDiscriminant) / (2.0f * a);
        if (t < t_min || t > t_max)
            return false;
    }

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = ray.pointAt(t);
        hitInfo->normal = (hitInfo->intersectionPoint - circle.m_center).normalize();
    }

    return true;
}

bool intersect(const Ray2D& ray, const IBaseShape2D& shape, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    switch(shape.type())
    {
        case ShapeType2D::SHAPE2D_TRIANGLE:
        case ShapeType2D::SHAPE2D_RECTANGLE:
        case ShapeType2D::SHAPE2D_POLYGON:
            return intersect(ray, dynamic_cast<const IPolygonalShape2D&>(shape), t_min, t_max, hitInfo);
         case ShapeType2D::SHAPE2D_CIRCLE:
            return intersect(ray, dynamic_cast<const Circle2D&>(shape), t_min, t_max, hitInfo);
        default:
            //Todo: we might just want to throw and error here
            assert("Unreachable Code reached!");
            return false;
    }
}

bool intersectSegmentWithSegment(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    return detail::intersectRayWithSegment(Ray2D(p1, p2 - p1), q1, q2, 0.f, 1.f, hitInfo);
}

bool intersect(const Segment2D& segment1, const Segment2D& segment2, HitInfo2D* hitInfo)
{
    return detail::intersectRayWithSegment(Ray2D(segment1.m_start, segment1.deltaVector()), segment2.m_start, segment2.m_end, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    real_t t, u;
    if (!detail::intersectParams(p1, p2 - p1, q1, q2 - q1, t, u))
        return false;

    if (!approximatelyGreater(t, real_t{0}) || !approximatelyLess(t, real_t{1}) ||
        !approximatelyGreater(u, real_t{0}) || !approximatelyLess(u, real_t{1}))
        return false;

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = p1 + (p2 - p1) * t;
        hitInfo->normal = (q2 - q1).perpendicular().normalize();

        if (hitInfo->normal.dot(p2 - p1) > 0) hitInfo->normal = -hitInfo->normal;
    }

    return true;
}

bool intersect(const Segment2D& segment, const BBox2D& bbox, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(segment.m_start, segment.deltaVector()), bbox, 0.f, 1.f, hitInfo);
}

bool intersect(const Segment2D& segment, const Triangle2D& triangle, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(segment.m_start, segment.deltaVector()), triangle, 0.f, 1.f, hitInfo);
}

bool intersect(const Segment2D& segment, const IPolygonalShape2D& polygon, HitInfo2D* hitInfo)
{
    return(intersect(Ray2D(segment.m_start, segment.direction()), polygon, 0.0f, segment.length(), hitInfo));
}

bool intersect(const Segment2D& segment, const Circle2D& circle, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(segment.m_start, segment.deltaVector()), circle, 0.f, 1.f, hitInfo);
}

bool intersect(const Segment2D& segment, const IBaseShape2D& shape, HitInfo2D* hitInfo)
{
    return intersect(Ray2D(segment.m_start, segment.direction()), shape, 0.0f, segment.length(), hitInfo);
}


// --- Bounding Box ---

bool intersect(const BBox2D& b1, const BBox2D& b2)
{
    //Separating Axis Theorem (SAT)
#if 0
    bool separated_x = (b1.m_min.x > b2.m_max.x) || (b2.m_min.x > b1.m_max.x);
    bool separated_y = (b1.m_min.y > b2.m_max.y) || (b2.m_min.y > b1.m_max.y);

    if (separated_x || separated_y)
        return false;

    return true;
#else
    return intervalsOverlap(b1.m_min.x, b1.m_max.x, b2.m_min.x, b2.m_max.x) &&
           intervalsOverlap(b1.m_min.y, b1.m_max.y, b2.m_min.y, b2.m_max.y);
#endif
}


bool intersect(const BBox2D& bbox, const IPolygonalShape2D& polygon)
{
    if (isConvexPolygonal(polygon.type()))
    {
        const ConvexPolygon2D* convex = static_cast<const ConvexPolygon2D*>(&polygon);
        return detail::intersectBBoxWithConvexShape(bbox, polygon.vertices());
    }
    
    Vector2D box[4] =
    {
        bbox.m_min,
        Vector2D(bbox.m_max.x, bbox.m_min.y),
        bbox.m_max,
        Vector2D(bbox.m_min.x, bbox.m_max.y)
    };

    for (const Vector2D& p : box)
        if (polygon.contains(p))
            return true;

    for (const Vector2D& p : polygon.vertices())
        if (bbox.contains(p))
            return true;

    Vector2D boxEdges[4][2] =
    {
        { box[0], box[1] },
        { box[1], box[2] },
        { box[2], box[3] },
        { box[3], box[0] }
    };

    const size_t n = polygon.vertices().size();

    for (size_t i = 0; i < n; ++i)
    {
        Vector2D a = polygon.vertices()[i];
        Vector2D b = polygon.vertices()[(i + 1) % n];

        for (int j = 0; j < 4; ++j)
            if (intersectSegmentWithSegment(a, b, boxEdges[j][0], boxEdges[j][1]))
                return true;
    }

    return false;
}

bool intersect(const BBox2D& bbox, const Circle2D& circle)
{
    real_t distSq = bbox.minDistanceSquared(circle.m_center);
    real_t rSq = circle.m_radius * circle.m_radius;
    return distSq <= rSq;
}

// --- Polygons --- 

bool intersectConvexPolygonWithConvexPolygon(const std::span<const Vector2D> v1, const std::span<const Vector2D> v2)
{
    const int n1 = (int)v1.size();
    const int n2 = (int)v2.size();

    if (n1 < 3 || n2 < 3)
        return false;

    real_t min1, max1, min2, max2;

    for (int i = 0; i < n1; ++i)
    {
        const Vector2D& a = v1[i];
        const Vector2D& b = v1[(i + 1) % n1];
        Vector2D edge = b - a;

        // Normal
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        detail::projectPointSetOntoAxis(v1, axis, min1, max1);
        detail::projectPointSetOntoAxis(v2, axis, min2, max2);

        if (!intervalsOverlap(min1, max1, min2, max2))
            return false; // separating axis found
    }

    for (int i = 0; i < n2; ++i)
    {
        const Vector2D& a = v2[i];
        const Vector2D& b = v2[(i + 1) % n2];
        Vector2D edge = b - a;

        // Normal
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        detail::projectPointSetOntoAxis(v1, axis, min1, max1);
        detail::projectPointSetOntoAxis(v2, axis, min2, max2);

        if (!intervalsOverlap(min1, max1, min2, max2))
            return false;
    }

    // No separating axis -> intersection
    return true;
}

bool intersectPolygonWithPolygon(const std::span<const Vector2D> v1, const std::span<const Vector2D> v2)
{
    const int n1 = (int)v1.size();
    const int n2 = (int)v2.size();

    if (n1 < 3 || n2 < 3)
        return false;

    for (const auto& v : v1)
    {
        if (isPointInsidePolygon(v, v2))
            return true;
    }

    for (const auto& v : v2)
    {
        if (isPointInsidePolygon(v, v1))
            return true;
    }

    // --- 2. Edge intersection tests ---
    for (size_t i = 0; i < n1; ++i)
    {
        Vector2D a1 = v1[i];
        Vector2D a2 = v1[(i + 1) % n1];

        for (size_t j = 0; j < n2; ++j)
        {
            Vector2D b1 = v2[j];
            Vector2D b2 = v2[(j + 1) % n2];

            if (intersectSegmentWithSegment(a1, a2, b1, b2))
                return true;
        }
    }

   return false;
}

bool intersectPolygonWithCircle(const std::span<const Vector2D> poly, const Circle2D& circle)
{
    if (isPointInsidePolygon(circle.m_center, poly))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    real_t r = circle.m_radius;

    for (int i = 0; i < poly.size(); ++i)
    {
        const Vector2D& a = poly[i];
        const Vector2D& b = poly[(i + 1) % poly.size()];

        if (detail::distancePointToSegment(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

// --- Polygonal ---

bool intersect(const IPolygonalShape2D& polygon1, const IPolygonalShape2D& polygon2)
{
    if (isConvexPolygonal(polygon1.type()) && isConvexPolygonal(polygon2.type()))
        return intersectConvexPolygonWithConvexPolygon(polygon1.vertices(), polygon2.vertices());

    return intersectPolygonWithPolygon(polygon1.vertices(), polygon2.vertices());
}

bool intersect(const IPolygonalShape2D& polygon, const Circle2D& circle)
{
    if (polygon.contains(circle.centroid()))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    size_t n = polygon.vertexCount();
    real_t r = circle.m_radius;

    for (int i = 0; i < polygon.vertexCount(); ++i)
    {
        const Vector2D& a = polygon[i];
        const Vector2D& b = polygon[(i + 1) % n];

        if (detail::distancePointToSegment(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

// --- Circles ---

bool intersect(const Circle2D& circle1, const Circle2D& circle2)
{
    return circle1.intersects(circle2);
}

// --- Generic ---

bool intersect(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2)
{
    const Circle2D* circle1 = shape_cast<Circle2D>(&shape1);
    const Circle2D* circle2 = shape_cast<Circle2D>(&shape2);
    const IPolygonalShape2D* polygon1 = shape1.polygonal();
    const IPolygonalShape2D* polygon2 = shape2.polygonal();

    if (circle1 && circle2)
        return intersect(*circle1, *circle2);

    if (polygon1 && polygon2)
        return intersect(*polygon1, *polygon2);

    if (polygon1 && circle2)
        return intersect(*polygon1, *circle2);

    if (circle1 && polygon2)
        return intersect(*polygon2, *circle1);

    //Todo: Consider behaviour
    return false;
}
    
}