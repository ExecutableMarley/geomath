/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

namespace Arns::geomath
{
    
struct HitInfo2D
{
    real_t t;
    Vector2D intersectionPoint;
    Vector2D normal;
    //IShape2D* shape;
};

namespace detail
{

    bool intersectRayWithSegment(const Ray2D &ray, const Vector2D &p1, const Vector2D &p2,
                                 real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    bool rayBoundsInterval(const Ray2D &ray, const BBox2D &bbox,
                           real_t t_min, real_t t_max, real_t &t_enter, real_t &t_exit);

    bool intersectRayWithPolygon(const Ray2D &ray, std::span<const Vector2D> vertices,
                                 real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    bool intersectSegmentWithSegment(const Vector2D &p1, const Vector2D &p2, const Vector2D &q1, const Vector2D &q2, HitInfo2D *hitInfo = nullptr);

    bool intersectSegmentWithSegmentStrict(const Vector2D &p1, const Vector2D &p2, const Vector2D &q1, const Vector2D &q2, HitInfo2D *hitInfo = nullptr);

    //    --- Rays ---

    bool intersectRayImpl(const Ray2D &ray, const Segment2D &segment, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    bool intersectRayImpl(const Ray2D &ray, const BBox2D &bbox, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    template <PolygonalShape Polygon>
    bool intersectRayImpl(const Ray2D &ray, const Polygon &polygon, real_t t_min = 0,
                          real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr)
    {
        return detail::intersectRayWithPolygon(ray, polygon.vertices(), t_min, t_max, hitInfo);
    }

    bool intersectRayImpl(const Ray2D &ray, const Circle2D &circle, real_t t_min = 0, real_t t_max = std::numeric_limits<real_t>::max(), HitInfo2D *hitInfo = nullptr);

    template <typename Object>
    concept HasRayIntersectionImpl =
        requires(
            const Ray2D &ray,
            const Object &object,
            real_t tMin,
            real_t tMax,
            HitInfo2D *hitInfo) {
            {
                intersectRayImpl(
                    ray, object, tMin, tMax, hitInfo)
            } -> std::convertible_to<bool>;
        };

    bool intersectImpl(const BBox2D &bbox1, const BBox2D &bbox2);
    bool intersectImpl(const BBox2D &bbox, const Circle2D &circle);
    bool intersectImpl(const Circle2D &circle1, const Circle2D &circle2);

    bool intersectBBoxWithPolygon(const BBox2D &bbox, std::span<const Vector2D> vertices, bool convex);
    bool intersectPolygonWithPolygon(std::span<const Vector2D> vertices1, bool convex1,
                                     std::span<const Vector2D> vertices2, bool convex2);
    bool intersectPolygonWithCircle(std::span<const Vector2D> vertices, bool convex, const Circle2D &circle);

    template <PolygonalShape Polygon>
    bool intersectImpl(const BBox2D &bbox, const Polygon &polygon)
    {
        return intersectBBoxWithPolygon(bbox, polygon.vertices(), polygon.isConvex());
    }

    template <PolygonalShape Polygon1, PolygonalShape Polygon2>
    bool intersectImpl(const Polygon1 &polygon1, const Polygon2 &polygon2)
    {
        return intersectPolygonWithPolygon(polygon1.vertices(), polygon1.isConvex(),
                                           polygon2.vertices(), polygon2.isConvex());
    }

    template <PolygonalShape Polygon>
    bool intersectImpl(const Polygon &polygon, const Circle2D &circle)
    {
        return intersectPolygonWithCircle(polygon.vertices(), polygon.isConvex(), circle);
    }
}

//=======================================
// Public Intersection dispatch functions
//=======================================

//Todo: Consider renaming ray,object and segment,object overloads due to different calling conventions

//    --- Ray intersection function ---

template <typename B> requires (detail::HasRayIntersectionImpl<B>)
bool raycast(
    const Ray2D& first, const B& second, real_t tMin = real_t{0},
        real_t tMax = std::numeric_limits<real_t>::max(),
        HitInfo2D* hitInfo = nullptr)
{
    return detail::intersectRayImpl(first, second, tMin, tMax, hitInfo);
}

//    --- Segment intersection function ---

template <typename Object> requires detail::HasRayIntersectionImpl<Object>
bool raycast(const Segment2D& segment, const Object& object, HitInfo2D* hitInfo = nullptr)
{
    const Ray2D ray(segment.m_start, segment.deltaVector());
    return raycast(ray, object, real_t{0}, real_t{1}, hitInfo);
}

//    --- Shape-Shape intersection function ---

template <typename A, typename B>
    requires (requires(const A& a, const B& b) { detail::intersectImpl(a, b); } ||
              requires(const A& a, const B& b) { detail::intersectImpl(b, a); })
bool intersect(const A& a, const B& b)
{
    if constexpr (requires { detail::intersectImpl(a, b); })
        return detail::intersectImpl(a, b);
    else
        return detail::intersectImpl(b, a);
}
    
} // namespace Arns::geomath