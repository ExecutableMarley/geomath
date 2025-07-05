/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms2D.hpp"

#include "CommonMath.hpp"
#include "../Ray2D.hpp"
#include "../Line2D.hpp"
#include "../BBox2D.hpp"
#include "../IShape2D.hpp"
#include "../Circle.hpp"
#include "../Triangle.hpp"
#include "../Rectangle.hpp"
#include "../Polygon.hpp"


namespace Utility
{

namespace Math
{

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd)
{
    Vector2D r = segmentEnd - segmentStart;
    Vector2D s = point - segmentStart;
    float cross = r.cross(s);
    if (!approximatelyZero(cross))
        return false;

    float dot = r.dot(s);
    if (dot < 0.0 || dot > r.lengthSquared())
        return false;

    return true;
}

bool isPointOnSegment(const Vector2D& point, const Line2D& line)
{
    return isPointOnSegment(point, line.m_start, line.m_end);
}

// Distance calculation algorithms

float distancePointToLine(const Vector2D& point, const Vector2D& lineStart, const Vector2D& lineEnd, Vector2D* closestPoint)
{
    const Vector2D r = lineEnd - lineStart;
    const Vector2D s = point - lineStart;
    const float rDotR = r.dot(r);
    if (rDotR == 0.0f)
    {
        if (closestPoint) *closestPoint = lineStart;
        return (point - lineStart).length();
    }
    float t = std::clamp(s.dot(r) / rDotR, 0.0f, 1.0f);
    const Vector2D projection = lineStart + r * t;
    if (closestPoint)
        *closestPoint = projection;
    return (point - projection).length();
}

float distancePointToLine(const Vector2D& point, const Line2D& line, Vector2D* closestPoint)
{
    return distancePointToLine(point, line.m_start, line.m_end, closestPoint);
}

float distanceLineToLine(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const Vector2D delta1 = s2 - s1;
    const Vector2D delta2 = k2 - k1;
    float lengthSquared1 = delta1.lengthSquared();
    float lengthSquared2 = delta2.lengthSquared();

    if (lengthSquared1 == 0 || lengthSquared2 == 0)
    {
        if (closestPoint1)
            *closestPoint1 = s1;
        if (closestPoint2)
            *closestPoint2 = k1;
        return (s1 - k1).length();
    }

    const Vector2D startDelta = s1 - k1;
    
    const float a = delta1.dot(startDelta);
    const float b = delta2.dot(startDelta);
    const float c = delta1.dot(delta2);

    const float d = lengthSquared1 * lengthSquared2 - c * c;

    float t = d != 0 ? (b * c - a * lengthSquared2) / d : 0.0f; 
    float u = d != 0 ? (lengthSquared1 * b - c * a) / d : 0.0f;

    t = clamp(t, 0.0f, 1.0f);
    u = clamp(u, 0.0f, 1.0f);

    const Vector2D closestPointOnLine1 = s1 + delta1 * t;
    const Vector2D closestPointOnLine2 = k1 + delta2 * u;

    if (closestPoint1)
        *closestPoint1 = closestPointOnLine1;
    if (closestPoint2)
        *closestPoint2 = closestPointOnLine2;
    
    return (closestPointOnLine1 - closestPointOnLine2).length();
}

float distanceLineToLine(const Line2D& line1, const Line2D& line2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distanceLineToLine(line1.m_start, line1.m_end, line2.m_start, line2.m_end, closestPoint1, closestPoint2);
}

// Intersection calculation algorithms

bool intersectRayWithBBox(const Ray2D& ray, const BBox2D& bbox, float t_min, float t_max, HitInfo2D* hitInfo)
{
    for (int i = 0; i < 2; i++)
    {
        if (approximatelyZero(ray.m_direction[i]))
        {
            if (ray.m_origin[i] < bbox.m_min[i] || ray.m_origin[i] > bbox.m_max[i])
                return false;
            continue;
        }

        float invD = 1.0f / ray.m_direction[i];
        float t0 = (bbox.m_min[i] - ray.m_origin[i]) * invD;
        float t1 = (bbox.m_max[i] - ray.m_origin[i]) * invD;

        if (invD < 0.0f)
            std::swap(t0, t1);
        
        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        if (t_max <= t_min)
            return false;
    }

    if (hitInfo)
    {
        hitInfo->t = t_min;
        hitInfo->intersectionPoint = ray.m_origin + ray.m_direction * t_min;
    }

    return true;
}

bool intersectRayWithCircle(const Ray2D& ray, const Circle& circle, float t_min, float t_max, HitInfo2D* hitInfo)
{
    const Vector2D oc = ray.m_origin - circle.m_center;
    const float a = ray.m_direction.dot(ray.m_direction);
    const float b = 2.0f * oc.dot(ray.m_direction);
    const float c = oc.dot(oc) - circle.m_radius * circle.m_radius;
    const float discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
        return false;

    const float sqrtDiscriminant = sqrt(discriminant);
    float t = (-b - sqrtDiscriminant) / (2.0f * a);
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
        //Vector2D normal = (hitInfo->intersectionPoint - circle.m_center).normalize();
    }

    return true;
}

bool intersectRayWithTriangle(const Ray2D& ray, const Triangle& triangle, float t_min, float t_max, HitInfo2D* hitInfo)
{
    const Vector2D edge1 = triangle.m_b - triangle.m_a;
    const Vector2D edge2 = triangle.m_c - triangle.m_a;

    // Calculate determinant
    const Vector2D h = { -ray.m_direction.y, ray.m_direction.x };
    const float det = edge1.dot(h);

    if (approximatelyZero(det))
        return false;

    const float invDet = 1.0f / det;

    const Vector2D s = ray.m_origin - triangle.m_a;
    const float u = s.dot(h) * invDet;

    if (u < 0.0f || u > 1.0f)
        return false;

    const Vector2D q = { -s.y, s.x };
    const float v = ray.m_direction.dot(q) * invDet;

    if (v < 0.0f || u + v > 1.0f)
        return false;

    const float t = edge2.dot(q) * invDet;

    if (t < t_min || t > t_max)
        return false;

    if (hitInfo) 
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = ray.pointAt(t);
    }

    return true;
}

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, float t_min, float t_max, HitInfo2D* hitInfo)
{
    Vector2D v1 = ray.m_origin - p1;
    Vector2D v2 = p2 - p1;
    Vector2D v3(-ray.m_direction.y, ray.m_direction.x);

    float dot = v2.dot(v3);
    
    if (approximatelyZero(dot))
        return false; // Parallel

    float t1 = v2.cross(v1) / dot;
    float t2 = v1.dot(v3) / dot;

    if (t1 < 0 || t2 < 0 || t2 > 1)
        return false;

    if (t1 >= t_min && t1 <= t_max)
    {
        if (hitInfo)
        {
            hitInfo->t = t1;
            hitInfo->intersectionPoint = ray.m_origin + ray.m_direction * t1;
        }
        return true;
    }

    return false;
}

bool intersectRayWithRectangleOptimized(const Ray2D& ray, const Rectangle& rectangle, float t_min, float t_max)
{
    for (int i = 0; i < 4; i++)
    {
        const Vector2D& p1 = rectangle.vertexAt(i);
        const Vector2D& p2 = rectangle.vertexAt((i + 1) % 4);

        if (intersectRayWithSegment(ray, p1, p2, t_min, t_max, NULL))
        {
            return true;
        }
    }
    return false;
}

bool intersectRayWithRectangle(const Ray2D& ray, const Rectangle& rectangle, float t_min, float t_max, HitInfo2D* hitInfo)
{
    if (hitInfo == nullptr)
        return intersectRayWithRectangleOptimized(ray, rectangle, t_min, t_max);

    bool hit = false;
    float t = t_max;

    for (int i = 0; i < 4; i++)
    {
        const Vector2D& p1 = rectangle.vertexAt(i);
        const Vector2D& p2 = rectangle.vertexAt((i + 1) % 4);

        if (intersectRayWithSegment(ray, p1, p2, t_min, t, hitInfo))
        {
            hit = true;
            t = hitInfo->t;
        }
    }
    return hit;
}

bool intersectRayWithPolygonOptimized(const Ray2D& ray, const Polygon& polygon, float t_min, float t_max)
{
    for (int i = 0; i < polygon.vertexCount(); i++)
    {
        const Vector2D& p1 = polygon.vertexAt(i);
        const Vector2D& p2 = polygon.vertexAt((i + 1) % polygon.vertexCount());

        if (intersectRayWithSegment(ray, p1, p2, t_min, t_max, NULL))
        {
            return true;
        }
    }
    return false;
}

bool intersectRayWithPolygon(const Ray2D& ray, const Polygon& polygon, float t_min, float t_max, HitInfo2D* hitInfo)
{
    if (hitInfo == nullptr)
        return intersectRayWithPolygonOptimized(ray, polygon, t_min, t_max);

    bool hit = false;
    float t = t_max;

    for (int i = 0; i < polygon.vertexCount(); i++)
    {
        const Vector2D& p1 = polygon.vertexAt(i);
        const Vector2D& p2 = polygon.vertexAt((i + 1) % polygon.vertexCount());

        if (intersectRayWithSegment(ray, p1, p2, t_min, t_max, hitInfo))
        {
            hit = true;
            t = hitInfo->t;
        }
    }
    return hit;
}

bool intersectRayWithShape(const Ray2D& ray, const IShape2D& shape, float t_min, float t_max, HitInfo2D* hitInfo)
{
    switch(shape.type())
    {
        case ShapeType2D::SHAPE2D_CIRCLE:
            return intersectRayWithCircle(ray, dynamic_cast<const Circle&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_RECTANGLE:
            return intersectRayWithRectangle(ray, dynamic_cast<const Rectangle&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_TRIANGLE:
            return intersectRayWithTriangle(ray, dynamic_cast<const Triangle&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_POLYGON:
            return intersectRayWithPolygon(ray, dynamic_cast<const Polygon&>(shape), t_min, t_max, hitInfo);
        default:
            //Todo: we might just want to throw and error here
            return false;
    }
}

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    Vector2D r = p2 - p1;
    Vector2D s = q2 - q1;
    Vector2D qp = q1 - p1;

    float rxs = r.cross(s);
    float qpxr = qp.cross(r);

    if (approximatelyZero(rxs))
        return false; // Parallel/collinear

    float t = qp.cross(s) / rxs;
    float u = qp.cross(r) / rxs;

    if (t <= 0.0f || t >= 1.0f || u <= 0.0f || u >= 1.0f)
        return false;

    if (hitInfo)
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = p1 + r * t;
    }

    return true;
}

bool intersectSegmentWithSegmentStrict(const Line2D& line1, const Line2D& line2, HitInfo2D* hitInfo)
{
    return intersectSegmentWithSegmentStrict(line1.m_start, line1.m_end, line2.m_start, line2.m_end, hitInfo);
}

bool intersectSegmentWithPolygon(const Line2D& line, const Polygon& polygon, HitInfo2D* hitInfo)
{
    return(intersectRayWithPolygon(Ray2D(line.m_start, line.direction()), polygon, 0.0f, line.length(), hitInfo));
}

} // namespace Math

} // namespace Utility