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
#include "../Circle2D.hpp"
#include "../Triangle2D.hpp"
#include "../Rectangle2D.hpp"
#include "../Polygon2D.hpp"
#include "../ConvexPolygon2D.hpp"

namespace Arns
{

namespace Math
{

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd)
{
    Vector2D r = segmentEnd - segmentStart;
    Vector2D s = point - segmentStart;
    real_t cross = r.cross(s);
    if (!approximatelyZero(cross))
        return false;

    real_t dot = r.dot(s);
    if (dot < 0.0 || dot > r.lengthSquared())
        return false;

    return true;
}

bool isPointOnSegment(const Vector2D& point, const Line2D& line)
{
    return isPointOnSegment(point, line.m_start, line.m_end);
}

bool isSegmentOnSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2)
{
    return isPointOnSegment(s1, k1, k2) && isPointOnSegment(s2, k1, k2);
}

bool isSegmentOnSegment(const Segment2D& segment1, const Segment2D& segment2)
{
    return isSegmentOnSegment(segment1.m_start, segment1.m_end, segment2.m_start, segment2.m_end);
}

// Point inside shape

bool isPointInsideBBox(const Vector2D& point, const BBox2D& bbox)
{
    return bbox.contains(point);
}

bool isPointInsideTriangle(const Vector2D& point, const Triangle2D& triangle)
{
    return triangle.contains(point);
}

bool isPointInsideRectangle(const Vector2D& point, const Rectangle2D& rectangle)
{
    return rectangle.contains(point);
}

//Todo: Consider the O(log n) method
bool isPointInsideConvexPolygon(const Vector2D& point, const std::vector<Vector2D>& vertices)
{
    const size_t n = vertices.size();
    if (n < 3) return false;

    bool hasPositive = false;
    bool hasNegative = false;

    for (size_t i = 0; i < n; ++i)
    {
        const auto& A = vertices[i];
        const auto& B = vertices[(i + 1) % n];

        real_t val = orient2D(A, B, point);

        if (approximatelyGreaterAbs(val, real_t{0}))   hasPositive = true;
        else if (approximatelyLessAbs(val, real_t{0})) hasNegative = true;

        if (hasPositive && hasNegative) return false;
    }

    return true; 
}

bool isPointInsideConvexPolygon(const Vector2D& point, const ConvexPolygon2D& polygon)
{
    return isPointInsideConvexPolygon(point, polygon.m_vertices);
}

bool isPointInsidePolygon(const Vector2D& point, const Polygon2D& polygon)
{
    return polygon.contains(point);
}

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle)
{
    return circle.contains(point);
}


bool isSegmentInsideBBox(const Segment2D& segment, const BBox2D& bbox)
{
    return isPointInsideBBox(segment.m_start, bbox) && isPointInsideBBox(segment.m_end, bbox);
}

bool isSegmentInsideTriangle(const Segment2D& segment, const Triangle2D& triangle)
{
    return isPointInsideTriangle(segment.m_start, triangle) && isPointInsideTriangle(segment.m_end, triangle);
}

bool isSegmentInsideRectangle(const Segment2D& segment, const Rectangle2D& rectangle)
{
    return isPointInsideRectangle(segment.m_start, rectangle) && isPointInsideRectangle(segment.m_end, rectangle);
}

bool isSegmentInsideConvexPolygon(const Segment2D& segment, const ConvexPolygon2D& polygon)
{
    return isPointInsideConvexPolygon(segment.m_start, polygon) && isPointInsideConvexPolygon(segment.m_end, polygon);
}

bool isSegmentInsidePolygon(const Segment2D& segment, const Polygon2D polygon)
{
    if (isPointInsidePolygon(segment.m_start, polygon) && isPointInsidePolygon(segment.m_end, polygon))
    {
        if (intersectSegmentWithPolygon(segment, polygon))
            return false;
        return true;
    }
    return false;
}

bool isSegmentInsideCircle(const Segment2D& segment, const Circle2D& circle)
{
    return isPointInsideCircle(segment.m_start, circle) && isPointInsideCircle(segment.m_end, circle);
}


//Todo: Unfinished stuff:

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2)
{
    return bbox2.contains(bbox1);
}

bool isConvexPolygonInsideConvexPolygon(const std::vector<Vector2D>& vertices1, const std::vector<Vector2D>& vertices2)
{
    for (const auto point : vertices1)
    {
        if (!isPointInsideConvexPolygon(point, vertices2))
        {
            return false;
        }
    }
    return true;
}

bool isConvexPolygonInsideConvexPolygon(const ConvexPolygon2D& polygon1, const ConvexPolygon2D& polygon2)
{
    for (const auto point : polygon1)
    {
        if (!isPointInsideConvexPolygon(point, polygon2))
        {
            return false;
        }
    }
    return true;
}

bool isShapeInsideShape(const Triangle2D& shape1, const Triangle2D& shape2)
{
    return isConvexPolygonInsideConvexPolygon(shape1.getVertices(), shape2.getVertices());
}

bool isPolygonInsideConvexPolygon(const Polygon2D& polygon1, const ConvexPolygon2D& polygon2)
{
    for (const auto point : polygon1)
    {
        if (!isPointInsideConvexPolygon(point, polygon2))
        {
            return false;
        }
    }
    return true;
}

bool edgeIntersectionPolygonWithPolygon(const std::vector<Vector2D>& vertices1, const std::vector<Vector2D>& vertices2, HitInfo2D* hitInfo = nullptr)
{
    const size_t countA = vertices1.size();
    const size_t countB = vertices2.size();

    if (countA < 3 || countB < 3)
        return false;

    for (size_t i = 0; i < countA; ++i)
    {
        const Vector2D& a1 = vertices1[i];
        const Vector2D& a2 = vertices1[(i + 1) % countA];

        for (size_t j = 0; j < countB; ++j)
        {
            const Vector2D& b1 = vertices2[j];
            const Vector2D& b2 = vertices2[(j + 1) % countB];

            if (intersectSegmentWithSegment(a1, a2, b1, b2, hitInfo))
                return true;
        }
    }
    return false;
}

bool isPolygonInsidePolygon(const std::vector<Vector2D>& vertices1, const std::vector<Vector2D>& vertices2)
{
    for (const auto& point : vertices1)
    {
        if (!isPointInsidePolygon(point, vertices2))
        {
            return false;
        }
    }

    if (edgeIntersectionPolygonWithPolygon(vertices1, vertices2))
        return false;

    return true;
}

bool isPolygonInsidePolygon(const Polygon2D& polygon1, const Polygon2D& polygon2)
{
    return isPolygonInsidePolygon(polygon1.getVertices(), polygon2.getVertices());
}

// -----------------



// Distance calculation algorithms

real_t distancePointToLine(const Vector2D& point, const Vector2D& lineStart, const Vector2D& lineEnd, Vector2D* closestPoint)
{
    const Vector2D r = lineEnd - lineStart;
    const Vector2D s = point - lineStart;
    const real_t rDotR = r.dot(r);
    if (rDotR == 0.0f)
    {
        if (closestPoint) *closestPoint = lineStart;
        return (point - lineStart).length();
    }
    real_t t = std::clamp(s.dot(r) / rDotR, 0.0f, 1.0f);
    const Vector2D projection = lineStart + r * t;
    if (closestPoint)
        *closestPoint = projection;
    return (point - projection).length();
}

real_t distancePointToLine(const Vector2D& point, const Line2D& line, Vector2D* closestPoint)
{
    return distancePointToLine(point, line.m_start, line.m_end, closestPoint);
}

real_t distanceLineToLine(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    const Vector2D delta1 = s2 - s1;
    const Vector2D delta2 = k2 - k1;
    real_t lengthSquared1 = delta1.lengthSquared();
    real_t lengthSquared2 = delta2.lengthSquared();

    if (lengthSquared1 == 0 || lengthSquared2 == 0)
    {
        if (closestPoint1)
            *closestPoint1 = s1;
        if (closestPoint2)
            *closestPoint2 = k1;
        return (s1 - k1).length();
    }

    const Vector2D startDelta = s1 - k1;
    
    const real_t a = delta1.dot(startDelta);
    const real_t b = delta2.dot(startDelta);
    const real_t c = delta1.dot(delta2);

    const real_t d = lengthSquared1 * lengthSquared2 - c * c;

    real_t t = d != 0 ? (b * c - a * lengthSquared2) / d : 0.0f; 
    real_t u = d != 0 ? (lengthSquared1 * b - c * a) / d : 0.0f;

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

real_t distanceLineToLine(const Line2D& line1, const Line2D& line2, Vector2D* closestPoint1, Vector2D* closestPoint2)
{
    return distanceLineToLine(line1.m_start, line1.m_end, line2.m_start, line2.m_end, closestPoint1, closestPoint2);
}

//Todo: Distance to shapes

real_t distancePointToBBox(const Vector2D& point, const BBox2D& bbox, Vector2D* closestPoint)
{
    Vector2D c = bbox.closestPoint(point);
    if (closestPoint) *closestPoint = c;
    return point.distance(c);
}

//Todo: Consider renaming slightly

//Helper
real_t distancePointToConvexVertices(const Vector2D& point, const std::vector<Vector2D>& polygon, Vector2D* closestPoint)
{
    real_t minDist = std::numeric_limits<real_t>::max();
    Vector2D best;

    for (size_t i = 0; i < polygon.size(); ++i)
    {
        const Vector2D& a = polygon[i];
        const Vector2D& b = polygon[(i + 1) % polygon.size()];

        Vector2D q;
        real_t dist = distancePointToLine(point, a, b, &q);

        if (dist < minDist)
        {
            minDist = dist;
            best = q;
        }
    }

    if (closestPoint) *closestPoint = best;
    return std::sqrt(minDist);
}

//Consider unrolling loops for better performance. For now this is fine

real_t distancePointToTriangle(const Vector2D& point, const Triangle2D& triangle, Vector2D* closestPoint)
{
    if (triangle.contains(point))
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }

    return distancePointToConvexVertices(point, triangle.getVertices(), closestPoint);
}

real_t distancePointToRectangle(const Vector2D& point, const Rectangle2D& rectangle, Vector2D* closestPoint)
{
    if (rectangle.contains(point))
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }

    return distancePointToConvexVertices(point, rectangle.getVertices(), closestPoint);
}

real_t distancePointToPolygon(const Vector2D& point, const ConvexPolygon2D& polygon, Vector2D* closestPoint)
{
    if (polygon.contains(point))
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }

    return distancePointToConvexVertices(point, polygon.getVertices(), closestPoint);
}

real_t distancePointToPolygon(const Vector2D& point, const Polygon2D& polygon, Vector2D* closestPoint)
{
    if (polygon.contains(point))
    {
        if (closestPoint) *closestPoint = point;
        return real_t{0};
    }

    return distancePointToConvexVertices(point, polygon.getVertices(), closestPoint);
}

real_t distancePointToCircle(const Vector2D& point, const Circle2D& circle, Vector2D* closestPoint)
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

//Todo: Implement
real_t distanceSegmentToBBox(const Vector2D& segmentStart, Vector2D& segmentEnd, const BBox2D& bbox, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceSegmentToBBox(const Line2D& segment, const BBox2D& bbox, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceSegmentToTriangle(const Line2D& segment, const Triangle2D& triangle, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceSegmentToRectangle2D(const Line2D& segment, const Rectangle2D& rectangle, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceSegmentToPolygon(const Line2D& segment, const ConvexPolygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceSegmentToPolygon(const Line2D& segment, const Polygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceSegmentToCircle(const Line2D& segment, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);


real_t distanceBBoxToBBox(const BBox2D& bbox1, const BBox2D& bbox2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceBBoxToTriangle(const BBox2D& bbox1, const Triangle2D& triangle, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceBBoxToRectangle(const BBox2D& bbox1, const Rectangle2D& rectangle, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceBBoxToPolygon(const BBox2D& bbox1, const ConvexPolygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceBBoxToPolygon(const BBox2D& bbox1, const Polygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceBBoxToCircle(const BBox2D& bbox1, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);


real_t distanceTriangleToTriangle(const Triangle2D& triangle1, const Triangle2D& triangle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceTriangleToRectangle(const Triangle2D& triangle, const Rectangle2D& rectangle, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceTriangleToPolygon(const Triangle2D& triangle, const ConvexPolygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceTriangleToPolygon(const Triangle2D& triangle, const Polygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceTriangleToCircle(const Triangle2D& triangle, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);


real_t distanceRectangleToRectangle(const Rectangle2D& rectangle1, const Rectangle2D& rectangle2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceRectangleToPolygon(const Rectangle2D& rectangle, const ConvexPolygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceRectangleToPolygon(const Rectangle2D& rectangle, const Polygon2D& polygon, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distanceRectangleToCircle(const Rectangle2D& rectangle, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);


real_t distancePolygonToPolygon(const Polygon2D& polygon1, const Polygon2D& polygon2, Vector2D* closestPoint1, Vector2D* closestPoint2);

real_t distancePolygonToPolygon(const Polygon2D& polygon, const Circle2D& circle, Vector2D* closestPoint1, Vector2D* closestPoint2);


real_t distanceCircleToCircle(const Circle2D& circle1, const Circle2D& circle2, Vector2D* closestPoint1, Vector2D* closestPoint2);


// Intersection calculation algorithms

//    --- Rays ---

bool intersectRayWithBBox(const Ray2D& ray, const BBox2D& bbox, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    for (int i = 0; i < 2; i++)
    {
        if (approximatelyZero(ray.m_direction[i]))
        {
            if (ray.m_origin[i] < bbox.m_min[i] || ray.m_origin[i] > bbox.m_max[i])
                return false;
            continue;
        }

        real_t invD = 1.0f / ray.m_direction[i];
        real_t t0 = (bbox.m_min[i] - ray.m_origin[i]) * invD;
        real_t t1 = (bbox.m_max[i] - ray.m_origin[i]) * invD;

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

bool intersectRayWithCircle(const Ray2D& ray, const Circle2D& circle, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
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
        //Vector2D normal = (hitInfo->intersectionPoint - circle.m_center).normalize();
    }

    return true;
}

bool intersectRayWithTriangle(const Ray2D& ray, const Triangle2D& triangle, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    const Vector2D edge1 = triangle.m_b - triangle.m_a;
    const Vector2D edge2 = triangle.m_c - triangle.m_a;

    // Calculate determinant
    const Vector2D h = { -ray.m_direction.y, ray.m_direction.x };
    const real_t det = edge1.dot(h);

    if (approximatelyZero(det))
        return false;

    const real_t invDet = 1.0f / det;

    const Vector2D s = ray.m_origin - triangle.m_a;
    const real_t u = s.dot(h) * invDet;

    if (u < 0.0f || u > 1.0f)
        return false;

    const Vector2D q = { -s.y, s.x };
    const real_t v = ray.m_direction.dot(q) * invDet;

    if (v < 0.0f || u + v > 1.0f)
        return false;

    const real_t t = edge2.dot(q) * invDet;

    if (t < t_min || t > t_max)
        return false;

    if (hitInfo) 
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = ray.pointAt(t);
    }

    return true;
}

bool intersectRayWithSegment(const Ray2D& ray, const Vector2D& p1, const Vector2D& p2, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    Vector2D v1 = ray.m_origin - p1;
    Vector2D v2 = p2 - p1;
    Vector2D v3(-ray.m_direction.y, ray.m_direction.x);

    real_t dot = v2.dot(v3);
    
    if (approximatelyZero(dot))
        return false; // Parallel

    real_t t1 = v2.cross(v1) / dot;
    real_t t2 = v1.dot(v3) / dot;

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

bool intersectRayWithRectangleOptimized(const Ray2D& ray, const Rectangle2D& rectangle, real_t t_min, real_t t_max)
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

bool intersectRayWithRectangle(const Ray2D& ray, const Rectangle2D& rectangle, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    if (hitInfo == nullptr)
        return intersectRayWithRectangleOptimized(ray, rectangle, t_min, t_max);

    bool hit = false;
    real_t t = t_max;

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

bool intersectRayWithPolygonOptimized(const Ray2D& ray, const Polygon2D& polygon, real_t t_min, real_t t_max)
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

bool intersectRayWithPolygon(const Ray2D& ray, const Polygon2D& polygon, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    if (hitInfo == nullptr)
        return intersectRayWithPolygonOptimized(ray, polygon, t_min, t_max);

    bool hit = false;
    real_t t = t_max;

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

bool intersectRayWithShape(const Ray2D& ray, const IBaseShape2D& shape, real_t t_min, real_t t_max, HitInfo2D* hitInfo)
{
    switch(shape.type())
    {
        case ShapeType2D::SHAPE2D_CIRCLE:
            return intersectRayWithCircle(ray, dynamic_cast<const Circle2D&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_RECTANGLE:
            return intersectRayWithRectangle(ray, dynamic_cast<const Rectangle2D&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_TRIANGLE:
            return intersectRayWithTriangle(ray, dynamic_cast<const Triangle2D&>(shape), t_min, t_max, hitInfo);
        case ShapeType2D::SHAPE2D_POLYGON:
            return intersectRayWithPolygon(ray, dynamic_cast<const Polygon2D&>(shape), t_min, t_max, hitInfo);
        default:
            //Todo: we might just want to throw and error here
            return false;
    }
}

bool intersectSegmentWithSegment(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    return intersectRayWithSegment(Ray2D(p1, p2 - p1), q1, q2, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithSegment(const Line2D& line1, const Line2D& line2, HitInfo2D* hitInfo)
{
    return intersectRayWithSegment(Ray2D(line1.m_start, line1.deltaVector()), line2.m_start, line2.m_end, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithSegmentStrict(const Vector2D& p1, const Vector2D& p2, const Vector2D& q1, const Vector2D& q2, HitInfo2D* hitInfo)
{
    Vector2D r = p2 - p1;
    Vector2D s = q2 - q1;
    Vector2D qp = q1 - p1;

    real_t rxs = r.cross(s);
    real_t qpxr = qp.cross(r);

    if (approximatelyZero(rxs))
        return false; // Parallel/collinear

    real_t t = qp.cross(s) / rxs;
    real_t u = qp.cross(r) / rxs;

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

bool intersectSegmentWithBBox(const Line2D& line, const BBox2D& bbox, HitInfo2D* hitInfo)
{
    return intersectRayWithBBox(Ray2D(line.m_start, line.deltaVector()), bbox, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithCircle(const Line2D& line, const Circle2D& circle, HitInfo2D* hitInfo)
{
    return intersectRayWithCircle(Ray2D(line.m_start, line.deltaVector()), circle, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithTriangle(const Line2D& line, const Triangle2D& triangle, HitInfo2D* hitInfo)
{
    return intersectRayWithTriangle(Ray2D(line.m_start, line.deltaVector()), triangle, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithRectangle(const Line2D& line, const Rectangle2D& rectangle, HitInfo2D* hitInfo)
{
    return intersectRayWithRectangle(Ray2D(line.m_start, line.deltaVector()), rectangle, 0.f, 1.f, hitInfo);
}

bool intersectSegmentWithPolygon(const Line2D& line, const Polygon2D& polygon, HitInfo2D* hitInfo)
{
    return(intersectRayWithPolygon(Ray2D(line.m_start, line.direction()), polygon, 0.0f, line.length(), hitInfo));
}

// --- Misc Helper ---

//Todo: Specialize for direct use on BBox2D and other shapes
inline void projectPointSetOntoAxis(const Vector2D* pts, int count,
                                    const Vector2D& axis,
                                    real_t& outMin, real_t& outMax)
{
    //real_t v = dot(pts[0], axis);
    real_t v = pts[0].dot(axis);
    outMin = outMax = v;
    for (int i = 1; i < count; ++i)
    {
        //v = dot(pts[i], axis);
        v = pts[i].dot(axis);
        if (v < outMin) outMin = v;
        if (v > outMax) outMax = v;
    }
}

// --- Bounding Box ---

bool intersectBBoxWithBBox(const BBox2D& b1, const BBox2D& b2)
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

bool intersectBBoxWithTriangle(const BBox2D& bbox, const Triangle2D& triangle)
{
    Vector2D tri[3] = { triangle.m_a, triangle.m_b, triangle.m_c };

    //Todo: Projecting Axis-Aligned boxes might not be required
    Vector2D box[4] =
    {
        bbox.m_min,
        Vector2D(bbox.m_max.x, bbox.m_min.y),
        bbox.m_max,
        Vector2D(bbox.m_min.x, bbox.m_max.y)
    };

    real_t bmin, bmax, tmin, tmax;

    {
        // X-axis
        Vector2D axisX(1, 0);
        projectPointSetOntoAxis(box, 4, axisX, bmin, bmax);
        projectPointSetOntoAxis(tri, 3, axisX, tmin, tmax);
        if (!intervalsOverlap(bmin, bmax, tmin, tmax)) return false;

        // Y-axis
        Vector2D axisY(0, 1);
        projectPointSetOntoAxis(box, 4, axisY, bmin, bmax);
        projectPointSetOntoAxis(tri, 3, axisY, tmin, tmax);
        if (!intervalsOverlap(bmin, bmax, tmin, tmax)) return false;
    }

    Vector2D edges[3] =
    {
        triangle.m_b - triangle.m_a,
        triangle.m_c - triangle.m_b,
        triangle.m_a - triangle.m_c
    };

    for (int i = 0; i < 3; ++i)
    {
        // Edge normal
        Vector2D axis = edges[i].createPerpendicular();

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        projectPointSetOntoAxis(box, 4, axis, bmin, bmax);
        projectPointSetOntoAxis(tri, 3, axis, tmin, tmax);

        if (!intervalsOverlap(bmin, bmax, tmin, tmax))
            return false;
    }

    // No separating axis -> intersection
    return true;
}

bool intersectBBoxWithRectangle(const BBox2D& bbox, const Rectangle2D& rectangle)
{
    Vector2D rect[4] = { rectangle.m_a, rectangle.m_b, rectangle.m_c, rectangle.m_d };

    //Todo: Projecting Axis-Aligned boxes might not be required
    Vector2D box[4] =
    {
        bbox.m_min,
        Vector2D(bbox.m_max.x, bbox.m_min.y),
        bbox.m_max,
        Vector2D(bbox.m_min.x, bbox.m_max.y)
    };

    real_t rmin, rmax, bmin, bmax;

    {
        // X-axis
        Vector2D axisX(1, 0);
        projectPointSetOntoAxis(rect, 4, axisX, rmin, rmax);
        projectPointSetOntoAxis(box, 4, axisX, bmin, bmax);
        if (!intervalsOverlap(rmin, rmax, bmin, bmax)) return false;

        // Y-axis
        Vector2D axisY(0, 1);
        projectPointSetOntoAxis(rect, 4, axisY, rmin, rmax);
        projectPointSetOntoAxis(box, 4, axisY, bmin, bmax);
        if (!intervalsOverlap(rmin, rmax, bmin, bmax)) return false;
    }

    Vector2D e1 = rectangle.m_b - rectangle.m_a;
    Vector2D e2 = rectangle.m_c - rectangle.m_b;

    Vector2D axes[2] =
    {
        e1.createPerpendicular(), // normal to edge AB
        e2.createPerpendicular()  // normal to edge BC
    };

    for (int i = 0; i < 2; ++i)
    {
        const Vector2D& axis = axes[i];

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        projectPointSetOntoAxis(rect, 4, axis, rmin, rmax);
        projectPointSetOntoAxis(box, 4, axis, bmin, bmax);

        if (!intervalsOverlap(rmin, rmax, bmin, bmax))
            return false;
    }

    // No separating axis -> intersection
    return true;
}

bool intersectBBoxWithConvexPolygon(const BBox2D& bbox, const ConvexPolygon2D& poly)
{
    const auto& pts = poly.m_vertices;
    const size_t n = pts.size();
    if (n < 3) return false; // degenerate

    //Todo: Projecting Axis-Aligned boxes might not be required
    Vector2D box[4] =
    {
        bbox.m_min,
        Vector2D(bbox.m_max.x, bbox.m_min.y),
        bbox.m_max,
        Vector2D(bbox.m_min.x, bbox.m_max.y)
    };

    real_t pmin, pmax, bmin, bmax;

    {
        // X-axis
        Vector2D axisX(1, 0);
        projectPointSetOntoAxis(pts.data(), (int)n, axisX, pmin, pmax);
        projectPointSetOntoAxis(box, 4, axisX, bmin, bmax);
        if (!intervalsOverlap(pmin, pmax, bmin, bmax)) return false;

        // Y-axis
        Vector2D axisY(0, 1);
        projectPointSetOntoAxis(pts.data(), (int)n, axisY, pmin, pmax);
        projectPointSetOntoAxis(box, 4, axisY, bmin, bmax);
        if (!intervalsOverlap(pmin, pmax, bmin, bmax)) return false;
    }

    for (size_t i = 0; i < n; ++i)
    {
        const Vector2D& a = pts[i];
        const Vector2D& b = pts[(i + 1) % n];
        Vector2D edge = b - a;

        // Normal
        Vector2D axis = edge.createPerpendicular();

        // Skip degenerate edges
        if (axis.isZero())
            continue;

        projectPointSetOntoAxis(pts.data(), (int)n, axis, pmin, pmax);
        projectPointSetOntoAxis(box, 4, axis, bmin, bmax);

        if (!intervalsOverlap(pmin, pmax, bmin, bmax))
            return false; // separating axis found
    }

    // No separating axis -> intersection
    return true;
}

bool intersectBBoxWithPolygon(const BBox2D& bbox, const Polygon2D& polygon)
{
    if (polygon.isConvex())
    {
        const ConvexPolygon2D* convex = static_cast<const ConvexPolygon2D*>(&polygon);
        return intersectBBoxWithConvexPolygon(bbox, *convex);
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

    for (const Vector2D& p : polygon.m_vertices)
        if (bbox.contains(p))
            return true;

    Vector2D boxEdges[4][2] =
    {
        { box[0], box[1] },
        { box[1], box[2] },
        { box[2], box[3] },
        { box[3], box[0] }
    };

    const size_t n = polygon.m_vertices.size();

    for (size_t i = 0; i < n; ++i)
    {
        Vector2D a = polygon.m_vertices[i];
        Vector2D b = polygon.m_vertices[(i + 1) % n];

        for (int j = 0; j < 4; ++j)
            if (intersectSegmentWithSegment(a, b, boxEdges[j][0], boxEdges[j][1]))
                return true;
    }

    return false;
}

bool intersectBBoxWithCircle(const BBox2D& bbox, const Circle2D& circle)
{
    real_t distSq = bbox.minDistanceSquared(circle.m_center);
    real_t rSq = circle.m_radius * circle.m_radius;
    return distSq <= rSq;
}

// --- Triangles ---

bool intersectTriangleWithTriangle(const Triangle2D& triangle1, const Triangle2D& triangle2)
{
    return intersectConvexPolygonWithConvexPolygon(triangle1.getVertices(), triangle2.getVertices());
}

bool intersectTriangleWithRectangle(const Triangle2D& triangle, const Rectangle2D& rectangle)
{
    return intersectConvexPolygonWithConvexPolygon(triangle.getVertices(), rectangle.getVertices());
}

bool intersectTriangleWithPolygon(const Triangle2D& triangle, const ConvexPolygon2D& polygon)
{
    return intersectConvexPolygonWithConvexPolygon(triangle.getVertices(), polygon.getVertices());
}

bool intersectTriangleWithCircle(const Triangle2D& triangle, const Circle2D& circle)
{
    if (triangle.contains(circle.centroid()))
        return true;
    
    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    real_t r = circle.m_radius;

    if (distancePointToLine(circle.m_center, triangle.m_a, triangle.m_b) <= r)
        return true;
    if (distancePointToLine(circle.m_center, triangle.m_b, triangle.m_c) <= r)
        return true;
    if (distancePointToLine(circle.m_center, triangle.m_c, triangle.m_a) <= r)
        return true;
    return false;
}

// --- Rectangles ---

bool intersectRectangleWithRectangle(const Rectangle2D& rectangle1, const Rectangle2D& rectangle2)
{
    //Todo: Can be optimized
    return intersectConvexPolygonWithConvexPolygon(rectangle1.getVertices(), rectangle2.getVertices());
}

bool intersectRectangleWithPolygon(const Rectangle2D& rectangle, const ConvexPolygon2D& polygon)
{
    return intersectConvexPolygonWithConvexPolygon(rectangle.getVertices(), polygon.getVertices());
}

bool intersectRectangleWithCircle(const Rectangle2D& rectangle, const Circle2D& circle)
{
    if (rectangle.contains(circle.centroid()))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    real_t r = circle.m_radius;

    for (int i = 0; i < 4; ++i)
    {
        const Vector2D& a = rectangle[i];
        const Vector2D& b = rectangle[(i + 1) % 4];

        if (distancePointToLine(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

// --- Polygons --- 

bool intersectConvexPolygonWithConvexPolygon(const std::vector<Vector2D>& v1, const std::vector<Vector2D>& v2)
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

        projectPointSetOntoAxis(v1.data(), n1, axis, min1, max1);
        projectPointSetOntoAxis(v2.data(), n2, axis, min2, max2);

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

        projectPointSetOntoAxis(v1.data(), n1, axis, min1, max1);
        projectPointSetOntoAxis(v2.data(), n2, axis, min2, max2);

        if (!intervalsOverlap(min1, max1, min2, max2))
            return false;
    }

    // No separating axis -> intersection
    return true;
}

bool intersectPolygonWithPolygon(const ConvexPolygon2D& polygon1, const ConvexPolygon2D& polygon2)
{
    return intersectConvexPolygonWithConvexPolygon(polygon1.m_vertices, polygon2.m_vertices);
}

bool intersectPolygonWithCircle(const ConvexPolygon2D& polygon, const Circle2D& circle)
{
    if (polygon.contains(circle.centroid()))
        return true;

    //Todo: Optimize with r^2
    //real_t r2 = circle.m_radius * circle.m_radius;
    real_t r = circle.m_radius;

    for (int i = 0; i < polygon.vertexCount(); ++i)
    {
        const Vector2D& a = polygon.wrappedVertexAt(i);
        const Vector2D& b = polygon.wrappedVertexAt(i+1);

        if (distancePointToLine(circle.m_center, a, b) <= r)
            return true;
    }
    return false;
}

// --- Circles ---

bool intersectCircleWithCircle(const Circle2D& circle1, const Circle2D& circle2)
{
    return circle1.intersects(circle2);
}

// --- Generic ---

bool intersect(const IFiniteShape2D& shape1, const IFiniteShape2D& shape2)
{
    ShapeType2D type1 = shape1.type();
    ShapeType2D type2 = shape2.type();

    if (type1 > type2)
        return intersect(shape2, shape1);

    switch (type1)
    {
        case SHAPE2D_TRIANGLE:
            if (type2 == SHAPE2D_TRIANGLE)
            {
                const Triangle2D* s1 = shape1.shape_cast<Triangle2D>();
                const Triangle2D* s2 = shape2.shape_cast<Triangle2D>();
                return intersectTriangleWithTriangle(*s1, *s2);
            }
            else if (type2 == SHAPE2D_RECTANGLE)
            {
                const Triangle2D*  s1 = shape1.shape_cast<Triangle2D>();
                const Rectangle2D* s2 = shape2.shape_cast<Rectangle2D>();
                return intersectTriangleWithRectangle(*s1, *s2);
            }
            else if (type2 == SHAPE2D_CONVEX_POLYGON)
            {
                const Triangle2D*      s1 = shape1.shape_cast<Triangle2D>();
                const ConvexPolygon2D* s2 = shape2.shape_cast<ConvexPolygon2D>();
                return intersectTriangleWithPolygon(*s1, *s2);
            }
            else if (type2 == SHAPE2D_CIRCLE)
            {
                const Triangle2D*  s1 = shape1.shape_cast<Triangle2D>();
                const Circle2D*    s2 = shape2.shape_cast<Circle2D>();
                return intersectTriangleWithCircle(*s1, *s2);
            }
            break;
        case SHAPE2D_RECTANGLE:
            if (type2 == SHAPE2D_RECTANGLE)
            {
                const Rectangle2D* s1 = shape1.shape_cast<Rectangle2D>();
                const Rectangle2D* s2 = shape2.shape_cast<Rectangle2D>();
                return intersectRectangleWithRectangle(*s1, *s2);
            }
            else if (type2 == SHAPE2D_CONVEX_POLYGON)
            {
                const Rectangle2D*     s1 = shape1.shape_cast<Rectangle2D>();
                const ConvexPolygon2D* s2 = shape2.shape_cast<ConvexPolygon2D>();
                return intersectRectangleWithPolygon(*s1, *s2);
            }
            else if (type2 == SHAPE2D_CIRCLE)
            {
                const Rectangle2D* s1 = shape1.shape_cast<Rectangle2D>();
                const Circle2D*    s2 = shape2.shape_cast<Circle2D>();
                return intersectRectangleWithCircle(*s1, *s2);
            }
            break;
        case SHAPE2D_CONVEX_POLYGON:
            if (type2 == SHAPE2D_CONVEX_POLYGON)
            {
                const ConvexPolygon2D* s1 = shape1.shape_cast<ConvexPolygon2D>();
                const ConvexPolygon2D* s2 = shape2.shape_cast<ConvexPolygon2D>();
                return intersectPolygonWithPolygon(*s1, *s2);
            }
            else if (type2 == SHAPE2D_CIRCLE)
            {
                const ConvexPolygon2D* s1 = shape1.shape_cast<ConvexPolygon2D>();
                const Circle2D*        s2 = shape2.shape_cast<Circle2D>();
                return intersectPolygonWithCircle(*s1, *s2);
            }
            break;
        case SHAPE2D_POLYGON:
            break;
        case SHAPE2D_CIRCLE:
            if (type2 == SHAPE2D_CIRCLE)
            {
                const Circle2D* s1 = shape1.shape_cast<Circle2D>();
                const Circle2D* s2 = shape2.shape_cast<Circle2D>();
                return intersectCircleWithCircle(*s1, *s2);
            }
            break;
    }
    return false;
}

} // namespace Math

} // namespace Arns