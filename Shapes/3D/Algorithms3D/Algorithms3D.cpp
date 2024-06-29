/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms3D.hpp"
#include "../Ray3D.hpp"
#include "../Line3D.hpp"
#include "../BBox3D.hpp"
#include "../Plane.hpp"
#include "../Triangle3D.hpp"
#include "../Sphere.hpp"
#include "../Cylinder.hpp"
#include "../Capsule.hpp"

namespace Utility
{

namespace Math
{

float distancePointToLine(const Vector3D& point, const Vector3D& lineStart, const Vector3D& lineEnd, Vector3D* closestPoint)
{
    Vector3D v = lineEnd - lineStart;
    Vector3D w = point - lineStart;

    float c1 = w.dot(v);
    if (c1 <= 0)
    {
        if (closestPoint)
            *closestPoint = lineStart;
        return (point - lineStart).length();
    }

    float c2 = v.dot(v);
    if (c2 <= c1)
    {
        if (closestPoint)
            *closestPoint = lineEnd;
        return (point - lineEnd).length();
    }

    float b = c1 / c2;
    Vector3D tmp = lineStart + v * b;
    if (closestPoint)
        *closestPoint = tmp;
    return (point - tmp).length();   
}

float distancePointToLine(const Vector3D& point, const Line3D& line, Vector3D* closestPoint)
{
    return distancePointToLine(point, line.m_start, line.m_end, closestPoint);
}

float distanceLineToLine(const Vector3D& s1, const Vector3D& s2, const Vector3D& k1, const Vector3D& k2, Vector3D* closestPoint1, Vector3D* closestPoint2)
{
    const Vector3D   u = s2 - s1; // Line 1 Delta
	const Vector3D   v = k2 - k1; // Line 2 Delta
	const Vector3D   w = s1 - k1; // Line 1 Start - Line 2 Start

    float    a = u.dot(u); // Delta 1 squared
    float    b = u.dot(v); // Delta 1 * Delta 2
    float    c = v.dot(v); // Delta 2 squared
    float    d = u.dot(w); // Delta 1 * (Line 1 Start - Line 2 Start)
    float    e = v.dot(w); // Delta 2 * (Line 1 Start - Line 2 Start)

	float    D = a * c - b * b; // Delta 1 squared * Delta 2 squared - Delta 1 * Delta 2 squared
	float    sc, sN, sD = D;
	float    tc, tN, tD = D;

    // Parallel check
	if (approximatelyZero(D)) 
	{
		sN = 0.0f;
		sD = 1.0f;
		tN = e;
		tD = c;
	}
    //
	else 
    {
		sN = (b*e - c * d);
		tN = (a*e - b * d);
		if (sN < 0.0) 
        {
			sN = 0.0f;
			tN = e;
			tD = c;
		}
		else if (sN > sD) 
        {
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0.0f) 
    {
		tN = 0.0f;

		if (-d < 0.0f)
			sN = 0.0f;
		else if (-d > a)
			sN = sD;
		else 
        {
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) 
    {
		tN = tD;

		if ((-d + b) < 0.0f)
			sN = 0.0f;
		else if ((-d + b) > a)
			sN = sD;
		else 
        {
			sN = (-d + b);
			sD = a;
		}
	}

    sc = (approximatelyZero(sN) ? 0.0f : sN / sD);
    tc = (approximatelyZero(tN) ? 0.0f : tN / tD);

    if (closestPoint1)
        *closestPoint1 = s1 + u * sc;
    if (closestPoint2)
        *closestPoint2 = k1 + v * tc;

	Vector3D deltaP = w + (u * sc) - (v * tc);

	return deltaP.length();
}

float distanceLineToLine(const Line3D& line1, const Line3D& line2, Vector3D* closestPoint1, Vector3D* closestPoint2)
{
    return distanceLineToLine(line1.m_start, line1.m_end, line2.m_start, line2.m_end, closestPoint1, closestPoint2);
}

// Intersection calculation algorithms

struct HitInfo3D
{
    float t;
    Vector3D intersectionPoint;
    //Vector3D normal;
    //IShape3D* shape;
};

bool intersectRayWithBBox(const Ray3D& ray, const BBox3D& bbox, float t_min, float t_max, HitInfo3D* hitInfo)
{
    for (int i = 0; i < 3; i++)
    {
        float invD = 1.0f / ray.m_direction[i];
        float t0 = (bbox.m_min[i] - ray.m_origin[i]) * invD;
        float t1 = (bbox.m_max[i] - ray.m_origin[i]) * invD;

        if (invD < 0.0f)
        {
            std::swap(t0, t1);
        }

        t_min = t0 > t_min ? t0 : t_min;
        t_max = t1 < t_max ? t1 : t_max;

        if (t_max <= t_min)
        {
            return false;
        }
    }

	if (hitInfo)
	{
		hitInfo->t = t_min;
    	hitInfo->intersectionPoint = ray.pointAt(t_min);
	}

    return true;
}

bool intersectRayWithSphere(const Ray3D& ray, const Sphere& sphere, float t_min, float t_max, HitInfo3D* hitInfo)
{
    Vector3D oc = ray.m_origin - sphere.m_center;
    float a = ray.m_direction.lengthSquared();
    float b = oc.dot(ray.m_direction);
    float c = oc.lengthSquared() - sphere.m_radius * sphere.m_radius;
    float discriminant = b * b - a * c;

    if (discriminant > 0)
    {
        float temp = (-b - sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min)
        {
			if (hitInfo)
			{
				hitInfo->t = temp;
            	hitInfo->intersectionPoint = ray.pointAt(temp);
			}
            return true;
        }

        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min)
        {
			if (hitInfo)
			{
				hitInfo->t = temp;
            	hitInfo->intersectionPoint = ray.pointAt(temp);
            	//Vector3D normal = (hitInfo.intersectionPoint - sphere.m_center) / sphere.m_radius;
			}
            return true;
        }
    }
    return false;
}

bool intersectRayWithPlane(const Ray3D& ray, const Plane& plane, float t_min, float t_max, HitInfo3D* hitInfo)
{
    float denominator = plane.m_normal.dot(ray.m_direction);

    if (approximatelyZero(denominator))
    {
        return false;
    }

    float t = (plane.m_normal.dot(plane.m_planePoint) - plane.m_normal.dot(ray.m_origin)) / denominator;

    if (t > t_min && t < t_max)
    {
		if (hitInfo)
		{
			hitInfo->t = t;
        	hitInfo->intersectionPoint = ray.pointAt(t);
        	//Vector3D normal = plane.m_normal;
		}
        return true;
    }
    return false;
}

bool intersectRayWithTriangle(const Ray3D& ray, const Triangle3D& triangle, float t_min, float t_max, HitInfo3D* hitInfo)
{
    Vector3D edge1 = triangle.m_b - triangle.m_a;
    Vector3D edge2 = triangle.m_c - triangle.m_a;
    Vector3D h = ray.m_direction.cross(edge2);
    float a = edge1.dot(h);

    if (approximatelyZero(a))
    {
        return false;
    }

    float f = 1.0f / a;
    Vector3D s = ray.m_origin - triangle.m_a;
    float u = f * s.dot(h);

    if (u < 0.0f || u > 1.0f)
    {
        return false;
    }

    Vector3D q = s.cross(edge1);
    float v = f * ray.m_direction.dot(q);

    if (v < 0.0f || u + v > 1.0f)
        return false;

    float t = f * edge2.dot(q);

    if (t > t_min && t < t_max)
    {
		if (hitInfo)
		{
			hitInfo->t = t;
        	hitInfo->intersectionPoint = ray.pointAt(t);
        	//Vector3D normal = triangle.normal();
		}
        return true;
    }
    return false;
}

bool intersectRayWithCylinder(const Ray3D& ray, const Cylinder& cylinder, float t_min, float t_max, HitInfo3D* hitInfo)
{
    Vector3D d = cylinder.m_endPoint - cylinder.m_startPoint;
    Vector3D m = ray.m_origin - cylinder.m_startPoint;
    Vector3D n = ray.m_direction;
    
    float dd = d.dot(d);
    float nd = n.dot(d);
    float md = m.dot(d);
    
    float a = dd - nd * nd;
    float b = dd * m.dot(n) - md * nd;
    float c = dd * m.dot(m) - md * md - cylinder.m_radius * cylinder.m_radius * dd;
    
    float discriminant = b * b - a * c;
    
    if (discriminant < 0) return false;
    
    float sqrtDiscriminant = sqrt(discriminant);
    float t0 = (-b - sqrtDiscriminant) / a;
    float t1 = (-b + sqrtDiscriminant) / a;
    
    if (t0 > t1) std::swap(t0, t1);
    
    if (t1 < t_min || t0 > t_max) return false;
    
    float t = t0;
    if (t < t_min) t = t1;
    if (t < t_min || t > t_max) return false;
    
    Vector3D hitPoint = ray.pointAt(t);
    Vector3D projOnCylinderAxis = cylinder.m_startPoint + d * ((hitPoint - cylinder.m_startPoint).dot(d) / dd);
    
    if ((projOnCylinderAxis - cylinder.m_startPoint).dot(d) < 0 || 
        (projOnCylinderAxis - cylinder.m_endPoint).dot(d) > 0) 
    {
        return false;
    }
    
    if (hitInfo) 
    {
        hitInfo->t = t;
        hitInfo->intersectionPoint = hitPoint;
    }
    
    return true;
}

bool intersectRayWithCapsule(const Ray3D& ray, const Capsule& capsule, float t_min, float t_max)
{
    Vector3D startPoint = ray.pointAt(t_min);
    Vector3D endPoint = ray.pointAt(t_max);
    float distance = distanceLineToLine(startPoint, endPoint, capsule.m_startPoint, capsule.m_endPoint);
    if (distance <= capsule.m_radius)
    {
        return true;
    }
    return false;
}

bool intersectRayWithCapsule(const Ray3D& ray, const Capsule& capsule, float t_min, float t_max, HitInfo3D* hitInfo)
{
    Sphere sphere1(capsule.m_startPoint, capsule.m_radius);
    Sphere sphere2(capsule.m_endPoint  , capsule.m_radius);
    Cylinder cylinder(capsule.m_startPoint, capsule.m_endPoint, capsule.m_radius);

    bool hit = false;
    if (intersectRayWithSphere(ray, sphere1, t_min, t_max, hitInfo))
    {
        t_max = hitInfo->t;
        hit = true;
    }
    if (intersectRayWithSphere(ray, sphere2, t_min, t_max, hitInfo))
    {
        t_max = hitInfo->t;
        hit = true;
    }
    if (intersectRayWithCylinder(ray, cylinder, t_min, t_max, hitInfo))
    {
        hit = true;
    }
    return hit;
}

} // namespace Math

} // namespace Utility