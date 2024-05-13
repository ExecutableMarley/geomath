/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Algorithms3D.hpp"
#include "../Line3D.hpp"

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



} // namespace Math

} // namespace Utility