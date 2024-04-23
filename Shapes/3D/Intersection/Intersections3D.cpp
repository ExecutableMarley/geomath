/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "Intersections3D.hpp"
#include "Algorithms3D.hpp"

namespace Utility
{

namespace Math
{

bool intersects(const Line3D &line1, const Line3D &line2, Vector3D *intersection)
{
    Vector3D p1 = line1.m_start;
    Vector3D q1 = line1.m_end;
    Vector3D p2 = line2.m_start;
    Vector3D q2 = line2.m_end;

    Vector3D p1q1 = q1 - p1;
    Vector3D p2q2 = q2 - p2;
    Vector3D p1p2 = p2 - p1;

    float det = p1q1.cross(p2q2).length();

    if (det == 0)
    {
        return false;
    }

    float t = p1p2.cross(p2q2).length() / det;
    float u = p1p2.cross(p1q1).length() / det;

    if (t >= 0 && t <= 1 && u >= 0 && u <= 1)
    {
        if (intersection)
        {
            *intersection = p1 + p1q1 * t;
        }
        return true;
    }
    return false;
}


bool linePlaneIntersection(const Line3D& line, const Vector3D& planeNormal, const Vector3D& planePoint, Vector3D& intersection)
{
    Vector3D lineDirection = line.direction();
    Vector3D lineStart = line.m_start;

    float denominator = planeNormal.dot(lineDirection);

    if (denominator == 0)
    {
        return false;
    }

    float t = (planeNormal.dot(planePoint) - planeNormal.dot(lineStart)) / denominator;

    if (t > 0.f && t < 1.f)
    {
        intersection = lineStart + lineDirection * t;
        return true;
    }
    return false;
}


//Todo: There is lots of room for improvement here

bool intersects(const Line3D &line, const BBox3D &bbox, Vector3D* intersection)
{
    if (bbox.contains(line.m_start) || bbox.contains(line.m_end)) 
    {
        return true;
    }

    for (int i = 0; i < 6; i++) 
    {
        Vector3D faceNormal;
        Vector3D facePoint;

        switch (i) 
        {
            case 0: // Left face
                faceNormal = Vector3D(-1, 0, 0);
                facePoint = bbox.m_min;
                break;
            case 1: // Right face
                faceNormal = Vector3D(1, 0, 0);
                facePoint = bbox.m_max;
                break;
            case 2: // Bottom face
                faceNormal = Vector3D(0, -1, 0);
                facePoint = bbox.m_min;
                break;
            case 3: // Top face
                faceNormal = Vector3D(0, 1, 0);
                facePoint = bbox.m_max;
                break;
            case 4: // Back face
                faceNormal = Vector3D(0, 0, -1);
                facePoint = bbox.m_min;
                break;
            case 5: // Front face
                faceNormal = Vector3D(0, 0, 1);
                facePoint = bbox.m_max;
                break;
        }

        Vector3D intersectionPoint; 
        if (linePlaneIntersection(line, faceNormal, facePoint, intersectionPoint))
        {
            if (bbox.contains(intersectionPoint))
            {
                if (intersection)
                    *intersection = intersectionPoint;
                return true;
            }
        }
    }
    return false;
}

bool intersects(const Line3D &line, const Sphere &sphere, Vector3D *intersection)
{
    Vector3D cloestPoint;
    float distance = distancePointToLine(sphere.m_center, line, &cloestPoint);
    if (distance <= sphere.m_radius)
    {
        if (intersection)
        {
            *intersection = cloestPoint;
        }
        return true;
    }
    return false;
}


} // namespace Math

} // namespace Utility