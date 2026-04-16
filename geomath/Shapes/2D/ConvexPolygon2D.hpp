/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "CommonMath.hpp"
#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"
#include "Polygon2D.hpp"

namespace Arns
{

namespace Math
{

class ConvexPolygon2D : public Polygon2D
{
public:
    ConvexPolygon2D() = default;

    //Todo: This is bad
    ConvexPolygon2D(const std::vector<Vector2D>& points)
    {
        this->m_vertices = convex_hull(points);
    }

    ShapeType2D type() const { return SHAPE2D_CONVEX_POLYGON; }
    static constexpr ShapeType2D shapeType = SHAPE2D_CONVEX_POLYGON;

private:

    static std::vector<Vector2D> convex_hull(const std::vector<Vector2D>& v)
    {
        auto points = v;

        sort(points.begin(), points.end(), [](const Vector2D& a, const Vector2D& b)
            {
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            });

        points.erase(std::unique(points.begin(), points.end(), [](Vector2D a, Vector2D b)
            { return a == b; }), points.end());

        size_t hullSize = 0;
        size_t n = points.size();

        if (n <= 1)
            return points;
        
	    std::vector<Vector2D> hull(n*2);

	    for (size_t i = 0; i < n; ++i)
        {
	    	while (hullSize >= 2 && orient2D(hull[hullSize-2], hull[hullSize-1], points[i]) <= 0)
                hullSize--;
	    	hull[hullSize++] = points[i];
	    }

	    for (size_t i = n-1, t = hullSize+1; i > 0; --i)
        {
	    	while (hullSize >= t && orient2D(hull[hullSize-2], hull[hullSize-1], points[i-1]) <= 0)
                hullSize--;
	    	hull[hullSize++] = points[i-1];
	    }

	    hull.resize(hullSize-1);
	    return hull;
    }

public:
    static ConvexPolygon2D fromPoints(const std::vector<Vector2D>& points)
    {
        ConvexPolygon2D convexPolygon;
        convexPolygon.m_vertices = convex_hull(points);
        return convexPolygon;
    }
};

} // namespace Math

} // namespace Arns