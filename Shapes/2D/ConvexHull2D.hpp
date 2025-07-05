/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>

#include "Geometry/Vector2D.hpp"
#include "BBox2D.hpp"
#include "IShape2D.hpp"
#include "Polygon.hpp"

namespace Utility
{

namespace Math
{

//Todo: Consider renaming to convexPolygon
class ConvexPolygon2D : public Polygon2D
{
public:
    ConvexPolygon2D() = default;

    ConvexPolygon2D(const std::vector<Vector2D> points)
    {
        this->m_vertices = convex_hull(points);
    }

private:
    //Todo: Consider moving this to Vector2D file
    float cross(const Vector2D& a, const Vector2D& b, const Vector2D& c)
    {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    std::vector<Vector2D> convex_hull(std::vector<Vector2D> P)
    {
        size_t n = P.size(), k = 0;
	    
        //Check for duplicates
        //Handle n <= 3 better

        if (n <= 3) 
            return P;
        
	    std::vector<Vector2D> H(2*n);

	    // Sort points lexicographically
        sort(P.begin(), P.end(), [](Vector2D a, Vector2D b)
            {
                return a.x < b.x || (a.x == b.x && a.y < b.y);
            });

	    // Build lower hull
	    for (size_t i = 0; i < n; ++i)
        {
	    	while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
	    	H[k++] = P[i];
	    }

	    // Build upper hull
	    for (size_t i = n-1, t = k+1; i > 0; --i)
        {
	    	while (k >= t && cross(H[k-2], H[k-1], P[i-1]) <= 0) k--;
	    	H[k++] = P[i-1];
	    }

	    H.resize(k-1);
	    return H;
    }
};

}

}

/*
https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
https://github.com/JernejPuc/convex-hull
*/