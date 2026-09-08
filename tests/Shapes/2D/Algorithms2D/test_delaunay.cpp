#include "../../../third_party/doctest.h"

#include "Shapes/2D/Algorithms2D/Triangulation/Triangulation.hpp"

using namespace Arns::geomath;

TEST_CASE("isDelaunay accepts a valid triangle mesh")
{
    const std::vector<Vector2D> vertices{
        {0, 0},
        {4, 0},
        {0, 4}
    };
    const TriangleMesh2D mesh(vertices, {TriangleIndices(0, 1, 2)});

    CHECK(isDelaunay(mesh));
}

TEST_CASE("isDelaunay rejects a vertex inside a circumcircle")
{
    const std::vector<Vector2D> vertices{
        {0, 0},
        {4, 0},
        {0, 4},
        {1, 1}
    };
    const TriangleMesh2D mesh(vertices, {TriangleIndices(0, 1, 2)});

    CHECK(!isDelaunay(mesh));
}

TEST_CASE("isDelaunay accepts points lying EXACTLY on the circumcircle boundary")
{
    // Unit square: All 4 points lie on the circle x^2 + y^2 = 2
    const std::vector<Vector2D> vertices{
        {-1, -1}, 
        {1, -1}, 
        {1, 1}, 
        {-1, 1}
    };
    
    const TriangleMesh2D mesh(vertices, {
        TriangleIndices(0, 1, 2),
        TriangleIndices(0, 2, 3)
    });

    // Vertex 3 is on the boundary of circle(0,1,2), vertex 1 is on circle(0,2,3).
    CHECK(isDelaunay(mesh));
}

TEST_CASE("isDelaunay rejects non-Delaunay edge diagonal selection")
{
    const std::vector<Vector2D> vertices{
        {0, -2},
        {1,  0},
        {0,  0.25f},
        {-1, 0}
    };

    const TriangleMesh2D badMesh(vertices, {
        TriangleIndices(0, 1, 3),
        TriangleIndices(2, 3, 1)
    });

    // Vertex 2 is inside circumcircle of triangle(0,1,3)
    CHECK(!isDelaunay(badMesh));
}

TEST_CASE("isDelaunay rejects degenerate flat triangles (collinear vertices)")
{
    const std::vector<Vector2D> vertices{
        {0, 0},
        {1, 0},
        {2, 0}  // Collinear point, zero area
    };
    const TriangleMesh2D flatMesh(vertices, {TriangleIndices(0, 1, 2)});

    CHECK(!isDelaunay(flatMesh));
}

/*
TEST_CASE("isDelaunay rejects inverted / clockwise winding order triangles")
{
    const std::vector<Vector2D> vertices{
        {0, 0},
        {1, 0},
        {0, 1}
    };
    // Reversed winding (clockwise instead of counter-clockwise)
    const TriangleMesh2D invertedMesh(vertices, {TriangleIndices(0, 2, 1)});

    CHECK(!isDelaunay(invertedMesh));
}
*/

TEST_CASE("delaunay triangulates a square")
{
    const std::vector<Vector2D> points{
        {0, 0},
        {1, 0},
        {1, 1},
        {0, 1}
    };

    const TriangleMesh2D mesh = delaunay(points);

    CHECK(mesh.m_vertices == points);
    CHECK(mesh.m_triangles.size() == 2);
    CHECK(mesh.validate());
    CHECK(isDelaunay(mesh));
}

TEST_CASE("delaunay returns no triangles for fewer than three points")
{
    const std::vector<Vector2D> points{{0, 0}, {1, 0}};

    const TriangleMesh2D mesh = delaunay(points);

    CHECK(mesh.m_vertices == points);
    CHECK(mesh.m_triangles.empty());
    CHECK(mesh.validate());
}

TEST_CASE("delaunay handles high-density co-circular ring")
{
    const int N = 12;
    std::vector<Vector2D> points;
    for (int i = 0; i < N; ++i) 
    {
        real_t angle = real_t(2) * PI * i / N;
        points.push_back({std::cos(angle), std::sin(angle)});
    }

    const TriangleMesh2D mesh = delaunay(points);

    // Any convex polygon with N vertices decomposes into exactly N - 2 triangles
    CHECK(mesh.m_triangles.size() == static_cast<size_t>(N - 2));
    CHECK(mesh.validate());
    CHECK(isDelaunay(mesh));
}

TEST_CASE("delaunay processes structured grid")
{
    const int W = 5;
    const int H = 5;
    std::vector<Vector2D> points;
    for (int x = 0; x < W; ++x) {
        for (int y = 0; y < H; ++y) {
            points.push_back({static_cast<real_t>(x), static_cast<real_t>(y)});
        }
    }

    const TriangleMesh2D mesh = delaunay(points);

    // A rectangular grid of WxH points produces exactly 2 * (W - 1) * (H - 1) triangles
    const size_t expectedTriangles = 2 * (W - 1) * (H - 1);
    CHECK(mesh.m_triangles.size() == expectedTriangles);
    CHECK(mesh.validate());
    CHECK(isDelaunay(mesh));
}

TEST_CASE("delaunay handles near-collinear points under perturbation")
{
    const std::vector<Vector2D> points{
        {0.0, 0.0},
        {1.0, 1e-9},   // Tiny epsilon deviation from y = 0
        {2.0, -1e-9},
        {3.0, 0.0},
        {1.5, 1.0}     // One clearly off-line point to force triangulation
    };

    const TriangleMesh2D mesh = delaunay(points);

    CHECK(!mesh.m_triangles.empty());
    CHECK(mesh.validate());
    CHECK(isDelaunay(mesh));
}

TEST_CASE("delaunay ignores or safely handles duplicate points")
{
    const std::vector<Vector2D> points{
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 1.0},
        {0.0, 0.0},  // Exact duplicate of point 0
        {1.0, 0.0}   // Exact duplicate of point 1
    };

    const TriangleMesh2D mesh = delaunay(points);

    CHECK(mesh.validate());
    CHECK(isDelaunay(mesh));
}

TEST_CASE("delaunay correctly flips edges for interior point")
{
    const std::vector<Vector2D> points{
        {-2.0, -2.0},
        { 2.0, -2.0},
        { 2.0,  2.0},
        {-2.0,  2.0},
        { 0.0,  0.0}   // Center point forces 4 radial triangles
    };

    const TriangleMesh2D mesh = delaunay(points);

    CHECK(mesh.m_triangles.size() == 4);
    CHECK(mesh.validate());
    CHECK(isDelaunay(mesh));
}