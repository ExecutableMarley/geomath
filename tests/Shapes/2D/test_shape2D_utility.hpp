#pragma once

#include "../third_party/doctest.h"

#include "../../Shapes/2D/IShape2D.hpp"

using namespace Arns::Math;

template <typename ShapeT>
void check_is_polygonal(const IFiniteShape2D& shape, size_t expectedVertexCount)
{
    SUBCASE("Detected as polygonal")
    {
        const IPolygonalShape2D* poly = shape.polygonal();
        REQUIRE(poly != nullptr);
    }

    SUBCASE("Vertex access works")
    {
        const IPolygonalShape2D* poly = shape.polygonal();
        REQUIRE(poly != nullptr);

        CHECK(poly->vertexCount() == expectedVertexCount);

        // operator[] must be usable
        CHECK_NOTHROW((*poly)[0]);
        CHECK_NOTHROW((*poly)[expectedVertexCount - 1]);
    }
}

void check_is_not_polygonal(const IFiniteShape2D& shape)
{
    SUBCASE("Not detected as polygonal")
    {
        CHECK(shape.polygonal() == nullptr);
    }
}