#include "../../../third_party/doctest.h"
#include "Shapes/2D/Ray2D.hpp"
#include "Shapes/2D/Structures/Spatial/BVH2D.hpp"
#include "Shapes/2D/Structures/Spatial/RTree.hpp"
#include "Shapes/2D/Structures/Spatial/TransientQuadTree2D.hpp"

using namespace Arns::geomath;

TYPE_TO_STRING(BVH2D);
TYPE_TO_STRING(RTree2D);
TYPE_TO_STRING(TransientQuadTree2D);

template <typename SpatialIndex>
struct SpatialIndexFactory
{
    static SpatialIndex create()
    {
        return SpatialIndex{};
    }
};


template <>
struct SpatialIndexFactory<TransientQuadTree2D>
{
    static TransientQuadTree2D create()
    {
        return TransientQuadTree2D(BBox2D({0, 0}, {10, 10}), 8, 8);
    }
};

TEST_CASE_TEMPLATE_DEFINE("ISpatialIndex2D Contract Tests", SpatialIndex, SpatialIndexContract)
{
    SUBCASE("empty index produces no query results")
    {
        SpatialIndex index = SpatialIndexFactory<SpatialIndex>::create();
        size_t callbackCount = 0;

        index.query_range(BBox2D({0, 0}, {10, 10}), [&](const ShapeID&) { ++callbackCount; return true; });
        index.query_point(Vector2D{1, 1}, [&](const ShapeID&) { ++callbackCount; return true; });
        index.query_knn(Vector2D{1, 1}, 3, [&](const ShapeID&) { ++callbackCount; return true; });
        index.query_ray(Ray2D({0, 0}, {1, 0}), 0, 100, [&](const ShapeID&, real_t) { ++callbackCount; return true; });

        CHECK(callbackCount == 0);
        CHECK(index.query_nearest(Vector2D{1, 1}) == ShapeID::invalid());
    }

    SUBCASE("clear() resets index state completely")
    {
        SpatialIndex index = SpatialIndexFactory<SpatialIndex>::create();
        index.build({
            { BBox2D({0, 0}, {2, 2}), ShapeID{1, 1} },
            { BBox2D({3, 3}, {5, 5}), ShapeID{2, 1} }
        });

        index.clear();

        size_t hits = 0;
        index.query_range(BBox2D({0, 0}, {10, 10}), [&](const ShapeID&) { ++hits; return true; });
        
        CHECK(hits == 0);
        CHECK(index.query_nearest(Vector2D{1, 1}) == ShapeID::invalid());
    }

    SUBCASE("range and point queries respect spatial boundaries and filters")
    {
        SpatialIndex index = SpatialIndexFactory<SpatialIndex>::create();
        const ShapeID first{1, 10};
        const ShapeID second{2, 20};
        const ShapeID third{3, 30};

        index.build({
            { BBox2D({0, 0}, {2, 2}), first },
            { BBox2D({3, 3}, {5, 5}), second },
            { BBox2D({8, 8}, {10, 10}), third }
        });

        // Range Query
        std::vector<ShapeID> rangeResults;
        index.query_range(BBox2D({1, 1}, {4, 4}), [&](const ShapeID& id) {
            rangeResults.push_back(id);
            return true;
        });

        REQUIRE(rangeResults.size() == 2);
        CHECK_MESSAGE(std::find(rangeResults.begin(), rangeResults.end(), first) != rangeResults.end(), "Should contain first entry");
        CHECK_MESSAGE(std::find(rangeResults.begin(), rangeResults.end(), second) != rangeResults.end(), "Should contain second entry");

        // Filtered Query
        std::vector<ShapeID> filteredResults;
        index.query_range(BBox2D({0, 0}, {10, 10}),
            [&](const ShapeID& id) { filteredResults.push_back(id); return true; },
            [](const ShapeID& id) { return id.generation == 30; }
        );

        REQUIRE(filteredResults.size() == 1);
        CHECK(filteredResults[0] == third);

        // Point Query
        std::vector<ShapeID> pointResults;
        index.query_point(Vector2D{3.5, 3.5}, [&](const ShapeID& id) {
            pointResults.push_back(id);
            return true;
        });

        REQUIRE(pointResults.size() == 1);
        CHECK(pointResults[0] == second);
    }

    SUBCASE("nearest and k-nearest return results ordered by distance")
    {
        SpatialIndex index = SpatialIndexFactory<SpatialIndex>::create();
        const ShapeID near{1, 1};
        const ShapeID middle{2, 1};
        const ShapeID far{3, 1};

        index.build({
            { BBox2D({1, 0}, {1, 0}), near },
            { BBox2D({4, 0}, {4, 0}), middle },
            { BBox2D({9, 0}, {9, 0}), far }
        });

        CHECK(index.query_nearest(Vector2D{0, 0}) == near);

        // Filtered Nearest
        CHECK(index.query_nearest(Vector2D{0, 0}, [](const ShapeID& id) { return id.index != 1; }) == middle);

        // KNN Order check
        std::vector<ShapeID> results;
        index.query_knn(Vector2D{0, 0}, 2, [&](const ShapeID& id) {
            results.push_back(id);
            return true;
        });

        REQUIRE(results.size() == 2);
        CHECK(results[0] == near);
        CHECK(results[1] == middle);
    }

    SUBCASE("queries support early termination via callback return value")
    {
        SpatialIndex index = SpatialIndexFactory<SpatialIndex>::create();
        index.build({
            { BBox2D({1, 1}, {1, 1}), ShapeID{1, 1} },
            { BBox2D({2, 2}, {2, 2}), ShapeID{2, 1} },
            { BBox2D({3, 3}, {3, 3}), ShapeID{3, 1} }
        });

        auto test_early_exit = [&](auto query_fn) {
            size_t count = 0;
            query_fn([&](const ShapeID&) {
                ++count;
                return false; // Stop after first match
            });
            CHECK(count == 1);
        };

        test_early_exit([&](auto cb) { index.query_range(BBox2D({0, 0}, {4, 4}), cb); });
        test_early_exit([&](auto cb) { index.query_point(Vector2D{2, 2}, cb); });
        test_early_exit([&](auto cb) { index.query_knn(Vector2D{0, 0}, 3, cb); });
    }

    SUBCASE("ray query traverses along ray direction, returns sorted t values, and respects t limits")
    {
        SpatialIndex index = SpatialIndexFactory<SpatialIndex>::create();
        const ShapeID near{1, 1};
        const ShapeID middle{2, 1};
        const ShapeID far{3, 1};

        index.build({
            { BBox2D({8, -1}, {9, 1}), far },
            { BBox2D({1, -1}, {2, 1}), near },
            { BBox2D({4, -1}, {5, 1}), middle }
        });

        Ray2D ray(Vector2D{0, 0}, Vector2D{1, 0});
        
        // 1. Check ordering by t distance
        std::vector<ShapeID> results;
        std::vector<real_t> distances;
        index.query_ray(ray, 0, 20, [&](ShapeID id, real_t t) {
            results.push_back(id);
            distances.push_back(t);
            return true;
        });

        REQUIRE(results.size() == 3);
        CHECK(results[0] == near);
        CHECK(results[1] == middle);
        CHECK(results[2] == far);
        CHECK(distances[0] < distances[1]);
        CHECK(distances[1] < distances[2]);

        // 2. Check t_min / t_max bounds filtering
        results.clear();
        index.query_ray(ray, 3.0, 6.0, [&](ShapeID id, real_t) {
            results.push_back(id);
            return true;
        });
        REQUIRE(results.size() == 1);
        CHECK(results[0] == middle);
    }
}

// Register all spatial structures against the test template suite
TEST_CASE_TEMPLATE_INVOKE(SpatialIndexContract, BVH2D, RTree2D, TransientQuadTree2D);