#include "../../../third_party/doctest.h"
#include "Shapes/2D/Ray2D.hpp"
#include "Shapes/2D/Structures/Spatial/RTree.hpp"

using namespace Arns::geomath;

TEST_SUITE("RTree2D")
{

TEST_CASE("empty tree produces no query results")
{
	RTree2D tree;
	size_t callbackCount = 0;

	tree.query_range(BBox2D(Vector2D{0, 0}, Vector2D{10, 10}),
		[&](const ShapeID&) { ++callbackCount; return true; });
	tree.query_point(Vector2D{1, 1},
		[&](const ShapeID&) { ++callbackCount; return true;});
	tree.query_knn(Vector2D{1, 1}, 3,
		[&](const ShapeID&) { ++callbackCount; return true;});

	CHECK(callbackCount == 0);
	CHECK(tree.query_nearest(Vector2D{1, 1}) == ShapeID::invalid());
}

TEST_CASE("range and point queries use callbacks and filters")
{
	RTree2D tree(2, 1);
	const ShapeID first{1, 10};
	const ShapeID second{2, 20};
	const ShapeID third{3, 30};

	tree.build({
		{ BBox2D(Vector2D{0, 0}, Vector2D{2, 2}), first },
		{ BBox2D(Vector2D{3, 3}, Vector2D{5, 5}), second },
		{ BBox2D(Vector2D{8, 8}, Vector2D{10, 10}), third }
	});

	std::vector<ShapeID> rangeResults;
	tree.query_range(BBox2D(Vector2D{1, 1}, Vector2D{4, 4}),
		[&](const ShapeID& id) { rangeResults.push_back(id); return true; });

	REQUIRE(rangeResults.size() == 2);
	CHECK(rangeResults[0] == first);
	CHECK(rangeResults[1] == second);

	std::vector<ShapeID> filteredResults;
	tree.query_range(BBox2D(Vector2D{0, 0}, Vector2D{10, 10}),
		[&](const ShapeID& id) { filteredResults.push_back(id); return true; },
		[](const ShapeID& id) { return id.generation == 30; });

	REQUIRE(filteredResults.size() == 1);
	CHECK(filteredResults[0] == third);

	std::vector<ShapeID> pointResults;
	tree.query_point(Vector2D{3, 4},
		[&](const ShapeID& id) { pointResults.push_back(id); return true; });

	REQUIRE(pointResults.size() == 1);
	CHECK(pointResults[0] == second);
}

TEST_CASE("nearest and k-nearest queries return closest shapes")
{
	RTree2D tree(2, 1);
	const ShapeID near{1, 1};
	const ShapeID middle{2, 1};
	const ShapeID far{3, 1};

	tree.build({
		{ BBox2D(Vector2D{1, 0}, Vector2D{1, 0}), near },
		{ BBox2D(Vector2D{4, 0}, Vector2D{4, 0}), middle },
		{ BBox2D(Vector2D{9, 0}, Vector2D{9, 0}), far }
	});

	CHECK(tree.query_nearest(Vector2D{0, 0}) == near);
	CHECK(tree.query_nearest(Vector2D{0, 0},
		[](const ShapeID& id) { return id.index != 1; }) == middle);

	std::vector<ShapeID> results;
	tree.query_knn(Vector2D{0, 0}, 2,
		[&](const ShapeID& id) { results.push_back(id); return true; });

	REQUIRE(results.size() == 2);
	CHECK(results[0] == near);
	CHECK(results[1] == middle);
}

TEST_CASE("shape query callbacks support early termination")
{
	RTree2D tree(2, 1);
	tree.build({
		{ BBox2D(Vector2D{1, 1}, Vector2D{1, 1}), ShapeID{1, 1} },
		{ BBox2D(Vector2D{2, 2}, Vector2D{2, 2}), ShapeID{2, 1} },
		{ BBox2D(Vector2D{3, 3}, Vector2D{3, 3}), ShapeID{3, 1} }
	});

	size_t callbackCount = 0;
	tree.query_range(BBox2D(Vector2D{0, 0}, Vector2D{4, 4}),
		[&](const ShapeID&) { ++callbackCount; return false; });
	CHECK(callbackCount == 1);

	callbackCount = 0;
	tree.query_point(Vector2D{2, 2},
		[&](const ShapeID&) { ++callbackCount; return false; });
	CHECK(callbackCount == 1);

	callbackCount = 0;
	tree.query_knn(Vector2D{0, 0}, 3,
		[&](const ShapeID&) { ++callbackCount; return false; });
	CHECK(callbackCount == 1);
}

TEST_CASE("ray query sorts hits by t and supports early termination")
{
	RTree2D tree(2, 1);
	const ShapeID far{1, 1};
	const ShapeID near{2, 1};
	const ShapeID middle{3, 1};

	tree.build({
		{ BBox2D(Vector2D{8, -1}, Vector2D{9, 1}), far },
		{ BBox2D(Vector2D{1, -1}, Vector2D{2, 1}), near },
		{ BBox2D(Vector2D{4, -1}, Vector2D{5, 1}), middle }
	});

	Ray2D ray(Vector2D{0, 0}, Vector2D{1, 0});
	std::vector<ShapeID> results;
	std::vector<real_t> distances;
	tree.query_ray(ray, 0, 20,
		[&](ShapeID id, real_t t)
		{
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

	size_t callbackCount = 0;
	tree.query_ray(ray, 0, 20,
		[&](ShapeID, real_t)
		{
			++callbackCount;
			return false;
		});

	CHECK(callbackCount == 1);
}

}