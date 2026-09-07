#include "test_shape2D_utility.hpp"
#include "Shapes/2D/Ray2D.hpp"
#include <sstream>

TEST_CASE("Ray2D constructors and geometry")
{
	SUBCASE("Default constructor")
	{
		Ray2D ray;
		CHECK(ray.m_origin == Vector2D{0, 0});
		CHECK(ray.m_direction == Vector2D{0, 0});
	}

	SUBCASE("Origin and direction constructor")
	{
		Ray2D ray(Vector2D{1, 2}, Vector2D{3, 4});
		CHECK(ray.m_origin == Vector2D{1, 2});
		CHECK(ray.m_direction == Vector2D{3, 4});
	}

	Ray2D ray(Vector2D{1, 2}, Vector2D{3, 4});

	SUBCASE("Point and closest parameter calculations")
	{
		CHECK(ray.pointAt(0.0f) == Vector2D{1, 2});
		CHECK(ray.pointAt(0.5f) == Vector2D{2.5f, 4.0f});
		CHECK(ray.pointAt(1.0f) == Vector2D{4, 6});
		CHECK(ray.closestParameter(Vector2D{4, 6}) == doctest::Approx(1.0f));
	}
}

TEST_CASE("Ray2D transformations and comparison")
{
	Ray2D ray(Vector2D{1, 2}, Vector2D{3, 4});

	SUBCASE("Translate moves the origin")
	{
		ray.translate(Vector2D{2, -1});
		CHECK(ray.origin() == Vector2D{3, 1});
		CHECK(ray.direction() == Vector2D{3, 4});
	}

	SUBCASE("Scale around a point")
	{
		ray.scale(real_t{2}, Vector2D{0, 0});
		CHECK(ray.origin() == Vector2D{2, 4});
		CHECK(ray.direction() == Vector2D{3, 4});
	}

	SUBCASE("Transform point and direction")
	{
		ray.transform(Matrix3x3::createTranslation2D(Vector2D{2, 3}));
		CHECK(ray.origin() == Vector2D{3, 5});
		CHECK(ray.direction() == Vector2D{3, 4});
	}

	SUBCASE("Copy and comparison")
	{
		Ray2D copy = ray.copy();
		CHECK(copy == ray);
		copy.translate(Vector2D{1, 0});
		CHECK(copy != ray);
	}	
}

TEST_CASE("Ray2D derived values")
{
	Ray2D ray(Vector2D{1, 2}, Vector2D{3, 4});

	SUBCASE("Reflected direction")
	{
		Vector2D normal(0, 1);
		Ray2D reflectedRay = ray.reflected(Vector2D{1, 2}, normal);
		CHECK(reflectedRay.direction() == Vector2D{3, -4});
	}
}

TEST_CASE("Ray2D stream and format output")
{
	Ray2D ray(Vector2D{1, 2}, Vector2D{3, 4});

	std::ostringstream stream;
	stream << ray;
	CHECK(stream.str() == "Ray2D(origin: [1, 2], direction: [3, 4])");

	CHECK(std::format("{}", ray) == "Ray2D(origin: [1, 2], direction: [3, 4])");
	CHECK(std::format("{:.2f}", ray) == "Ray2D(origin: [1.00, 2.00], direction: [3.00, 4.00])");
	CHECK(to_string(ray) == "Ray2D(origin: [1, 2], direction: [3, 4])");
}