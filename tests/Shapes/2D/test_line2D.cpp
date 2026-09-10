#include "test_shape2D_utility.hpp"
#include "Shapes/2D/Line2D.hpp"
#include <sstream>

TEST_CASE("Segment2D constructors and geometry")
{
	SUBCASE("Default constructor")
	{
		Segment2D segment;
		CHECK(segment.m_start == Vector2D{0, 0});
		CHECK(segment.m_end == Vector2D{0, 0});
		CHECK(segment.length() == doctest::Approx(real_t{0}));
	}

	SUBCASE("Endpoint constructor")
	{
		Segment2D segment(Vector2D{1, 2}, Vector2D{4, 6});

		CHECK(segment.origin() == Vector2D{1, 2});
		CHECK(segment.deltaVector() == Vector2D{3, 4});
		CHECK(segment.length() == doctest::Approx(real_t{5}));
		CHECK(segment.direction() == Vector2D{real_t{0.6}, real_t{0.8}});
		CHECK(segment.normal() == Vector2D{real_t{0.8}, real_t{-0.6}});
	}

	SUBCASE("Direction and length constructor")
	{
		Segment2D segment(Vector2D{1, 2}, Vector2D{3, 4}, real_t{5});

		CHECK(segment.m_start == Vector2D{1, 2});
		CHECK(segment.m_end == Vector2D{4, 6});
		CHECK(segment.length() == doctest::Approx(real_t{5}));
	}

	SUBCASE("Point and closest parameter")
	{
		Segment2D segment(Vector2D{0, 0}, Vector2D{4, 0});

		CHECK(segment.pointAt(real_t{0}) == Vector2D{0, 0});
		CHECK(segment.pointAt(real_t{0.5}) == Vector2D{2, 0});
		CHECK(segment.pointAt(real_t{1}) == Vector2D{4, 0});
		CHECK(segment.closestParameter(Vector2D{2, 3}) == doctest::Approx(real_t{0.5}));
	}
}

TEST_CASE("Segment2D transformations and comparison")
{
	Segment2D segment(Vector2D{0, 0}, Vector2D{2, 2});

	SUBCASE("Translate")
	{
		segment.translate(Vector2D{1, -1});
		CHECK(segment.m_start == Vector2D{1, -1});
		CHECK(segment.m_end == Vector2D{3, 1});
	}

	SUBCASE("Scale around a point")
	{
		segment.scale(real_t{2}, Vector2D{0, 0});
		CHECK(segment.m_start == Vector2D{0, 0});
		CHECK(segment.m_end == Vector2D{4, 4});
	}

	SUBCASE("Rotate around a point")
	{
		segment.rotate(real_t{90}, Vector2D{0, 0});
		CHECK(segment.m_start == Vector2D{0, 0});
		CHECK(segment.m_end.x == doctest::Approx(real_t{-2}));
		CHECK(segment.m_end.y == doctest::Approx(real_t{2}));
	}

	SUBCASE("Transform")
	{
		segment.transform(Matrix3x3::createTranslation2D(Vector2D{3, 4}));
		CHECK(segment.m_start == Vector2D{3, 4});
		CHECK(segment.m_end == Vector2D{5, 6});
	}

	SUBCASE("Copy and comparison")
	{
		Segment2D copy = segment.copy();
		CHECK(copy == segment);
		copy.translate(Vector2D{1, 0});
		CHECK(copy != segment);
	}
}

TEST_CASE("Segment2D stream and format output")
{
	Segment2D segment(Vector2D{1, 2}, Vector2D{3, 4});

	std::ostringstream stream;
	stream << segment;
	CHECK(stream.str() == "Segment2D(start: [1, 2], end: [3, 4])");

	CHECK(std::format("{}", segment) == "Segment2D(start: [1, 2], end: [3, 4])");
	CHECK(std::format("{:.2f}", segment) == "Segment2D(start: [1.00, 2.00], end: [3.00, 4.00])");
	CHECK(to_string(segment) == "Segment2D(start: [1, 2], end: [3, 4])");
}