/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

namespace Arns::geomath
{

// Forward declarations
class Ray2D;
class Segment2D;
class BBox2D;
class IBaseShape2D;
class IFiniteShape2D;
class IPolygonalShape2D;
class Circle2D;
class Triangle2D;
class Rectangle2D;
class Polygon2D;
class ConvexPolygon2D;

namespace detail
{

bool isPointInsideConvexPolygon(const Vector2D& point, const std::span<const Vector2D> vertices);

bool isPointInsidePolygon(const Vector2D& point, const std::span<const Vector2D> vertices);

bool isPolygonInsideConvexPolygon(const std::span<const Vector2D> poly, const std::span<const Vector2D> convexPoly);

}

bool isPointOnSegment(const Vector2D& point, const Vector2D& segmentStart, const Vector2D& segmentEnd);

bool isPointOnSegment(const Vector2D& point, const Segment2D& segment);

bool isSegmentOnSegment(const Vector2D& s1, const Vector2D& s2, const Vector2D& k1, const Vector2D& k2);

bool isSegmentOnSegment(const Segment2D& segment1, const Segment2D& segment2);

// --- Point-in-Shape containment

bool isPointInsideBBox(const Vector2D& point, const BBox2D& bbox);

bool isPointInsideTriangle(const Vector2D& point, const Triangle2D& triangle);

bool isPointInsideRectangle(const Vector2D& point, const Rectangle2D& rectangle);

bool isPointInsideConvexPolygon(const Vector2D& point, const ConvexPolygon2D& polygon);

bool isPointInsidePolygon(const Vector2D& point, const Polygon2D& polygon);

//bool isPointInsidePolygon(const Vector2D& point, const IPolygonalShape2D& polygon);

bool isPointInsideCircle(const Vector2D& point, const Circle2D& circle);

//bool isPointInShape(const Vector2D& point, const IBaseShape2D& shape);

// --- Segment-in-Shape containment

bool isSegmentInsideBBox(const Segment2D &segment, const BBox2D& bbox);

bool isSegmentInsideTriangle(const Segment2D &segment, const Triangle2D& triangle);

bool isSegmentInsideRectangle(const Segment2D &segment, const Rectangle2D& rectangle);

bool isSegmentInsideConvexPolygon(const Segment2D &segment, const ConvexPolygon2D& polygon);

bool isSegmentInsidePolygon(const Segment2D &segment, const Polygon2D& polygon);

bool isSegmentInsideCircle(const Segment2D &segment, const Circle2D& circle);

//bool isSegmentInsideShape(const Segment2D& segment, const IBaseShape2D& shape);

//

bool isBBoxInsideBBox(const BBox2D& bbox1, const BBox2D& bbox2);

//bool isPolygonalInsideBBox(const IPolygonalShape2D& polygon, const BBox2D& bbox);

bool isCircleInsideBBox(const Circle2D& circle, const BBox2D& bbox);

//

bool isCircleInsideCircle(const Circle2D& circle1, const Circle2D& circle2);


} // namespace Arns::geomath