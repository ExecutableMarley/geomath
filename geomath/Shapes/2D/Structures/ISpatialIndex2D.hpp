#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "../IShape2D.hpp"
#include "../BBox2D.hpp"
#include "ShapeStore2D.hpp"

namespace Arns
{

namespace Math
{

class ISpatialIndex2D
{
public:

    //Building/Clearing

    virtual void build(const std::vector<std::pair<BBox2D, ShapeID>>& shapeBounds) = 0;

    virtual void insert(ShapeID shapeIndex, const BBox2D& bbox) = 0;

    void insert(ShapeID shapeIndex, const IFiniteShape2D& shape)
    {
        return insert(shapeIndex, shape.boundingBox());
    }

    virtual void clear() = 0;

    //Updating


    //Rebalancing/Optimizing


    //[Query]

    virtual void rangeQuery(const BBox2D& queryArea, std::vector<ShapeID>& results) const = 0;

    virtual void rangeQuery(const BBox2D& query, std::vector<ShapeID>& result, const std::vector<bool>& inclusionMask) const = 0;

    virtual ShapeID nearestNeighbour(const Vector2D& queryPoint) const = 0;

    virtual ShapeID nearestNeighbour(const Vector2D& queryPoint, const std::vector<bool>& inclusionMask) const = 0;

    virtual void kNearest(const Vector2D& queryPoint, size_t k, std::vector<ShapeID>& results) const = 0;

    virtual void kNearest(const Vector2D& queryPoint, size_t k, std::vector<ShapeID>& results, const std::vector<bool>& inclusionMask) const = 0;

    using CandidateCallback = std::function<bool(ShapeID)>;

    virtual void queryRange(const BBox2D& queryArea, const CandidateCallback& callback) const = 0;
    virtual void queryNearest(const Vector2D& queryPoint, const CandidateCallback& callback) const = 0;

    //virtual void bboxIntersectionQuery(const BBox2D& box, std::vector<size_t>& results) const = 0;

    //virtual void shapeIntersectionQuery(const IFiniteShape2D& shape, const ShapeStore& shapeStore, std::vector<size_t>& results) const = 0;

    //virtual void pointContainmentQuery(const Vector2D& point, const ShapeStore& shapeStore, std::vector<size_t>& results) const = 0;

    //Ray trace query


    //virtual size_t size() const = 0;
    //Depth?

    virtual ~ISpatialIndex2D() = default;
};

//Todo: More filter options
//Convenient intersection checks and ray traces

//Todo: Updating/Rebalancing/Optimizing functions




} // namespace Math

} // namespace Arns