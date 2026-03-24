#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "../IShape2D.hpp"
#include "../BBox2D.hpp"

namespace Arns
{

namespace Math
{

class ISpatialIndex2D
{
public:

    //Building/Clearing

    virtual void build(const std::vector<std::pair<BBox2D, size_t>>& shapeBounds) = 0;

    virtual void insert(size_t shapeIndex, const BBox2D& bbox) = 0;

    void insert(size_t shapeIndex, const IFiniteShape2D& shape)
    {
        return insert(shapeIndex, shape.boundingBox());
    }

    virtual void clear() = 0;

    //Updating


    //Rebalancing/Optimizing


    //[Query]

    virtual void rangeQuery(const BBox2D& queryArea, std::vector<size_t>& results) const = 0;

    virtual void rangeQuery(const BBox2D& query, std::vector<size_t>& result, const std::vector<bool>& inclusionMask) const = 0;

    virtual size_t nearestNeighbour(const Vector2D& queryPoint) const = 0;

    virtual size_t nearestNeighbour(const Vector2D& queryPoint, const std::vector<bool>& inclusionMask) const = 0;

    virtual void kNearest(const Vector2D& queryPoint, size_t k, std::vector<size_t>& results) const = 0;

    virtual void kNearest(const Vector2D& queryPoint, size_t k, std::vector<size_t>& results, const std::vector<bool>& inclusionMask) const = 0;

    //virtual void bboxIntersectionQuery(const BBox2D& box, std::vector<size_t>& results) const = 0;

    //virtual void shapeIntersectionQuery(const IFiniteShape2D& shape, const std::vector<IFiniteShape2D*>& shapeStore, std::vector<size_t>& results) const = 0;

    //virtual void pointContainmentQuery(const Vector2D& point, const std::vector<IFiniteShape2D*>& shapeStore, std::vector<size_t>& results) const = 0;


    //virtual size_t size() const = 0;
    //Depth?

    virtual ~ISpatialIndex2D() = default;
};

//Todo: More filter options
//Convenient intersection checks and ray traces

//Todo: Updating/Rebalancing/Optimizing functions




} // namespace Math

} // namespace Arns