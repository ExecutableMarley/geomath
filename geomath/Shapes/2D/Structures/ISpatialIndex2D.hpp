#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "../Interfaces/IFiniteShape2D.hpp"
#include "../BBox2D.hpp"
#include "../Ray2D.hpp"
#include "ShapeStore2D.hpp"

namespace Arns
{

namespace geomath
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

    using ShapeFilter = std::function<bool(const ShapeID&)>;
    using ShapeCallback = std::function<void(const ShapeID&)>;
    using RayHitCallback = std::function<bool(ShapeID id, real_t t)>;

    virtual void query_range(const BBox2D& queryArea, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const = 0;

    virtual ShapeID query_nearest(const Vector2D& queryPoint, const ShapeFilter& filter = {}) const = 0;

    virtual void query_knn(const Vector2D& queryPoint, size_t k, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const = 0;

    virtual void query_point(const Vector2D& point, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const = 0;

    virtual void query_ray(const Ray2D& ray, real_t t_min, real_t t_max, const RayHitCallback& callback,
        const ShapeFilter& filter = {}) const = 0;


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