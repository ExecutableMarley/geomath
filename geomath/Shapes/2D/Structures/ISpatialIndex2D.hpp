#pragma once

#include "../Interfaces/IFiniteShape2D.hpp"
#include "../BBox2D.hpp"
#include "../Ray2D.hpp"
#include "ShapeStore2D.hpp"

#include <initializer_list>
#include <memory>
#include <vector>
#include <span>
#include <functional>

namespace Arns
{

namespace geomath
{


struct SpatialEntry
{
    BBox2D bounds;
    ShapeID shape;
};

class ISpatialIndex2D
{
public:

    //--- [Building/Clearing] ---

    //virtual void build(const std::vector<SpatialEntry>& shapeBounds) = 0;

    virtual void build(const std::span<SpatialEntry> entries) = 0;

    void build(std::initializer_list<SpatialEntry> entries)
    {
        std::vector<SpatialEntry> _entries(entries);
        build(std::span<SpatialEntry>(_entries));
    }

    /*
    virtual void insert(ShapeID shapeIndex, const BBox2D& bbox) = 0;

    void insert(ShapeID shapeIndex, const IFiniteShape2D& shape)
    {
        return insert(shapeIndex, shape.boundingBox());
    }
    */

    virtual void clear() = 0;

    //--- [Updating] ---

    //--- [Rebalancing/Optimizing] ---

    //--- [Query] ---

    using ShapeFilter    = std::function<bool(const ShapeID& id)>;
    using ShapeCallback  = std::function<bool(const ShapeID& id)>;
    using RayHitCallback = std::function<bool(const ShapeID& id, real_t t)>;

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