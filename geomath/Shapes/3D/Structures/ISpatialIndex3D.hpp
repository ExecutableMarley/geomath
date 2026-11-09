/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "Containers/SlotID.hpp"
#include "../BBox3D.hpp"
#include "../Ray3D.hpp"

#include <functional>
#include <initializer_list>
#include <span>
#include <vector>

namespace Arns::geomath
{

struct SpatialEntry3D
{
    BBox3D bounds;
    ShapeID shape;
};

class ISpatialIndex3D
{
public:
    using ShapeFilter = std::function<bool(const ShapeID&)>;
    using ShapeCallback = std::function<bool(const ShapeID&)>;
    using RayHitCallback = std::function<bool(const ShapeID&, real_t)>;

    //--- [Building/Clearing] ---

    virtual void build(std::span<SpatialEntry3D> entries) = 0;

    void build(std::initializer_list<SpatialEntry3D> entries)
    {
        std::vector<SpatialEntry3D> entryCopy(entries);
        build(std::span<SpatialEntry3D>(entryCopy));
    }

    virtual void clear() = 0;

    //--- [Updating] ---

    //--- [Rebalancing/Optimizing] ---

    //--- [Query] ---

    virtual void query_range(const BBox3D& queryArea, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const = 0;

    virtual ShapeID query_nearest(const Vector3D& queryPoint,
        const ShapeFilter& filter = {}) const = 0;

    virtual void query_knn(const Vector3D& queryPoint, size_t k, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const = 0;

    virtual void query_point(const Vector3D& point, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const = 0;

    virtual void query_ray(const Ray3D& ray, real_t t_min, real_t t_max,
        const RayHitCallback& callback, const ShapeFilter& filter = {}) const = 0;

    virtual ~ISpatialIndex3D() = default;
};

} // namespace Arns::geomath