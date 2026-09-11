/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "../../BBox2D.hpp"
#include "../ISpatialIndex2D.hpp"

#include <cassert>
#include <queue>


namespace Arns
{

namespace geomath
{

using RTreeNodeIndex = size_t;

struct RTreeLeafEntry
{
    BBox2D bounds;
    ShapeID shape;
};

struct RTreeChildEntry
{
    BBox2D bounds;
    RTreeNodeIndex child;
};

struct RTreeDistanceEntryCompare
{
    bool operator()(const std::pair<real_t, ShapeID>& left,
        const std::pair<real_t, ShapeID>& right) const
    {
        return left.first < right.first;
    }
};

using RTreeDistanceQueue = std::priority_queue<
    std::pair<real_t, ShapeID>,
    std::vector<std::pair<real_t, ShapeID>>,
    RTreeDistanceEntryCompare>;

struct RTreeNode
{
    bool isLeaf = false;
    std::vector<RTreeLeafEntry> leafEntries;
    std::vector<RTreeChildEntry> childEntries;
};

struct RTreeFilter
{
    BBox2D region;
    std::vector<size_t> nodeIndices;
    std::vector<BBox2D> bounds;
};

class RTree : public ISpatialIndex2D
{
public:
    RTree(size_t maxEntries = 8, size_t minEntries = 3)
        : _maxEntries(maxEntries), _minEntries(minEntries), _rootIndex(std::numeric_limits<size_t>::max()) {
        assert(_minEntries <= _maxEntries && _minEntries >= 1);
    }

    // Bulk-loading method
    void build(const std::vector<std::pair<BBox2D, ShapeID>>& elementBoundsWithIndices);

    // Dynamic insertion
    void insert(ShapeID shapeIndex, const BBox2D& bounds) override {}

    // Reset method
    void clear()
    {
        this->_rootIndex = InvalidNode;
        this->_nodes.clear();
    }

    // Range search
    void query_range(const BBox2D& query, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override;

    // k-nearest Neighbours with k = 1
    ShapeID query_nearest(const Vector2D& queryPoint, const ShapeFilter& filter = {}) const override;

    // k-Nearest Neighbours
    void query_knn(const Vector2D& queryPoint, size_t k, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override;

    // Point query
    void query_point(const Vector2D& point, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override;

    // Ray query
    void query_ray(const Ray2D& ray, real_t t_min, real_t t_max, const RayHitCallback& callback,
        const ShapeFilter& filter = {}) const override;

    RTreeFilter createFilter(const BBox2D& region, size_t maxDepth) const;

    static constexpr size_t InvalidNode = std::numeric_limits<size_t>::max();
private:
    size_t _minEntries;
    size_t _maxEntries;
    size_t _rootIndex;
    BBox2D _rootBounds;
    std::vector<RTreeNode> _nodes;
    //Height

    size_t createNode(bool isLeaf)
    {
        _nodes.emplace_back();
        _nodes.back().isLeaf = isLeaf;
        return _nodes.size() - 1;
    }

    // Returns the index of subtree root
    size_t buildRecursive(std::vector<RTreeLeafEntry>& currentElements, int axis);

    template <typename Entry>
    static BBox2D calculateMBR(const std::vector<Entry>& entries)
    {
        if (entries.empty())
            return BBox2D();

        BBox2D mbr = entries[0].bounds;
        for (size_t i = 1; i < entries.size(); i++)
            mbr = mbr.encapsulate(entries[i].bounds);
        return mbr;
    }

    BBox2D calculateNodeMBR(size_t nodeIndex) const
    {
        const RTreeNode& node = _nodes[nodeIndex];
        return node.isLeaf ? calculateMBR(node.leafEntries) : calculateMBR(node.childEntries);
    }

    // Recursive helper for range search
    void rangeSearchRecursive(size_t nodeIndex, const BBox2D& query, const ShapeCallback& callback,
        const ShapeFilter& filter) const;

    void pointSearchRecursive(size_t nodeIndex, const Vector2D& point,
        const ShapeCallback& callback, const ShapeFilter& filter) const;

    void raySearchRecursive(size_t nodeIndex, const Ray2D& ray, real_t t_min, real_t t_max,
        const ShapeFilter& filter, std::vector<std::pair<real_t, ShapeID>>& hits) const;

    // retrieve k final results from priority queue
    void extractKNearestResults(size_t k, RTreeDistanceQueue& queue_to_explore,
        const ShapeCallback& callback) const;

    void collectFilterNodes(size_t nodeIndex,
        const BBox2D& region,
        size_t depth,
        size_t maxDepth,
        std::vector<size_t>& outNodes,
        std::vector<BBox2D>& outBounds) const;
};

} // namespace Math

} // namespace Arns