/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "../../BBox2D.hpp"
#include "../ISpatialIndex2D.hpp"

#include <cassert>
#include <limits>
#include <queue>
#include <span>
#include <vector>


namespace Arns
{

namespace geomath
{

// Todo: Bring this up to the quality of BVH2D

//using RTreeNodeIndex = size_t;

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



struct RTreeFilter
{
    BBox2D region;
    std::vector<size_t> nodeIndices;
    std::vector<BBox2D> bounds;
};

class RTree2D : public ISpatialIndex2D
{
private:
    using NodeIndex = size_t;
	static constexpr NodeIndex InvalidNode = std::numeric_limits<NodeIndex>::max();

public:
    RTree2D(size_t maxEntries = 8, size_t minEntries = 3)
        : _maxEntries(maxEntries), _minEntries(minEntries), _rootIndex(std::numeric_limits<size_t>::max()) {
        assert(_minEntries <= _maxEntries && _minEntries >= 1);
    }

    using ISpatialIndex2D::build;

    // Bulk-loading method
    void build(const std::span<SpatialEntry2D> elementBoundsWithIndices) override;

    // Dynamic insertion
    // void insert(ShapeID shapeIndex, const BBox2D& bounds) override {}

    // Reset method
    void clear() override
    {
        this->_rootIndex = InvalidNode;
        this->_nodes.clear();
        this->_entries.clear();
    }

    //

    size_t size() const
    {
        return size(_rootIndex);
    }

    size_t size(NodeIndex idx) const
    {
        const Node& node = _nodes[idx];
        if (node.isLeaf)
        {
            return node.entryCount;
        }
        else
        {
            size_t subTreeSize = 0;
            for (size_t i = 0; i < node.childCount; ++i)
                subTreeSize += size(node.firstChild + i);
            return subTreeSize;
        }
    }



    //--- [Query] ---

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

    /*
    void query_ray(const Ray2D& ray, real_t t_min, real_t t_max, const IShapeResolver2D& resolver,
        const RayHitCallback& callback, const ShapeFilter& filter = {}) const;*/

    RTreeFilter createFilter(const BBox2D& region, size_t maxDepth) const;

    //static constexpr size_t InvalidNode = std::numeric_limits<size_t>::max();
private:
    //using NodeIndex = size_t;
	//static constexpr NodeIndex InvalidNode = std::numeric_limits<NodeIndex>::max();

    struct Node
    {
        //Todo: We prob want to compress this down.
        bool isLeaf = false;
        BBox2D bounds;
        size_t firstEntry = 0;
        size_t entryCount = 0;
        NodeIndex firstChild = 0;
        size_t childCount = 0;
    };

    size_t createNode()
    {
        _nodes.emplace_back();
        return _nodes.size() - 1;
    }

    // Returns the index of subtree root
    void buildRecursive(std::span<SpatialEntry2D> currentElements, int axis, size_t nodeIndex);

    static BBox2D calculateMBR(std::span<const SpatialEntry2D> entries)
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
        return _nodes[nodeIndex].bounds;
    }

    std::span<const SpatialEntry2D> leafEntries(const Node& node) const
    {
        return std::span<const SpatialEntry2D>(_entries).subspan(node.firstEntry, node.entryCount);
    }

    // Recursive helper for range search
    bool rangeSearchRecursive(size_t nodeIndex, const BBox2D& query, const ShapeCallback& callback,
        const ShapeFilter& filter) const;

    bool pointSearchRecursive(size_t nodeIndex, const Vector2D& point,
        const ShapeCallback& callback, const ShapeFilter& filter) const;

    // retrieve k final results from priority queue
    void extractKNearestResults(size_t k, RTreeDistanceQueue& queue_to_explore,
        const ShapeCallback& callback) const;

    void collectFilterNodes(size_t nodeIndex,
        const BBox2D& region,
        size_t depth,
        size_t maxDepth,
        std::vector<size_t>& outNodes,
        std::vector<BBox2D>& outBounds) const;

private:
    size_t _minEntries;
    size_t _maxEntries;
    size_t _rootIndex;
    BBox2D _rootBounds;
    std::vector<Node> _nodes;
    std::vector<SpatialEntry2D> _entries;
};

} // namespace Math

} // namespace Arns