/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "../../BBox2D.hpp"
#include "../../Algorithms2D/Intersection2D.hpp"
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

struct RTreeFilter
{
    BBox2D region;
    std::vector<size_t> nodeIndices;
    std::vector<BBox2D> bounds;
};

class RTree2D : public ISpatialIndex2D
{
public:
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
            return node.count;
        }
        else
        {
            size_t subTreeSize = 0;
            for (size_t i = 0; i < node.count; ++i)
                subTreeSize += size(node.firstIndex + i);
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

    RTreeFilter createFilter(const BBox2D& region, size_t maxDepth) const;

private:

    struct Node
    {
        bool isLeaf = false;
        BBox2D bounds;
        size_t firstIndex = 0;
        size_t count = 0;
    };

    struct Candidate
	{
		real_t priority;
		bool isEntry;
		NodeIndex index;
	};

    struct CandidateCompare
	{
		bool operator()(const Candidate& left, const Candidate& right) const noexcept
		{
			if (left.priority != right.priority)
				return left.priority > right.priority;
			return left.index < right.index; // tie-breaker
		}
	};

    template <typename Queue>
	void pushRayNode(Queue& queue, const Ray2D& ray, real_t t_min, real_t t_max, NodeIndex nodeIndex) const
	{
		if (nodeIndex == InvalidNode)
			return;
		real_t enter;
		real_t exit;
		if (detail::rayBoundsInterval(ray, _nodes[nodeIndex].bounds, t_min, t_max, enter, exit))
			queue.push({ enter, false, nodeIndex });
	}

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
        return std::span<const SpatialEntry2D>(_entries).subspan(node.firstIndex, node.count);
    }

    // Recursive helper for range search
    bool rangeSearchRecursive(size_t nodeIndex, const BBox2D& query, const ShapeCallback& callback,
        const ShapeFilter& filter) const;

    bool pointSearchRecursive(size_t nodeIndex, const Vector2D& point,
        const ShapeCallback& callback, const ShapeFilter& filter) const;

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
    std::vector<Node> _nodes;
    std::vector<SpatialEntry2D> _entries;
};

} // namespace Math

} // namespace Arns