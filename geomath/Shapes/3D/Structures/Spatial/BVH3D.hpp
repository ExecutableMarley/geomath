/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "../../BBox3D.hpp"
#include "../../Algorithms3D/Intersection3D.hpp"
#include "../ISpatialIndex3D.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <vector>

namespace Arns::geomath
{


class BVH3D : public ISpatialIndex3D
{
public:
    explicit BVH3D(size_t maxLeafEntries = 4)
        : _maxLeafEntries(std::max<size_t>(1, maxLeafEntries)) {}

    using ISpatialIndex3D::build;

    //--- [Building/Clearing] ---

    void build(std::span<SpatialEntry3D> entries) override
    {
        clear();
        _entries.assign(entries.begin(), entries.end());
        if (!_entries.empty())
            _rootIndex = buildNode(0, _entries.size());
    }

    void clear() override
    {
        _entries.clear();
        _nodes.clear();
        _rootIndex = InvalidNode;
    }

    //

    size_t size() const noexcept { return _entries.size(); }

    //--- [Query] ---

    void query_range(const BBox3D& query, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override
    {
        if (callback && _rootIndex != InvalidNode)
            rangeSearchRecursive(_rootIndex, query, callback, filter);
    }

    ShapeID query_nearest(const Vector3D& point, const ShapeFilter& filter = {}) const override
    {
        if (_rootIndex == InvalidNode)
            return ShapeID::invalid();

        std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare> queue;
        queue.push({_nodes[_rootIndex].bounds.minDistanceSquared(point), false, _rootIndex});
        while (!queue.empty())
        {
            const Candidate candidate = queue.top();
            queue.pop();
            if (candidate.isEntry)
                return _entries[candidate.index].shape;

            const Node& node = _nodes[candidate.index];
            for (size_t i = 0; i < node.count; ++i)
                pushEntry(queue, node.firstEntry + i, point, filter);
            pushNode(queue, node.left, point);
            pushNode(queue, node.right, point);
        }
        return ShapeID::invalid();
    }

    void query_knn(const Vector3D& point, size_t k, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override
    {
        if (!callback || k == 0 || _rootIndex == InvalidNode)
            return;

        std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare> queue;
        queue.push({_nodes[_rootIndex].bounds.minDistanceSquared(point), false, _rootIndex});
        size_t emitted = 0;
        while (!queue.empty() && emitted < k)
        {
            const Candidate candidate = queue.top();
            queue.pop();
            if (candidate.isEntry)
            {
                if (!callback(_entries[candidate.index].shape))
                    return;
                ++emitted;
                continue;
            }

            const Node& node = _nodes[candidate.index];
            for (size_t i = 0; i < node.count; ++i)
                pushEntry(queue, node.firstEntry + i, point, filter);
            pushNode(queue, node.left, point);
            pushNode(queue, node.right, point);
        }
    }

    void query_point(const Vector3D& point, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override
    {
        if (callback && _rootIndex != InvalidNode)
            pointSearchRecursive(_rootIndex, point, callback, filter);
    }

    void query_ray(const Ray3D& ray, real_t t_min, real_t t_max,
        const RayHitCallback& callback, const ShapeFilter& filter = {}) const override
    {
        if (!callback || _rootIndex == InvalidNode || t_min > t_max)
            return;

        std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare> queue;
        pushRayNode(queue, ray, t_min, t_max, _rootIndex);
        while (!queue.empty())
        {
            const Candidate candidate = queue.top();
            queue.pop();
            if (candidate.isEntry)
            {
                if (!callback(_entries[candidate.index].shape, candidate.priority))
                    return;
                continue;
            }

            const Node& node = _nodes[candidate.index];
            for (uint32_t i = 0; i < node.count; ++i)
            {
                const uint32_t entryIndex = node.firstEntry + i;
                const SpatialEntry3D& entry = _entries[entryIndex];
                if (filter && !filter(entry.shape))
                    continue;

                real_t enter;
                real_t exit;
                if (detail::rayBoundsInterval(ray, entry.bounds, t_min, t_max, enter, exit))
                    queue.push({enter, true, entryIndex});
            }
            pushRayNode(queue, ray, t_min, t_max, node.left);
            pushRayNode(queue, ray, t_min, t_max, node.right);
        }
    }

private:
    using NodeIndex = uint32_t;
    static constexpr NodeIndex InvalidNode = std::numeric_limits<NodeIndex>::max();

    struct Node
    {
        BBox3D bounds;
        NodeIndex left = InvalidNode;
        NodeIndex right = InvalidNode;
        uint32_t firstEntry = 0;
        uint32_t count = 0;
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
            return left.priority > right.priority;
        }
    };

    //Todo: Figure out SAH
    NodeIndex buildNode(size_t first, size_t last)
    {
        Node node;
        node.bounds = _entries[first].bounds;
        node.firstEntry = first;
        node.count = last - first;
        BBox3D centroidBounds(_entries[first].bounds.centroid(), _entries[first].bounds.centroid());
        for (size_t i = first + 1; i < last; ++i)
        {
            node.bounds.encapsulate(_entries[i].bounds);
            centroidBounds.encapsulate(_entries[i].bounds.centroid());
        }

        const NodeIndex nodeIndex = _nodes.size();
        _nodes.push_back(node);
        if (node.count <= _maxLeafEntries)
            return nodeIndex;

        const Vector3D centroidSize = centroidBounds.size();
        size_t splitAxis = 0;
        if (centroidSize.y > centroidSize.x)
            splitAxis = 1;
        if (centroidSize.z > centroidSize[splitAxis])
            splitAxis = 2;

        const size_t middle = first + node.count / 2;
        std::nth_element(_entries.begin() + first, _entries.begin() + middle, _entries.begin() + last,
            [splitAxis](const SpatialEntry3D& left, const SpatialEntry3D& right) {
                return left.bounds.centroid()[splitAxis] < right.bounds.centroid()[splitAxis];
            });

        _nodes[nodeIndex].count = 0;
        _nodes[nodeIndex].left = buildNode(first, middle);
        _nodes[nodeIndex].right = buildNode(middle, last);
        return nodeIndex;
    }

    template <typename Queue>
    void pushEntry(Queue& queue, size_t entryIndex, const Vector3D& point,
        const ShapeFilter& filter) const
    {
        const SpatialEntry3D& entry = _entries[entryIndex];
        if (!filter || filter(entry.shape))
            queue.push({entry.bounds.minDistanceSquared(point), true, entryIndex});
    }

    template <typename Queue>
    void pushNode(Queue& queue, NodeIndex nodeIndex, const Vector3D& point) const
    {
        if (nodeIndex != InvalidNode)
            queue.push({_nodes[nodeIndex].bounds.minDistanceSquared(point), false, nodeIndex});
    }

    template <typename Queue>
    void pushRayNode(Queue& queue, const Ray3D& ray, real_t t_min, real_t t_max,
        NodeIndex nodeIndex) const
    {
        if (nodeIndex == InvalidNode)
            return;
        real_t enter;
        real_t exit;
        if (detail::rayBoundsInterval(ray, _nodes[nodeIndex].bounds, t_min, t_max, enter, exit))
            queue.push({enter, false, nodeIndex});
    }

    bool rangeSearchRecursive(NodeIndex nodeIndex, const BBox3D& query,
        const ShapeCallback& callback, const ShapeFilter& filter) const
    {
        const Node& node = _nodes[nodeIndex];
        if (!query.intersects(node.bounds))
            return true;
        for (size_t i = 0; i < node.count; ++i)
        {
            const SpatialEntry3D& entry = _entries[node.firstEntry + i];
            if ((!filter || filter(entry.shape)) && query.intersects(entry.bounds) && !callback(entry.shape))
                return false;
        }
        return (node.left == InvalidNode || rangeSearchRecursive(node.left, query, callback, filter)) &&
            (node.right == InvalidNode || rangeSearchRecursive(node.right, query, callback, filter));
    }

    bool pointSearchRecursive(NodeIndex nodeIndex, const Vector3D& point,
        const ShapeCallback& callback, const ShapeFilter& filter) const
    {
        const Node& node = _nodes[nodeIndex];
        if (!node.bounds.contains(point))
            return true;
        for (size_t i = 0; i < node.count; ++i)
        {
            const SpatialEntry3D& entry = _entries[node.firstEntry + i];
            if ((!filter || filter(entry.shape)) && entry.bounds.contains(point) && !callback(entry.shape))
                return false;
        }
        return (node.left == InvalidNode || pointSearchRecursive(node.left, point, callback, filter)) &&
            (node.right == InvalidNode || pointSearchRecursive(node.right, point, callback, filter));
    }

    size_t _maxLeafEntries;
    NodeIndex _rootIndex = InvalidNode;
    std::vector<Node> _nodes;
    std::vector<SpatialEntry3D> _entries;
};

} // namespace Arns::geomath