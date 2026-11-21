/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "RTree.hpp"

#include "../../Algorithms2D/Intersection2D.hpp"
#include <queue>


namespace Arns
{

namespace geomath
{

namespace
{

/*
struct RTreeRayQueueEntry
{
    real_t t;
    size_t nodeIndex;
    ShapeID shape;
    bool isNode;
    size_t sequence;
};

struct RTreeRayQueueCompare
{
    bool operator()(const RTreeRayQueueEntry& left, const RTreeRayQueueEntry& right) const
    {
        if (left.t != right.t)
            return left.t > right.t;
        return left.sequence > right.sequence;
    }
};

using RTreeRayQueue = std::priority_queue<
    RTreeRayQueueEntry,
    std::vector<RTreeRayQueueEntry>,
    RTreeRayQueueCompare>;
*/
}


void RTree2D::build(const std::span<SpatialEntry2D> elementBoundsWithIndices)
{
    this->clear();

    //Consider reserving space in _nodes
    if (elementBoundsWithIndices.empty())
    {
        _rootIndex = InvalidNode;
        return;
    }

    // mutable workspace
    std::vector<SpatialEntry2D> buildWorkspace(elementBoundsWithIndices.begin(), elementBoundsWithIndices.end());
    _rootIndex = createNode();
    buildRecursive(std::span<SpatialEntry2D>(buildWorkspace), 0, _rootIndex);
}

void RTree2D::buildRecursive(const std::span<SpatialEntry2D> currentElements, int axis, size_t nodeIndex)
{
    // Base Case
    if (currentElements.size() <= _maxEntries)
    {
        Node& node = _nodes[nodeIndex];
        node.isLeaf = true;
        node.bounds = calculateMBR(currentElements);
        node.firstIndex = _entries.size();
        node.count = currentElements.size();
        _entries.insert(_entries.end(), currentElements.begin(), currentElements.end());
        return;
    }

    size_t numElements = currentElements.size();
    size_t numSlices = static_cast<size_t>(std::ceil(std::sqrt(static_cast<real_t>(numElements) / _maxEntries)));
    if (numSlices == 0) numSlices = 1;

    size_t elementsPerSlice = static_cast<size_t>(std::ceil(static_cast<real_t>(numElements) / numSlices));
    if (elementsPerSlice == 0) elementsPerSlice = 1;

    if (axis == 0)
    {
        // Sort by X-coordinate
        std::sort(currentElements.begin(), currentElements.end(),
            [](const SpatialEntry2D& a, const SpatialEntry2D& b)
            {
                return a.bounds.m_min.x < b.bounds.m_min.x;
            });
    }
    else
    {
        // Sort by Y-coordinate
        std::sort(currentElements.begin(), currentElements.end(),
            [](const SpatialEntry2D& a, const SpatialEntry2D& b)
            {
                return a.bounds.m_min.y < b.bounds.m_min.y;
            });
    }

    const size_t firstChild = _nodes.size();
    for (size_t i = 0; i < numSlices; ++i)
        createNode();

    Node& node = _nodes[nodeIndex];
    node.isLeaf = false;
    node.firstIndex = firstChild;
    node.count = numSlices;

    for (size_t i = 0; i < numSlices; ++i)
    {
        size_t sliceStart = i * elementsPerSlice;
        size_t sliceEnd = std::min(sliceStart + elementsPerSlice, numElements);

        if (sliceStart >= sliceEnd) continue; // Should not happen

        std::span<SpatialEntry2D> slice = currentElements.subspan(sliceStart, sliceEnd - sliceStart);
        buildRecursive(slice, 1 - axis, firstChild + i);
    }

    BBox2D bounds = _nodes[firstChild].bounds;
    for (size_t i = 1; i < numSlices; ++i)
        bounds.encapsulate(_nodes[firstChild + i].bounds);
    _nodes[nodeIndex].bounds = bounds;
}

// --- Queries ---

void RTree2D::query_range(const BBox2D& query, 
    const ShapeCallback& callback, const ShapeFilter& filter) const
{
    if (!callback)
        return;
    if (_rootIndex == InvalidNode)
    {
        // Tree is empty
        return;
    }
    if (!_nodes[_rootIndex].bounds.intersects(query))
        return;

    rangeSearchRecursive(_rootIndex, query, callback, filter);
}

bool RTree2D::rangeSearchRecursive(size_t nodeIndex, const BBox2D& query,
    const ShapeCallback& callback, const ShapeFilter& filter) const
{
    const Node& node = _nodes[nodeIndex];

    if (node.isLeaf)
    {
        // Check elements directly
        for (const auto& entry : leafEntries(node))
        {
            if (entry.bounds.intersects(query) && (!filter || filter(entry.shape)))
                if (!callback(entry.shape))
                    return false;
        }
    }
    else
    {
        // Recursive search children
        for (size_t i = 0; i < node.count; ++i)
        {
            const size_t childIndex = node.firstIndex + i;
            if (_nodes[childIndex].bounds.intersects(query))
            {
                if (!rangeSearchRecursive(childIndex, query, callback, filter))
                    return false;
            }
        }
    }

    return true;
}

bool RTree2D::pointSearchRecursive(size_t nodeIndex, const Vector2D &point,
    const ShapeCallback &callback, const ShapeFilter &filter) const
{
    const Node &node = _nodes[nodeIndex];

    if (node.isLeaf)
    {
        for (const auto& entry : leafEntries(node))
        {
            if (entry.bounds.contains(point) && (!filter || filter(entry.shape)))
                if (!callback(entry.shape))
                    return false;
        }
        return true;
    }

    for (size_t i = 0; i < node.count; ++i)
    {
        const size_t childIndex = node.firstIndex + i;
        if (_nodes[childIndex].bounds.contains(point))
            if (!pointSearchRecursive(childIndex, point, callback, filter))
                return false;
    }

    return true;
}

ShapeID RTree2D::query_nearest(const Vector2D& queryPoint, const ShapeFilter& filter) const
{
    if (_rootIndex == InvalidNode)
        return ShapeID::invalid();

    // Min-heap for nodes and elements to visit
    std::priority_queue<std::pair<real_t, size_t>,
        std::vector<std::pair<real_t, size_t>>,
        std::greater<std::pair<real_t, size_t>>> nodesToVisitQueue;

    nodesToVisitQueue.push({ _nodes[_rootIndex].bounds.minDistanceSquared(queryPoint), _rootIndex });

    ShapeID nearestShape = ShapeID::invalid();
    real_t nearestDistanceSq = std::numeric_limits<real_t>::max();

    while (!nodesToVisitQueue.empty())
    {
        const auto current = nodesToVisitQueue.top();
        nodesToVisitQueue.pop();

        const real_t currentMinDistSq = current.first;
        const size_t currentIndex = current.second;
        if (currentMinDistSq >= nearestDistanceSq)
            break;

        const Node& node = _nodes[currentIndex];
        if (node.isLeaf)
        {
            for (const auto& entry : leafEntries(node))
            {
                const Vector2D elementPoint{ entry.bounds.m_min.x, entry.bounds.m_min.y };
                const real_t actualDistSq = elementPoint.distanceSquared(queryPoint);
                if ((!filter || filter(entry.shape)) && actualDistSq < nearestDistanceSq)
                {
                    nearestDistanceSq = actualDistSq;
                    nearestShape = entry.shape;
                }
            }
        }
        else
        {
            for (size_t i = 0; i < node.count; ++i)
            {
                const size_t childIndex = node.firstIndex + i;
                const real_t childMinDistSq = _nodes[childIndex].bounds.minDistanceSquared(queryPoint);
                if (childMinDistSq < nearestDistanceSq)
                    nodesToVisitQueue.push({ childMinDistSq, childIndex });
            }
        }
    }

    return nearestShape;
}

void RTree2D::query_knn(const Vector2D& queryPoint, size_t k,
    const ShapeCallback& callback, const ShapeFilter& filter) const
{
    if (!callback || k == 0 || _rootIndex == InvalidNode)
        return;

    std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare> queue;
    queue.push({ _nodes[_rootIndex].bounds.minDistanceSquared(queryPoint), false, _rootIndex });

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
        if (node.isLeaf)
        {
            for (size_t i = 0; i < node.count; ++i)
            {
                const size_t entryIndex = node.firstIndex + i;
                const SpatialEntry2D& entry = _entries[entryIndex];
                if (!filter || filter(entry.shape))
                {
                    queue.push({ entry.bounds.minDistanceSquared(queryPoint), true,
                        static_cast<NodeIndex>(entryIndex) });
                }
            }
        }
        else
        {
            for (size_t i = 0; i < node.count; ++i)
            {
                const size_t childIndex = node.firstIndex + i;
                queue.push({ _nodes[childIndex].bounds.minDistanceSquared(queryPoint),
                    false, static_cast<NodeIndex>(childIndex) });
            }
        }
    }
}

// --- Queries for Point and Ray ---

// Point query
void RTree2D::query_point(const Vector2D& point, const ShapeCallback& callback,
    const ShapeFilter& filter) const
{
    if (!callback || _rootIndex == InvalidNode || !_nodes[_rootIndex].bounds.contains(point))
        return;

    pointSearchRecursive(_rootIndex, point, callback, filter);
}

// Ray query
void RTree2D::query_ray(const Ray2D& ray, real_t t_min, real_t t_max, 
    const RayHitCallback& callback, const ShapeFilter& filter) const
{
    if (!callback || _rootIndex == InvalidNode)
        return;

    real_t rootEnter;
    real_t rootExit;
    if (!detail::rayBoundsInterval(ray, _nodes[_rootIndex].bounds, t_min, t_max, rootEnter, rootExit))
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
        if (node.isLeaf)
        {
            for (NodeIndex i = 0; i < node.count; ++i)
            {
                const NodeIndex entryIndex = node.firstIndex + i;
                const SpatialEntry2D& entry = _entries[entryIndex];

                if (filter && !filter(entry.shape))
                    continue;

                HitInfo2D hitInfo;
                if (raycast(ray, entry.bounds, t_min, t_max, &hitInfo))
                    queue.push({ hitInfo.t, true, entryIndex });
            }
        }
        else
        {
            for (size_t i = 0; i < node.count; ++i)
            {
                pushRayNode(queue, ray, t_min, t_max, node.firstIndex + i);
            }
        }
        
    }
}

// --- Misc ---

RTreeFilter RTree2D::createFilter(const BBox2D& region, size_t maxDepth) const
{
    RTreeFilter filter;
    filter.region = region;

    if (_rootIndex == InvalidNode)
        return filter;

    collectFilterNodes(_rootIndex, region, 0, maxDepth, filter.nodeIndices, filter.bounds);
    return filter;
}

void RTree2D::collectFilterNodes(size_t nodeIndex,
    const BBox2D& region,
    size_t depth,
    size_t maxDepth,
    std::vector<size_t>& outNodes,
    std::vector<BBox2D>& outBounds) const
{
    const Node& node = _nodes[nodeIndex];
    BBox2D nodeMbr = calculateNodeMBR(nodeIndex);

    // Check if still relevant
    if (!nodeMbr.intersects(region))
        return;

    // Target Depth or Leaf
    if (depth >= maxDepth || node.isLeaf)
    {
        outNodes.push_back(nodeIndex);
        outBounds.push_back(nodeMbr);
        return;
    }

    // Otherwise, go deeper
    for (size_t i = 0; i < node.count; ++i)
    {
        const size_t childIndex = node.firstIndex + i;
        if (_nodes[childIndex].bounds.intersects(region))
        {
            collectFilterNodes(childIndex, region, depth + 1, maxDepth, outNodes, outBounds);
        }
    }
}



} // namespace Math

} // namespace Utility