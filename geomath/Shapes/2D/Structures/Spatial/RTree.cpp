/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#include "RTree.hpp"

#include <queue>


namespace Arns
{

namespace geomath
{

void RTree::build(const std::vector<std::pair<BBox2D, ShapeID>>& elementBoundsWithIndices)
{
    clear();

    //Consider reserving space in _nodes
    if (elementBoundsWithIndices.empty())
    {
        _rootIndex = InvalidNode;
        return;
    }

    std::vector<RTreeLeafEntry> initialEntries;
    initialEntries.reserve(elementBoundsWithIndices.size());
    for (const auto& pair : elementBoundsWithIndices)
    {
        initialEntries.push_back({ pair.first, pair.second });
    }

    _rootIndex = buildRecursive(initialEntries, 0);
    _rootBounds = calculateNodeMBR(_rootIndex);
}

size_t RTree::buildRecursive(std::vector<RTreeLeafEntry>& currentElements, int axis)
{
    // Base Case
    if (currentElements.size() <= _maxEntries)
    {
        size_t newNodeIdx = createNode(true);
        _nodes[newNodeIdx].leafEntries = currentElements;
        return newNodeIdx;
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
            [](const RTreeLeafEntry& a, const RTreeLeafEntry& b)
            {
                return a.bounds.m_min.x < b.bounds.m_min.x;
            });
    }
    else
    {
        // Sort by Y-coordinate
        std::sort(currentElements.begin(), currentElements.end(),
            [](const RTreeLeafEntry& a, const RTreeLeafEntry& b)
            {
                return a.bounds.m_min.y < b.bounds.m_min.y;
            });
    }

    size_t newNodeIdx = createNode(false);

    for (size_t i = 0; i < numSlices; ++i)
    {
        size_t sliceStart = i * elementsPerSlice;
        size_t sliceEnd = std::min(sliceStart + elementsPerSlice, numElements);

        if (sliceStart >= sliceEnd) continue; // Should not happen

        //vector view instead maybe?
        std::vector<RTreeLeafEntry> slice(currentElements.begin() + sliceStart,
            currentElements.begin() + sliceEnd);

        size_t childNodeIndex = buildRecursive(slice, 1 - axis);

        _nodes[newNodeIdx].childEntries.push_back({ calculateNodeMBR(childNodeIndex), childNodeIndex });
    }

    return newNodeIdx;
}

// --- Queries ---

void RTree::query_range(const BBox2D& query, const ShapeCallback& callback,
    const ShapeFilter& filter) const
{
    if (!callback)
        return;
    if (_rootIndex == InvalidNode)
    {
        // Tree is empty
        return;
    }
    if (!_rootBounds.intersects(query))
        return;

    rangeSearchRecursive(_rootIndex, query, callback, filter);
}

void RTree::rangeSearchRecursive(size_t nodeIndex, const BBox2D& query,
    const ShapeCallback& callback, const ShapeFilter& filter) const
{
    const RTreeNode& node = _nodes[nodeIndex];

    if (node.isLeaf)
    {
        // Check elements directly
        for (const auto& entry : node.leafEntries)
        {
            if (entry.bounds.intersects(query) && (!filter || filter(entry.shape)))
                callback(entry.shape);
        }
    }
    else
    {
        // Recursive search children
        for (const auto& entry : node.childEntries)
        {
            if (entry.bounds.intersects(query))
            {
                rangeSearchRecursive(entry.child, query, callback, filter);
            }
        }
    }
}

void RTree::pointSearchRecursive(size_t nodeIndex, const Vector2D& point,
    const ShapeCallback& callback, const ShapeFilter& filter) const
{
    const RTreeNode& node = _nodes[nodeIndex];

    if (node.isLeaf)
    {
        for (const auto& entry : node.leafEntries)
        {
            if (entry.bounds.contains(point) && (!filter || filter(entry.shape)))
                callback(entry.shape);
        }
        return;
    }

    for (const auto& entry : node.childEntries)
    {
        if (entry.bounds.contains(point))
            pointSearchRecursive(entry.child, point, callback, filter);
    }
}

ShapeID RTree::query_nearest(const Vector2D& queryPoint, const ShapeFilter& filter) const
{
    if (_rootIndex == InvalidNode)
        return ShapeID::invalid();

    // Min-heap for nodes and elements to visit
    std::priority_queue<std::pair<real_t, size_t>,
        std::vector<std::pair<real_t, size_t>>,
        std::greater<std::pair<real_t, size_t>>> nodesToVisitQueue;

    BBox2D rootMbr = _rootBounds;
    nodesToVisitQueue.push({ rootMbr.minDistanceSquared(queryPoint), _rootIndex });

    ShapeID nearestShape = ShapeID::invalid();
    real_t nearestDistanceSq = std::numeric_limits<real_t>::max();

    while (!nodesToVisitQueue.empty())
    {
        auto current = nodesToVisitQueue.top();
        nodesToVisitQueue.pop();

        real_t currentMinDistSq = current.first;
        size_t currentIndex = current.second;

        // Pruning
        if (currentMinDistSq >= nearestDistanceSq)
            break;

        const RTreeNode& node = _nodes[currentIndex];

        if (node.isLeaf)
        {
            for (const auto& entry : node.leafEntries)
            {
                Vector2D elementPoint = { entry.bounds.m_min.x, entry.bounds.m_min.y }; // Assuming point data
                real_t actualDistSq = elementPoint.distanceSquared(queryPoint);

                if ((!filter || filter(entry.shape)) && actualDistSq < nearestDistanceSq)
                {
                    nearestDistanceSq = actualDistSq;
                    nearestShape = entry.shape;
                }
            }
        }
        else
        {
            for (const auto& entry : node.childEntries)
            {
                real_t childMinDistSq = entry.bounds.minDistanceSquared(queryPoint);

                // Pruning check
                if (childMinDistSq < nearestDistanceSq)
                    nodesToVisitQueue.push({ childMinDistSq, entry.child });
            }
        }
    }

    return nearestShape;
}

void RTree::query_knn(const Vector2D& queryPoint, size_t k, const ShapeCallback& callback,
    const ShapeFilter& filter) const
{
    if (!callback || k == 0 || _rootIndex == InvalidNode)
    {
        return;
    }

    // Min-heap for nodes and elements to visit
    std::priority_queue<std::pair<real_t, size_t>,
        std::vector<std::pair<real_t, size_t>>,
        std::greater<std::pair<real_t, size_t>>> nodesToVisitQueue;

    // Max-heap for 'k' best results found
    RTreeDistanceQueue bestKResultsQueue;

    BBox2D rootMbr = _rootBounds;
    nodesToVisitQueue.push({ rootMbr.minDistanceSquared(queryPoint), _rootIndex });

    while (!nodesToVisitQueue.empty())
    {
        auto current = nodesToVisitQueue.top();
        nodesToVisitQueue.pop();

        real_t currentMinDistSq = current.first;
        size_t currentIndex = current.second;

        // Pruning
        if (bestKResultsQueue.size() == k && currentMinDistSq > bestKResultsQueue.top().first)
        {
            // All remaining will be further away than bestKResultsQueue.top().first
            break;
        }

        const RTreeNode& node = _nodes[currentIndex];

        if (node.isLeaf)
        {
            for (const auto& entry : node.leafEntries)
            {
                Vector2D elementPoint = { entry.bounds.m_min.x, entry.bounds.m_min.y }; // Assuming point data
                real_t actualDistSq = elementPoint.distanceSquared(queryPoint);

                if ((!filter || filter(entry.shape)) && bestKResultsQueue.size() < k)
                {
                    //Append
                    bestKResultsQueue.push({ actualDistSq, entry.shape });
                }
                else if ((!filter || filter(entry.shape)) && actualDistSq < bestKResultsQueue.top().first)
                {
                    //Replace with closer element
                    bestKResultsQueue.pop();
                    bestKResultsQueue.push({ actualDistSq, entry.shape });
                }
            }
        }
        else
        {
            for (const auto& entry : node.childEntries)
            {
                real_t childMinDistSq = entry.bounds.minDistanceSquared(queryPoint);

                //Pruning check
                if (bestKResultsQueue.size() == k && childMinDistSq > bestKResultsQueue.top().first)
                {
                    continue;
                }
                nodesToVisitQueue.push({ childMinDistSq, entry.child });
            }
        }
    }

    // Extract and reverse results from bestKResultsQueue.
    extractKNearestResults(k, bestKResultsQueue, callback);
}

void RTree::extractKNearestResults(size_t k,
    RTreeDistanceQueue& bestKResultsQueue,
    const ShapeCallback& callback) const
{
    // I wish we could iterate backwards through the priority queue
    // Todo: Consider using a different data structure to avoid this extra step 
    std::vector<std::pair<real_t, ShapeID>> tempResults;
    while (!bestKResultsQueue.empty())
    {
        tempResults.push_back(bestKResultsQueue.top());
        bestKResultsQueue.pop();
    }
    // Reverse
    std::reverse(tempResults.begin(), tempResults.end());

    // Take k results (min check not required atm)
    for (size_t i = 0; i < std::min(k, tempResults.size()); ++i)
        callback(tempResults[i].second);
}

// --- Queries for Point and Ray ---

// Point query
void RTree::query_point(const Vector2D& point, const ShapeCallback& callback,
    const ShapeFilter& filter) const
{
    if (!callback || _rootIndex == InvalidNode || !_rootBounds.contains(point))
        return;

    pointSearchRecursive(_rootIndex, point, callback, filter);
}

// Todo: This will collect all hits. Optimize later or offer a specialized version
void RTree::raySearchRecursive(size_t nodeIndex, const Ray2D& ray, real_t t_min, real_t t_max,
    const ShapeFilter& filter, std::vector<std::pair<real_t, ShapeID>>& hits) const
{
    const RTreeNode& node = _nodes[nodeIndex];

    if (node.isLeaf)
    {
        for (const auto& entry : node.leafEntries)
        {
            if (filter && !filter(entry.shape))
                continue;

            HitInfo2D hitInfo;
            if (intersect(ray, entry.bounds, t_min, t_max, &hitInfo))
                hits.emplace_back(hitInfo.t, entry.shape);
        }
        return;
    }

    for (const auto& entry : node.childEntries)
    {
        HitInfo2D hitInfo;
        if (intersect(ray, entry.bounds, t_min, t_max, &hitInfo))
            raySearchRecursive(entry.child, ray, t_min, t_max, filter, hits);
    }
}

// Ray query
void RTree::query_ray(const Ray2D& ray, real_t t_min, real_t t_max, const RayHitCallback& callback,
    const ShapeFilter& filter) const
{
    if (!callback || _rootIndex == InvalidNode)
        return;

    HitInfo2D rootHit;
    if (!intersect(ray, _rootBounds, t_min, t_max, &rootHit))
        return;

    std::vector<std::pair<real_t, ShapeID>> hits;
    raySearchRecursive(_rootIndex, ray, t_min, t_max, filter, hits);

    std::stable_sort(hits.begin(), hits.end(),
        [](const auto& left, const auto& right)
        {
            return left.first < right.first;
        });

    for (const auto& [hitT, shape] : hits)
    {
        if (!callback(shape, hitT))
            break;
    }
}

// --- Misc ---

RTreeFilter RTree::createFilter(const BBox2D& region, size_t maxDepth) const
{
    RTreeFilter filter;
    filter.region = region;

    if (_rootIndex == InvalidNode)
        return filter;

    collectFilterNodes(_rootIndex, region, 0, maxDepth, filter.nodeIndices, filter.bounds);
    return filter;
}

void RTree::collectFilterNodes(size_t nodeIndex,
    const BBox2D& region,
    size_t depth,
    size_t maxDepth,
    std::vector<size_t>& outNodes,
    std::vector<BBox2D>& outBounds) const
{
    const RTreeNode& node = _nodes[nodeIndex];
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
    for (const auto& entry : node.childEntries)
    {
        if (entry.bounds.intersects(region))
        {
            collectFilterNodes(entry.child, region, depth + 1, maxDepth, outNodes, outBounds);
        }
    }
}



} // namespace Math

} // namespace Utility