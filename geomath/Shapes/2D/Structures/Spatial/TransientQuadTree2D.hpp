/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "../../BBox2D.hpp"
#include "../../Algorithms2D/Algorithms2D.hpp"
#include "../ISpatialIndex2D.hpp"

#include <algorithm>
#include <initializer_list>
#include <limits>
#include <queue>
#include <utility>
#include <array>
#include <vector>
#include <span>

namespace Arns::geomath
{

using NodeIndex = uint32_t;

static constexpr NodeIndex INVALID_NODE = std::numeric_limits<NodeIndex>::max();

/*
struct QuadTreeEntry
{
    ShapeID shape;
    BBox2D bounds;
};
*/

struct FastQuadTreeNode
{
    NodeIndex firstChild = INVALID_NODE; // Index of Top-Left child; next 3 are +1, +2, +3
    uint32_t firstEntry = 0;             // Index into global QuadTreeEntry vector
    uint32_t entryCount = 0;

    bool isLeaf() const noexcept { return firstChild == INVALID_NODE; }
};


class TransientQuadTree2D : public ISpatialIndex2D
{
public:
    TransientQuadTree2D(const BBox2D& bounds, size_t maxObjects, size_t maxDepth)
        : _rootBounds(bounds),
          _maxEntries(maxObjects),
          _maxDepth(maxDepth),
          _rootIndex(0)
    {
        //_nodes.emplace_back(QuadTreeNode(bounds));
    }

    // Dynamic insertion
    //void insert(ShapeID shapeIndex, const BBox2D& bounds) override {}

    using ISpatialIndex2D::build;

    // Bulk-loading method
    void build(const std::span<SpatialEntry> entries)
    {
        this->clear();

        _entries.reserve(entries.size());
        _nodes.reserve(entries.size());

        _nodes.emplace_back(FastQuadTreeNode{});
        _rootIndex = 0;

        // 2. Partition out-of-bounds objects upfront
        const BBox2D& rootBBox = _rootBounds;
        auto oobStart = std::partition(entries.begin(), entries.end(), [&rootBBox](const SpatialEntry& e) {
        return e.bounds.m_min.x >= rootBBox.m_min.x && e.bounds.m_max.x <= rootBBox.m_max.x &&
               e.bounds.m_min.y >= rootBBox.m_min.y && e.bounds.m_max.y <= rootBBox.m_max.y;
        });

        std::span<SpatialEntry> validEntries(entries.begin(), oobStart);
        std::span<SpatialEntry> oobEntries(oobStart, entries.end());

        //Store out of bounds entries separately
        this->firstEntry = static_cast<uint32_t>(_entries.size());
        this->entryCount = static_cast<uint32_t>(oobEntries.size());
        if (!oobEntries.empty())
        {
            _entries.insert(_entries.end(), oobEntries.begin(), oobEntries.end());
        }

        if (!validEntries.empty())
        {
            buildNodeTopDown(0, _rootBounds, validEntries, 0);
        }
    }

    void buildNodeTopDown(NodeIndex nodeIdx, const BBox2D& bounds, 
                         std::span<SpatialEntry> entries, size_t depth)
    {
        if (entries.size() <= _maxEntries || depth >= _maxDepth)
        {
            _nodes[nodeIdx].firstEntry = static_cast<uint32_t>(_entries.size());
            _nodes[nodeIdx].entryCount = static_cast<uint32_t>(entries.size());
            _entries.insert(_entries.end(), entries.begin(), entries.end());
            return;
        }

        // 1. Partition entries into 5 buckets: Stays in Parent, NW, NE, SW, SE
        Vector2D c = bounds.centroid();

        auto parentEnd = std::partition(entries.begin(), entries.end(), [c](const SpatialEntry& e) {
            bool straddlesX = (e.bounds.m_min.x < c.x) && (e.bounds.m_max.x > c.x);
            bool straddlesY = (e.bounds.m_min.y < c.y) && (e.bounds.m_max.y > c.y);
            return straddlesX || straddlesY;
        });

        // Save parent-level items to the global flat vector
        _nodes[nodeIdx].firstEntry = static_cast<uint32_t>(_entries.size());
        _nodes[nodeIdx].entryCount = static_cast<uint32_t>(std::distance(entries.begin(), parentEnd));
        _entries.insert(_entries.end(), entries.begin(), parentEnd);

        std::span<SpatialEntry> childEntries(parentEnd, entries.end());
        if (childEntries.empty()) {
            return; // Nothing in the childs -> Early return
        }

        auto southStart = std::partition(childEntries.begin(), childEntries.end(), [c](const SpatialEntry& e) {
            return e.bounds.m_min.y >= c.y; // Top / North halves
        });

        std::span<SpatialEntry> northEntries(childEntries.begin(), southStart);
        std::span<SpatialEntry> southEntries(southStart, childEntries.end());

        auto neStart = std::partition(northEntries.begin(), northEntries.end(), [c](const SpatialEntry& e) {
            return e.bounds.m_max.x <= c.x; // NW: completely left of center
        });

        auto seStart = std::partition(southEntries.begin(), southEntries.end(), [c](const SpatialEntry& e) {
            return e.bounds.m_max.x <= c.x; // SW: completely left of center
        });

        std::array<std::span<SpatialEntry>, 4> quadrants = {
            std::span<SpatialEntry>(northEntries.begin(), neStart), // 0: NW
            std::span<SpatialEntry>(neStart, northEntries.end()),   // 1: NE
            std::span<SpatialEntry>(southEntries.begin(), seStart), // 2: SW
            std::span<SpatialEntry>(seStart, southEntries.end())    // 3: SE
        };

        // 2. Allocate 4 children contiguously in the arena vector
        uint32_t firstChildIdx = static_cast<uint32_t>(_nodes.size());

        _nodes[nodeIdx].firstChild = firstChildIdx;
        for (int i = 0; i < 4; ++i)
        {
            _nodes.emplace_back(FastQuadTreeNode{});
        }

        // 3. Recurse into children with sub-spans
        for (int i = 0; i < 4; ++i) 
        {
            BBox2D childBounds = getChildBounds(bounds, i);
            buildNodeTopDown(firstChildIdx + i, childBounds, quadrants[i], depth + 1);
        }

    }

    // Reset method
    void clear()
    {
        _nodes.clear();
        _entries.clear();
        _rootIndex = INVALID_NODE;
        firstEntry = 0;
        entryCount = 0;
    }

    //

    size_t size() const
    {
        return _entries.size();
    }


    //--- [Query] ---

    // Range search
    void query_range(const BBox2D& query, 
        const ShapeCallback& callback, const ShapeFilter& filter = {}) const override
    {
        if (!callback)
            return;

        //Currently firstEntry and rootIndex will always be 0.

        //Check out of bounds overflow objects
        for (uint32_t i = 0; i < entryCount; i++)
        {
            const SpatialEntry& entry = _entries[firstEntry + i];
            if (filter && !filter(entry.shape))
                    continue;

            if (query.intersects(entry.bounds))
                if (!callback(entry.shape))
                    return;
        }

        //Check root bounds
        if (!query.intersects(_rootBounds) || _rootIndex == INVALID_NODE || _nodes.empty()) return;

        //recursive call
        rangeSearchRecursive(_rootIndex, query, _rootBounds, callback, filter);
    }

    // k-nearest Neighbors with k = 1
    ShapeID query_nearest(const Vector2D& queryPoint, const ShapeFilter& filter = {}) const override
    {
        CandidateQueue queue;

        for (uint32_t i = 0; i < entryCount; ++i)
            pushEntry(queue, firstEntry + i, queryPoint, filter);

        if (!_nodes.empty() && _rootIndex != INVALID_NODE)
            pushNode(queue, _rootIndex, _rootBounds, queryPoint);

        while (!queue.empty())
        {
            const Candidate candidate = queue.top();
            queue.pop();

            if (candidate.isEntry)
                return _entries[candidate.index].shape;

            const FastQuadTreeNode& node = _nodes[candidate.index];

            for (uint32_t i = 0; i < node.entryCount; ++i)
                pushEntry(queue, node.firstEntry + i, queryPoint, filter);

            pushChildren(queue, node, candidate.bounds, queryPoint);
        }

        return ShapeID::invalid();
    }

    // k-Nearest Neighbors
    void query_knn(const Vector2D &queryPoint, size_t k, const ShapeCallback &callback,
                   const ShapeFilter &filter = {}) const override
    {
        if (!callback || k == 0)
            return;

        CandidateQueue queue;

        for (uint32_t i = 0; i < entryCount; ++i)
            pushEntry(queue, firstEntry + i, queryPoint, filter);

        if (!_nodes.empty() && _rootIndex != INVALID_NODE)
            pushNode(queue, _rootIndex, _rootBounds, queryPoint);

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

            const FastQuadTreeNode& node = _nodes[candidate.index];

            for (uint32_t i = 0; i < node.entryCount; ++i)
                pushEntry(queue, node.firstEntry + i, queryPoint, filter);

            pushChildren(queue, node, candidate.bounds, queryPoint);
        }
    }

    // Point query
    void query_point(const Vector2D& point, const ShapeCallback& callback,
        const ShapeFilter& filter = {}) const override
    {
        if (!callback)
            return;

        //Currently firstEntry and rootIndex will always be 0.

        //Check out of bounds overflow objects
        for (uint32_t i = 0; i < entryCount; i++)
        {
            const SpatialEntry& entry = _entries[firstEntry + i];
            if (filter && !filter(entry.shape))
                continue;

            if (entry.bounds.contains(point))
                if (!callback(entry.shape))
                    return;
        }

        //Check root bounds
        if (!_rootBounds.contains(point) || _rootIndex == INVALID_NODE || _nodes.empty()) return;

        pointSearchRecursive(_rootIndex, point, _rootBounds, callback, filter);
    }

    // Ray query
    void query_ray(const Ray2D& ray, real_t t_min, real_t t_max,
        const RayHitCallback& callback, const ShapeFilter& filter = {}) const override
    {
        if (!callback)
            return;

        CandidateQueue queue;

        for (uint32_t i = 0; i < entryCount; ++i)
        {
            pushRayEntry(queue, firstEntry + i, ray, t_min, t_max, filter);
        }

        if (!_nodes.empty() && _rootIndex != INVALID_NODE)
            pushRayNode(queue, _rootIndex, _rootBounds, ray, t_min, t_max);

        while (!queue.empty())
        {
            const Candidate candidate = queue.top();
            queue.pop();

            if (candidate.isEntry)
            {
                if (!callback(_entries[candidate.index].shape, candidate.priority))
                {
                    return;
                }

                continue;
            }

            const FastQuadTreeNode& node = _nodes[candidate.index];

            for (uint32_t i = 0; i < node.entryCount; ++i)
            {
                pushRayEntry(queue, node.firstEntry + i, ray, t_min, t_max, filter);
            }

            if (node.firstChild == INVALID_NODE)
                continue;

            for (uint32_t i = 0; i < 4; ++i)
            {
                const BBox2D childBounds =
                    getChildBounds(candidate.bounds, static_cast<int>(i));

                pushRayNode(queue, node.firstChild + i, childBounds, ray, t_min, t_max);
            }
        }
    }

    void query_ray(const Ray2D& ray, real_t t_min, real_t t_max, const IShapeResolver2D& resolver,
        const RayHitCallback& callback, const ShapeFilter& filter = {}) const 
    {

    }

    void query_intersecting_pairs(const PairCallback &callback) const
    {
        if (!callback)
            return;

        // overflow/out-of-bounds entries against each other
        for (uint32_t i = 0; i < entryCount; ++i)
        {
            const auto &entryA = _entries[firstEntry + i];
            for (uint32_t j = i + 1; j < entryCount; ++j)
            {
                const auto &entryB = _entries[firstEntry + j];
                if (entryA.bounds.intersects(entryB.bounds))
                {
                    if (!callback(entryA.shape, entryB.shape))
                        return;
                }
            }
        }

        // overflow entries against the main tree
        if (_rootIndex != INVALID_NODE && !_nodes.empty())
        {
            for (uint32_t i = 0; i < entryCount; ++i)
            {
                const auto &overflowEntry = _entries[firstEntry + i];
                const auto idA = overflowEntry.shape;
                query_range(overflowEntry.bounds, [&](ShapeID idB)
                {
                    /*
                    if (overflowEntry.shape.index < idB.index)
                    {
                        return callback(overflowEntry.shape, idB);
                    }
                    return true;
                    */
                    const bool reverse = idA.index > idB.index;
                    const auto p1 = reverse ? idB : idA;
                    const auto p2 = reverse ? idA : idB;

                    if (!callback(p1, p2))
                        return true;
                });
            }

            // Run full recursive self-collision on the tree
            selfCollideRecursive(_rootIndex, _rootBounds, callback);
        }
    }

    struct DebugNode
    {
        BBox2D bounds;
        size_t depth;
        uint32_t entryCount;
        bool isLeaf;
    };

    template <typename Visitor>
    void visitDebugNodes(Visitor&& visitor) const
    {
        if (_rootIndex == INVALID_NODE || _nodes.empty())
            return;

        visitDebugNodesRecursive(_rootIndex, _rootBounds, 0, visitor);
    }

private:
    struct Candidate
    {
        real_t priority;
        bool isEntry;
        NodeIndex index;
        BBox2D bounds{}; // Used only for node candidates
    };

    struct CandidateCompare
    {
        bool operator()(const Candidate& left, const Candidate& right) const noexcept
        {
            return left.priority > right.priority;
        }
    };

    using CandidateQueue = std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare>;

    BBox2D getChildBounds(const BBox2D& parent, int childIndex) const noexcept 
    {
        Vector2D c = parent.centroid();
        switch (childIndex) {
            case 0: return { {parent.m_min.x, c.y}, {c.x, parent.m_max.y} }; // 0: NW
            case 1: return { c, parent.m_max };                              // 1: NE
            case 2: return { parent.m_min, c };                              // 2: SW
            case 3: return { {c.x, parent.m_min.y}, {parent.m_max.x, c.y} }; // 3: SE
            default: return parent;
        }
    }

    void pushEntry(CandidateQueue& queue, NodeIndex entryIndex, const Vector2D& point, const ShapeFilter& filter) const
    {
        const SpatialEntry& entry = _entries[entryIndex];
        if (!filter || filter(entry.shape))
            queue.push({entry.bounds.minDistanceSquared(point), true, entryIndex, {}});
    }

    void pushNode(CandidateQueue& queue, NodeIndex nodeIndex, const BBox2D& bounds, const Vector2D& point) const
    {
        if (nodeIndex == INVALID_NODE)
            return;

        queue.push({bounds.minDistanceSquared(point), false, nodeIndex, bounds});
    }

    void pushRayEntry(CandidateQueue& queue, NodeIndex entryIndex, const Ray2D& ray,
        real_t tMin, real_t tMax, const ShapeFilter& filter) const
    {
        const SpatialEntry& entry = _entries[entryIndex];

        if (filter && !filter(entry.shape))
            return;

        HitInfo2D hitInfo;
        if (intersect(ray, entry.bounds, tMin, tMax, &hitInfo))
        {
            queue.push({hitInfo.t, true, entryIndex, {}});
        }
    }

    void pushRayNode(CandidateQueue& queue, NodeIndex nodeIndex, const BBox2D& bounds,
        const Ray2D& ray, real_t tMin, real_t tMax) const
    {
        if (nodeIndex == INVALID_NODE)
            return;
        real_t enter;
        real_t exit;
        if (detail::rayBoundsInterval(ray, bounds, tMin, tMax, enter, exit))
            queue.push({enter, false, nodeIndex, bounds});
    }

    void pushChildren(CandidateQueue& queue, const FastQuadTreeNode& node, const BBox2D& bounds,
        const Vector2D& point) const
    {
        if (node.firstChild == INVALID_NODE)
            return;

        for (uint32_t i = 0; i < 4; ++i)
        {
            const BBox2D childBounds = getChildBounds(bounds, static_cast<int>(i));
            pushNode(queue, node.firstChild + i, childBounds, point);
        }
    }

    // Recursive helper for range search
    bool rangeSearchRecursive(NodeIndex nodeIdx, const BBox2D& query, const BBox2D& curBounds,
        const ShapeCallback& callback, const ShapeFilter& filter) const
    {
        //Check this nodes entries
        for (NodeIndex i = 0; i < _nodes[nodeIdx].entryCount; i++)
        {
            const SpatialEntry& entry = _entries[_nodes[nodeIdx].firstEntry + i];
            if (filter && !filter(entry.shape))
                    continue;

            if (query.intersects(entry.bounds))
                if (!callback(entry.shape))
                    return false;
        }

        //Check children
        NodeIndex firstChild = _nodes[nodeIdx].firstChild;
        if (firstChild != INVALID_NODE)
        {
            for (uint32_t i = 0; i < 4; i++)
            {
                BBox2D childBounds = getChildBounds(curBounds, i);
                if (query.intersects(childBounds))
                    if (!rangeSearchRecursive(firstChild + i, query, childBounds, callback, filter))
                        return false;
            }
        }

        return true;
    }

    bool pointSearchRecursive(NodeIndex nodeIdx, const Vector2D& point, const BBox2D& curBounds,
        const ShapeCallback& callback, const ShapeFilter& filter) const
    {
        //Check this nodes entries
        for (NodeIndex i = 0; i < _nodes[nodeIdx].entryCount; i++)
        {
            const SpatialEntry& entry = _entries[_nodes[nodeIdx].firstEntry + i];
            if (filter && !filter(entry.shape))
                continue;

            if (entry.bounds.contains(point))
                if (!callback(entry.shape))
                    return false;
        }

        //Check children
        NodeIndex firstChild = _nodes[nodeIdx].firstChild;
        if (firstChild != INVALID_NODE)
        {
            for (uint32_t i = 0; i < 4; i++)
            {
                BBox2D childBounds = getChildBounds(curBounds, i);
                if (childBounds.contains(point))
                    if (!pointSearchRecursive(firstChild + i, point, childBounds, callback, filter))
                        return false;
            }
        }

        return true;
    }

    // Recursive pair intersection helpers

    void collideNodeWithNode(NodeIndex nodeAIdx, NodeIndex nodeBIdx, const PairCallback &callback) const
    {
        const auto &nodeA = _nodes[nodeAIdx];
        const auto &nodeB = _nodes[nodeBIdx];

        for (uint32_t i = 0; i < nodeA.entryCount; ++i)
        {
            const auto &entryA = _entries[nodeA.firstEntry + i];
            for (uint32_t j = 0; j < nodeB.entryCount; ++j)
            {
                const auto &entryB = _entries[nodeB.firstEntry + j];

                // Skip self-pairs
                if (entryA.shape.index == entryB.shape.index)
                    continue;

                if (entryA.bounds.intersects(entryB.bounds))
                {
                    // Apparently im stupid
                    /*
                    if (entryA.shape.index < entryB.shape.index)
                    {
                        if (!callback(entryA.shape, entryB.shape))
                            return;
                    }
                    */
                    const bool reverse = entryA.shape.index > entryB.shape.index;
                    const auto p1 = reverse ? entryB.shape : entryA.shape;
                    const auto p2 = reverse ? entryA.shape : entryB.shape;

                    if (!callback(p1, p2))
                        return;
                }
            }
        }
    }

    void collideNodeWithSubtree(NodeIndex nodeAIdx, NodeIndex nodeBIdx, const PairCallback &callback) const
    {
        if (nodeBIdx == INVALID_NODE)
            return;

        // Check nodeA entries against nodeB entries
        collideNodeWithNode(nodeAIdx, nodeBIdx, callback);

        // Recurse down nodeB's children
        NodeIndex firstChild = _nodes[nodeBIdx].firstChild;
        if (firstChild != INVALID_NODE)
        {
            for (uint32_t i = 0; i < 4; ++i)
            {
                collideNodeWithSubtree(nodeAIdx, firstChild + i, callback);
            }
        }
    }

    void collideSubtrees(NodeIndex nodeAIdx, NodeIndex nodeBIdx, const PairCallback &callback) const
    {
        if (nodeAIdx == INVALID_NODE || nodeBIdx == INVALID_NODE)
            return;

        // Cross entries of A with B and B's descendants
        collideNodeWithSubtree(nodeAIdx, nodeBIdx, callback);

        // Recurse down A's children against B
        NodeIndex firstChildA = _nodes[nodeAIdx].firstChild;
        if (firstChildA != INVALID_NODE)
        {
            for (uint32_t i = 0; i < 4; ++i)
            {
                collideSubtrees(firstChildA + i, nodeBIdx, callback);
            }
        }
    }

    void selfCollideRecursive(NodeIndex nodeIdx, const BBox2D &curBounds, const PairCallback &callback) const
    {
        if (nodeIdx == INVALID_NODE)
            return;
        const auto &node = _nodes[nodeIdx];

        // A. Collide entries within this node against each other
        for (uint32_t i = 0; i < node.entryCount; ++i)
        {
            const auto &entryA = _entries[node.firstEntry + i];
            for (uint32_t j = i + 1; j < node.entryCount; ++j)
            {
                const auto &entryB = _entries[node.firstEntry + j];
                if (entryA.bounds.intersects(entryB.bounds))
                {
                    if (!callback(entryA.shape, entryB.shape))
                        return;
                }
            }
        }

        NodeIndex firstChild = _nodes[nodeIdx].firstChild;
        if (firstChild != INVALID_NODE)
        {
            // B. Collide entries in this node against all entries in its descendant subtrees
            for (uint32_t i = 0; i < 4; ++i)
            {
                collideNodeWithSubtree(nodeIdx, firstChild + i, callback);
            }

            // C. Recurse down children for internal self-collisions
            for (uint32_t i = 0; i < 4; ++i)
            {
                selfCollideRecursive(firstChild + i, getChildBounds(curBounds, i), callback);
            }

            // D. Collide distinct child subtrees against each other if their bounds intersect
            for (uint32_t i = 0; i < 4; ++i)
            {
                BBox2D boundsI = getChildBounds(curBounds, i);
                for (uint32_t j = i + 1; j < 4; ++j)
                {
                    BBox2D boundsJ = getChildBounds(curBounds, j);
                    if (boundsI.intersects(boundsJ))
                    {
                        collideSubtrees(firstChild + i, firstChild + j, callback);
                    }
                }
            }
        }
    }

    template <typename Visitor>
    void visitDebugNodesRecursive(NodeIndex nodeIndex, const BBox2D& bounds,
        size_t depth, Visitor& visitor) const
    {
        const FastQuadTreeNode& node = _nodes[nodeIndex];
        visitor(DebugNode{bounds, depth, node.entryCount, node.isLeaf()});

        if (node.isLeaf())
            return;

        for (uint32_t child = 0; child < 4; ++child)
        {
            const BBox2D childBounds = getChildBounds(bounds, child);
            visitDebugNodesRecursive(
                node.firstChild + child,
                childBounds,
                depth + 1,
                visitor);
        }
    }

private:
    BBox2D _rootBounds;
    size_t _maxEntries;
    size_t _maxDepth;
    NodeIndex _rootIndex = INVALID_NODE;
    NodeIndex firstEntry = 0;
    NodeIndex entryCount = 0;
    std::vector<FastQuadTreeNode> _nodes;
    std::vector<SpatialEntry> _entries;
};











}