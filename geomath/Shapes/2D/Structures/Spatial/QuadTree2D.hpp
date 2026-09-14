/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "../../BBox2D.hpp"
#include "../ISpatialIndex2D.hpp"

namespace Arns::geomath
{

using NodeIndex = uint32_t;

static constexpr NodeIndex INVALID_NODE = std::numeric_limits<NodeIndex>::max();

struct QuadTreeEntry
{
    ShapeID shape;
    BBox2D bounds;
};

struct QuadTreeNode
{
    BBox2D bounds;

    std::array<NodeIndex, 4> children{
        INVALID_NODE,
        INVALID_NODE,
        INVALID_NODE,
        INVALID_NODE
    };

    std::vector<QuadTreeEntry> objects;

    QuadTreeNode(BBox2D bounds) : bounds(bounds) {}

    bool isLeaf() const noexcept
    {
        return children[0] == INVALID_NODE;
    }
};


class QuadTree2D : public ISpatialIndex2D
{
public:

    QuadTree2D(const BBox2D& bounds, size_t maxObjects, size_t maxDepth)
        : _rootBounds(bounds),
          _maxEntries(maxObjects),
          _maxDepth(maxDepth),
          _rootIndex(0)
    {
        _nodes.emplace_back(QuadTreeNode(bounds));
    }

    void insert(ShapeID shapeIndex, const BBox2D& bounds) override
    {
        if (!_rootBounds.contains(bounds) || _maxEntries == 0)
        {
            return;
        }

        insertRecursive(_rootIndex, 0, shapeIndex, bounds);
    }

    void insertRecursive(NodeIndex currentNode, uint32_t depth, ShapeID shapeIndex, const BBox2D& bounds)
    {
        const bool isLeaf = _nodes[currentNode].isLeaf();

        if (!isLeaf)
        {
            const QuadTreeNode& node = _nodes[currentNode];
            int index = getIndexOfBounds(node.bounds.centroid(), bounds);
            if (index != -1)
                return insertRecursive(node.children[index], depth + 1, shapeIndex, bounds);
        }

        _nodes[currentNode].objects.emplace_back(QuadTreeEntry{shapeIndex, bounds});

        if (isLeaf && _nodes[currentNode].objects.size() > _maxEntries && depth < _maxDepth)
        {
            subdivide(currentNode);
            std::vector<QuadTreeEntry> entries = std::move(_nodes[currentNode].objects);
            _nodes[currentNode].objects.clear();

            for (const QuadTreeEntry& entry : entries)
            {
                insertRecursive(currentNode, depth, entry.shape, entry.bounds);
            }
        }

    }

    int getIndexOfBounds(const Vector2D& parentMidPoint, const BBox2D& queryBounds)
    {
        bool topQuadrant = (queryBounds.m_min.y < parentMidPoint.y && queryBounds.m_max.y < parentMidPoint.y);
        bool bottomQuadrant = (queryBounds.m_min.y > parentMidPoint.y);
        bool leftQuadrant = (queryBounds.m_min.x < parentMidPoint.x && queryBounds.m_max.x < parentMidPoint.x);
        bool rightQuadrant = (queryBounds.m_min.x > parentMidPoint.x);

        if (topQuadrant && leftQuadrant)
        {
            return 0; // Top-left
        }
        else if (topQuadrant && rightQuadrant)
        {
            return 1; // Top-right
        }
        else if (bottomQuadrant && leftQuadrant)
        {
            return 2; // Bottom-left
        }
        else if (bottomQuadrant && rightQuadrant)
        {
            return 3; // Bottom-right
        }

        return -1; // Does not fit in any child
    }

    void subdivide(NodeIndex currentNode)
    {
        const BBox2D bounds = _nodes[currentNode].bounds;
        real_t subWidth  = bounds.width()  * 0.5f;
        real_t subHeight = bounds.height() * 0.5f;

        real_t x = bounds.m_min.x;
        real_t y = bounds.m_min.y;

        const NodeIndex firstChild = static_cast<NodeIndex>(_nodes.size());

        _nodes.push_back(QuadTreeNode(BBox2D( Vector2D{x, y}, subWidth, subHeight)));
        _nodes.push_back(QuadTreeNode(BBox2D( Vector2D{x + subWidth, y}, subWidth, subHeight)));
        _nodes.push_back(QuadTreeNode(BBox2D( Vector2D{x, y + subHeight}, subWidth, subHeight)));
        _nodes.push_back(QuadTreeNode(BBox2D( Vector2D{x + subWidth, y + subHeight}, subWidth, subHeight)));
    
        _nodes[currentNode].children = {
            firstChild,
            firstChild + 1,
            firstChild + 2,
            firstChild + 3
        };
    }

private:
    BBox2D _rootBounds;
    size_t _maxEntries;
    size_t _maxDepth;
    NodeIndex _rootIndex = INVALID_NODE;
    std::vector<QuadTreeNode> _nodes;
};


}