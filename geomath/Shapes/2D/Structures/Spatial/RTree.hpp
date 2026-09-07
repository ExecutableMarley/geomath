#include "../../BBox2D.hpp"
#include "../ISpatialIndex2D.hpp"

#include <cassert>
#include <queue>

namespace Arns
{

namespace geomath
{

struct RTreeEntry
{
    BBox2D bounds;
    size_t elementIndex;
};

struct RTreeNode
{
    bool isLeaf = false;
    std::vector<RTreeEntry> entries;
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
    void insert(const BBox2D& elementBound, ShapeID elementIndex) {}

    // Reset method
    void clear()
    {
        this->_rootIndex = InvalidNode;
        this->_nodes.clear();
    }

    // Range search
    void rangeQuery(const BBox2D& query, std::vector<ShapeID>& result) const;

    void rangeQuery(const BBox2D& query, std::vector<ShapeID>& result, const std::vector<bool>& inclusionMask) const;

    void rangeSearchWithFilter(const RTreeFilter& filter, const BBox2D& query,
        std::vector<ShapeID>& result, const std::vector<bool>& inclusionMask) const;

    // k-nearest Neighbours with k = 1
    ShapeID nearestNeighbour(const Vector2D& queryPoint) const;

    ShapeID nearestNeighbour(const Vector2D& queryPoint, const std::vector<bool>& inclusionMask) const;

    // k-Nearest Neighbours
    void kNearest(const Vector2D& queryPoint, size_t k, std::vector<ShapeID>& result) const;

    void kNearest(const Vector2D& queryPoint, size_t k, std::vector<ShapeID>& result, const std::vector<bool>& inclusionMask) const;

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
    size_t buildRecursive(std::vector<RTreeEntry>& currentElements, int axis);

    static BBox2D calculateMBR(const std::vector<RTreeEntry>& entries)
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
        return calculateMBR(_nodes[nodeIndex].entries);
    }

    // Recursive helper for range search
    void rangeSearchRecursive(size_t nodeIndex, const BBox2D& query, std::vector<ShapeID>& result) const;

    void rangeSearchRecursive(size_t nodeIndex, const BBox2D& query, std::vector<ShapeID>& result, const std::vector<bool>& inclusionMask) const;

    // retrieve k final results from priority queue
    void extractKNearestResults(size_t k,
        std::priority_queue<std::pair<double, size_t>>& queue_to_explore,
        std::vector<ShapeID>& result) const;

    void collectFilterNodes(size_t nodeIndex,
        const BBox2D& region,
        size_t depth,
        size_t maxDepth,
        std::vector<size_t>& outNodes,
        std::vector<BBox2D>& outBounds) const;
};

} // namespace Math

} // namespace Arns