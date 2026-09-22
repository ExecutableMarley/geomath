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
#include <functional>
#include <limits>
#include <queue>
#include <vector>

namespace Arns::geomath
{


class BVH2D : public ISpatialIndex2D
{
public:
	explicit BVH2D(size_t maxLeafEntries = 4)
		: _maxLeafEntries(std::max<size_t>(1, maxLeafEntries)) {}

	using ISpatialIndex2D::build;

	void build(const std::span<SpatialEntry> entries) override
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

	//Todo: Somewhat remove overlap with RTree code

	void query_range(const BBox2D& query, const ShapeCallback& callback,
		const ShapeFilter& filter = {}) const override
	{
		if (callback && _rootIndex != InvalidNode)
			rangeSearchRecursive(_rootIndex, query, callback, filter);
	}

	ShapeID query_nearest(const Vector2D& point, const ShapeFilter& filter = {}) const override
	{
		if (_rootIndex == InvalidNode)
			return ShapeID::invalid();

		std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare> queue;
		queue.push({ _nodes[_rootIndex].bounds.minDistanceSquared(point), false, _rootIndex });
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

	void query_knn(const Vector2D& point, size_t k, const ShapeCallback& callback,
		const ShapeFilter& filter = {}) const override
	{
		if (!callback || k == 0 || _rootIndex == InvalidNode)
			return;

		std::priority_queue<Candidate, std::vector<Candidate>, CandidateCompare> queue;
		queue.push({ _nodes[_rootIndex].bounds.minDistanceSquared(point), false, _rootIndex });
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
			for (uint32_t i = 0; i < node.count; ++i)
				pushEntry(queue, node.firstEntry + i, point, filter);
			pushNode(queue, node.left, point);
			pushNode(queue, node.right, point);
		}
	}

	void query_point(const Vector2D& point, const ShapeCallback& callback,
		const ShapeFilter& filter = {}) const override
	{
		if (callback && _rootIndex != InvalidNode)
			pointSearchRecursive(_rootIndex, point, callback, filter);
	}

	void query_ray(const Ray2D& ray, real_t t_min, real_t t_max,
		const RayHitCallback& callback, const ShapeFilter& filter = {}) const override
	{
		if (!callback || _rootIndex == InvalidNode)
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
				const SpatialEntry& entry = _entries[node.firstEntry + i];
				if (filter && !filter(entry.shape))
					continue;
				HitInfo2D hitInfo;
				if (intersect(ray, entry.bounds, t_min, t_max, &hitInfo))
					queue.push({ hitInfo.t, true, (node.firstEntry + i) });
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
		BBox2D bounds;
		NodeIndex left = InvalidNode;
		NodeIndex right = InvalidNode;
		uint32_t firstEntry = 0;
		uint32_t count = 0;
	};

	//Possible query time improvements:
	//1. Use a container which can stack allocate k entries (boost::container::small_vector)
	//  -This delays or skips heap allocations depending on query size
	//2. Figure out a way to prune candidates. Does not matter for small queries. 
	//  -Could be done with an max heap for entries. Will however make heap allocations worse
	//   without a better container

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
		BBox2D centroidBounds(_entries[first].bounds.centroid(), _entries[first].bounds.centroid());
		for (size_t i = first + 1; i < last; ++i)
		{
			node.bounds.encapsulate(_entries[i].bounds);
			centroidBounds.encapsulate(_entries[i].bounds.centroid());
		}

		const NodeIndex nodeIndex = static_cast<NodeIndex>(_nodes.size());
		_nodes.push_back(node);
		if (node.count <= _maxLeafEntries)
			return nodeIndex;

		const bool splitX = centroidBounds.width() >= centroidBounds.height();
		const size_t middle = first + node.count / 2;
		std::nth_element(_entries.begin() + first, _entries.begin() + middle, _entries.begin() + last,
			[splitX](const SpatialEntry& left, const SpatialEntry& right) {
				const Vector2D leftCenter = left.bounds.centroid();
				const Vector2D rightCenter = right.bounds.centroid();
				return splitX ? leftCenter.x < rightCenter.x : leftCenter.y < rightCenter.y;
			});

		_nodes[nodeIndex].count = 0;
		_nodes[nodeIndex].left = buildNode(first, middle);
		_nodes[nodeIndex].right = buildNode(middle, last);
		return nodeIndex;
	}

	template <typename Queue>
	void pushEntry(Queue& queue, NodeIndex entryIndex, const Vector2D& point, const ShapeFilter& filter) const
	{
		const SpatialEntry& entry = _entries[entryIndex];
		if (!filter || filter(entry.shape))
			queue.push({ entry.bounds.minDistanceSquared(point), true, static_cast<NodeIndex>(entryIndex) });
	}

	template <typename Queue>
	void pushNode(Queue& queue, NodeIndex nodeIndex, const Vector2D& point) const
	{
		if (nodeIndex != InvalidNode)
			queue.push({ _nodes[nodeIndex].bounds.minDistanceSquared(point), false, nodeIndex });
	}

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

	bool rangeSearchRecursive(NodeIndex nodeIndex, const BBox2D& query,
		const ShapeCallback& callback, const ShapeFilter& filter) const
	{
		const Node& node = _nodes[nodeIndex];
		if (!query.intersects(node.bounds))
			return true;
		for (size_t i = 0; i < node.count; ++i)
		{
			const SpatialEntry& entry = _entries[node.firstEntry + i];
			if ((!filter || filter(entry.shape)) && query.intersects(entry.bounds) && !callback(entry.shape))
				return false;
		}
		return (node.left == InvalidNode || rangeSearchRecursive(node.left, query, callback, filter)) &&
			   (node.right == InvalidNode || rangeSearchRecursive(node.right, query, callback, filter));
	}

	bool pointSearchRecursive(NodeIndex nodeIndex, const Vector2D& point,
		const ShapeCallback& callback, const ShapeFilter& filter) const
	{
		const Node& node = _nodes[nodeIndex];
		if (!node.bounds.contains(point))
			return true;
		for (size_t i = 0; i < node.count; ++i)
		{
			const SpatialEntry& entry = _entries[node.firstEntry + i];
			if ((!filter || filter(entry.shape)) && entry.bounds.contains(point) && !callback(entry.shape))
				return false;
		}
		return (node.left == InvalidNode || pointSearchRecursive(node.left, point, callback, filter)) &&
			   (node.right == InvalidNode || pointSearchRecursive(node.right, point, callback, filter));
	}

private:
	size_t _maxLeafEntries;
	NodeIndex _rootIndex = InvalidNode;
	std::vector<Node> _nodes;
	std::vector<SpatialEntry> _entries;
};



    
}