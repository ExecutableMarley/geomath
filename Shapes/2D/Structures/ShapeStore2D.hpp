#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "../IShape2D.hpp"
#include "../BBox2D.hpp"

namespace Arns
{

namespace Math
{

//Prototype
class ShapeStore2D
{
public:
    size_t add(std::unique_ptr<IFiniteShape2D> shape)
    {
        size_t idx;
        if (!_freeIndices.empty())
        {
            idx = _freeIndices.back();
            _freeIndices.pop_back();
            _shapes[idx] = std::move(shape);
        }
        else
        {
            idx = _shapes.size();
            _shapes.push_back(std::move(shape));
        }
        return idx;
    }

    size_t add(const IFiniteShape2D& shape)
    {
        return add(shape.clone());
    }
    
    template <typename Container>
    std::vector<size_t> add(const Container& shapes)
    {
        std::vector<size_t> indices;
        indices.reserve(shapes.size());
        
        for (const auto& shape_ref : shapes)
        {
            indices.push_back(add(shape_ref.clone()));
        }
        
        return indices;
    }

    void remove(size_t index)
    {
        if (index < _shapes.size() && _shapes[index])
        {
            _shapes[index].reset();
            _freeIndices.push_back(index);
        }
    }

    void remove(const std::vector<size_t>& indices)
    {
        for (auto index : indices)
            remove(index);
    }

    void clear()
    {
        _shapes.clear();
        _freeIndices.clear();
    }

    IFiniteShape2D* get(size_t index) const
    {
        return (index < _shapes.size()) ? _shapes[index].get() : nullptr;
    }

    size_t count() const { return _shapes.size() - _freeIndices.size(); }
    size_t size()  const { return _shapes.size(); }
    //size_t maxIndex() const { return size(); }
    //Todo: Consider sorted active indices list

    bool isValid(size_t index) const
    {
        return index < _shapes.size() && _shapes[index] != nullptr;
    }

    //Todo: Additional Iteration methods

    void forEach(const std::function<void(size_t, IFiniteShape2D&)>& fn)
    {
        for (size_t i = 0; i < _shapes.size(); ++i)
        if (_shapes[i])
            fn(i, *_shapes[i]);
    }

    void forEachConst(const std::function<void(size_t, const IFiniteShape2D&)>& fn) const
    {
        for (size_t i = 0; i < _shapes.size(); ++i)
            if (_shapes[i])
                fn(i, *_shapes[i]);
    }

    //Used to construct spatial index structures
    void getAllShapeBoundsPairs(std::vector<std::pair<BBox2D, size_t>>& out) const
    {
        out.clear();
        out.reserve(size());
        for (size_t i = 0; i < _shapes.size(); i++)
        {
            if (_shapes[i])
                out.emplace_back(_shapes[i]->boundingBox(), i);
        }
    }

private:
    std::vector<std::unique_ptr<IFiniteShape2D>> _shapes;
    std::vector<size_t> _freeIndices;

    //Todo: cache bounding boxes?
    //Todo: cache combined bounding box?
};

} // namespace Math

} // namespace Arns