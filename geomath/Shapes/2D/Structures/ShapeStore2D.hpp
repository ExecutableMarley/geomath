#pragma once

//#include "../Interfaces/IFiniteShape2D.hpp"
#include "Shapes/Common/SlotID.hpp"
#include "../Variants/ShapeConcepts.hpp"
#include "../Variants/ShapeVariants.hpp"
#include "../BBox2D.hpp"
#include "../Triangle2D.hpp"
#include "../Rectangle2D.hpp"
#include "../Polygon2D.hpp"
#include "../Circle2D.hpp"

#include <memory>
#include <vector>
#include <functional>

namespace Arns
{

namespace geomath
{


struct ShapeHandle2D
{
    ShapeType2D type;
    uint32_t index;
    //uint32_t generation;
};

struct ShapeObject
{
    ShapeHandle2D handle;
    uint32_t generation;
    bool isAlive;
};

class IShapeResolver2D
{
public:
    virtual const FiniteShapePtrVariant2D resolve(ShapeID id) const = 0;
    virtual ~IShapeResolver2D() = default;
};

template <typename T>
class SlotVector2D
{
private:
    struct Slot
    {
        T value;
        //uint32_t generation = 0;
    };

public:
    struct Handle
    {
        uint32_t index;
        //uint32_t generation;
    };

    Handle add(const T& value)
    {
        uint32_t idx;

        if (!_free.empty())
        {
            // Use an empty slot
            idx = _free.back();
            _free.pop_back();

            _slots[idx].value = value;
        }
        else
        {
            // Create a new slot
            idx = static_cast<uint32_t>(_slots.size());
            _slots.push_back({ value });
        }

        return { idx };
    }

    void remove(Handle h)
    {
        if (!isValid(h)) return;

        auto& slot = _slots[h.index];
        //slot.alive = false;
        //slot.generation++;       // invalidate old handles
        _free.push_back(h.index);
    }

    void remove(uint32_t idx)
    {
        if (!isValid(idx)) return;

        auto& slot = _slots[idx];
        //slot.alive = false;
        //slot.generation++;       // invalidate old handles
        _free.push_back(idx);
    }

    bool isValid(Handle h) const
    {
        return h.index < _slots.size();
    }

    bool isValid(uint32_t idx) const
    {
        return idx < _slots.size();
    }

    T* get(Handle h)
    {
        return isValid(h) ? &_slots[h.index].value : nullptr;
    }

    const T* get(Handle h) const
    {
        return isValid(h) ? &_slots[h.index].value : nullptr;
    }

private:
    std::vector<Slot> _slots;
    std::vector<uint32_t> _free;
};

class ShapeStore2D : public IShapeResolver2D
{
public:
    ShapeID add(const Triangle2D& t)
    {
        auto h = _triangles.add(t);
        return insertHandle({ ShapeType2D::SHAPE2D_TRIANGLE, h.index } );
    }

    ShapeID add(const Rectangle2D& r)
    {
        auto h = _rectangles.add(r);
        return insertHandle({ ShapeType2D::SHAPE2D_RECTANGLE, h.index } );
    }

    ShapeID add(const Polygon2D& p)
    {
        auto h = _generalPoly.add(p);
        return insertHandle({ ShapeType2D::SHAPE2D_POLYGON, h.index } );
    }

    ShapeID add(const ConvexPolygon2D& p)
    {
        auto h = _convexPoly.add(p);
        return insertHandle({ ShapeType2D::SHAPE2D_CONVEX_POLYGON, h.index } );
    }

    ShapeID add(const Circle2D& c)
    {
        auto h = _circles.add(c);
        return insertHandle({ ShapeType2D::SHAPE2D_CIRCLE, h.index } );
    }

    ShapeID add(const FiniteShapeVariant2D& shape)
    {
        return std::visit([&](const auto &shape) -> ShapeID { return add(shape); }, shape);
    }

private:
    ShapeID insertHandle(ShapeHandle2D h)
    {
        uint32_t idx;
        uint32_t generation;
        if (!_freeIndices.empty())
        {
            // Use an empty slot
            idx = _freeIndices.back();
            _freeIndices.pop_back();
            generation = _handleMap[idx].generation;
            _handleMap[idx].handle = h;
            _handleMap[idx].isAlive = true;
        }
        else
        {
            // Create a new slot
            idx = static_cast<uint32_t>(_handleMap.size());
            generation = 0;
            _handleMap.push_back({ h, 0, true });
        }
        return { idx, generation };
    }

private:
    void _remove(ShapeHandle2D h)
    {
        switch (h.type)
        {
        case ShapeType2D::SHAPE2D_CIRCLE:    _circles.remove(h.index); break;
        case ShapeType2D::SHAPE2D_TRIANGLE:  _triangles.remove(h.index); break;
        case ShapeType2D::SHAPE2D_RECTANGLE: _rectangles.remove(h.index); break;
        case ShapeType2D::SHAPE2D_POLYGON:   _generalPoly.remove(h.index); break;
        case ShapeType2D::SHAPE2D_CONVEX_POLYGON: _convexPoly.remove(h.index); break;
        }
    }

public:
    void remove(ShapeID id)
    {
        if (isValid(id))
        {
            ShapeHandle2D h = _handleMap[id.index].handle;
            _handleMap[id.index].generation++;
            _handleMap[id.index].isAlive = false;
            _remove(h);
            _freeIndices.push_back(id.index);
        }
    }

    void removeUnsafe(uint32_t index)
    {
        if (isValid(index))
        {
            ShapeHandle2D h = _handleMap[index].handle;
            _handleMap[index].generation++;
            _handleMap[index].isAlive = false;
            _remove(h);
            _freeIndices.push_back(index);
        }
    }

    bool isValid(uint32_t index) const
    {
        return index < _handleMap.size() && _handleMap[index].isAlive;
    }

    bool isValid(ShapeID id) const
    {
        return isValid(id.index) && _handleMap[id.index].generation == id.generation;
    }

    template <typename Fn>
    void visit(ShapeID id, Fn&& fn)
    {
        if (!isValid(id)) return;
        auto h = _handleMap[id.index].handle;

        switch (h.type)
        {
        case ShapeType2D::SHAPE2D_CIRCLE:
            if (auto* p = _circles.get({h.index})) fn(*p);
            break;
        case ShapeType2D::SHAPE2D_TRIANGLE:
            if (auto* p = _triangles.get({h.index})) fn(*p);
            break;
        case ShapeType2D::SHAPE2D_RECTANGLE:
            if (auto* p = _rectangles.get({h.index})) fn(*p);
            break;
        case ShapeType2D::SHAPE2D_POLYGON:
            if (auto* p = _generalPoly.get({h.index})) fn(*p);
            break;
        case ShapeType2D::SHAPE2D_CONVEX_POLYGON:
            if (auto* p = _convexPoly.get({h.index})) fn(*p);
            break;
        }
    }

    template <typename Container, typename Fn>
    void forEachId(const Container& ids, Fn&& fn)
    {
        for (const ShapeID& id : ids)
        {
            visit(id, fn);
        }
    }

    template <typename T>
    const T* get_as(ShapeID id) const
    {
        if (!isValid(id)) return nullptr;
        auto h = _handleMap[id.index].handle;

        if (h.type != T::shapeType)
            return nullptr;

        switch(h.type)
        {
        case ShapeType2D::SHAPE2D_CIRCLE:    return _circles.get({h.index});
        case ShapeType2D::SHAPE2D_TRIANGLE:  return _triangles.get({h.index});
        case ShapeType2D::SHAPE2D_RECTANGLE: return _rectangles.get({h.index});
        case ShapeType2D::SHAPE2D_POLYGON:   return _generalPoly.get({h.index});
        case ShapeType2D::SHAPE2D_CONVEX_POLYGON: return _convexPoly.get({h.index});
        }
        return nullptr;
    }

    const FiniteShapePtrVariant2D get(ShapeID id) const
    {
        if (!isValid(id)) return std::monostate{};
        auto h = _handleMap[id.index].handle;

        switch(h.type)
        {
        case ShapeType2D::SHAPE2D_CIRCLE:    return _circles.get({h.index});
        case ShapeType2D::SHAPE2D_TRIANGLE:  return _triangles.get({h.index});
        case ShapeType2D::SHAPE2D_RECTANGLE: return _rectangles.get({h.index});
        case ShapeType2D::SHAPE2D_POLYGON:   return _generalPoly.get({h.index});
        case ShapeType2D::SHAPE2D_CONVEX_POLYGON: return _convexPoly.get({h.index});
        }
        return std::monostate{};
    }

    const FiniteShapePtrVariant2D resolve(ShapeID id) const override { return get(id); }

private:
    SlotVector2D<Circle2D>        _circles;
    SlotVector2D<Triangle2D>      _triangles;
    SlotVector2D<Rectangle2D>     _rectangles;
    SlotVector2D<Polygon2D>       _generalPoly;
    SlotVector2D<ConvexPolygon2D> _convexPoly;

private:
    std::vector<ShapeObject> _handleMap;
    std::vector<uint32_t> _freeIndices;
};


} // namespace Math

} // namespace Arns