/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "SlotID.hpp"

namespace Arns::geomath
{


template <typename T>
class SparseSlotVector
{
public:
    // --- Types ---

    struct Entry
    {
        ShapeID id;
        T &value;
    };

    struct ConstEntry
    {
        ShapeID id;
        const T &value;
    };

    // --- Capacity ---

    size_t size() const { return _slots.size() - _free.size(); }

    size_t capacity() const { return _slots.capacity(); }

    bool empty() const { return size() == 0; }

    void reserve(size_t count)
    {
        _slots.reserve(count);
    }

    // --- Modifiers ---

    ShapeID add(const T &value)
    {
        uint32_t idx;
        uint32_t generation;

        if (!_free.empty())
        {
            // Use an empty slot
            idx = _free.back();
            _free.pop_back();

            generation = _slots[idx].generation;
            _slots[idx].value = value;
            _slots[idx].is_active = true;
        }
        else
        {
            // Create a new slot
            idx = static_cast<uint32_t>(_slots.size());
            generation = 0;
            _slots.push_back({value, generation, true});
        }

        return {idx, generation};
    }

    void remove(ShapeID h)
    {
        if (!isValid(h))
            return;

        auto &slot = _slots[h.index];
        slot.is_active = false;
        slot.generation++; // invalidate old handles
        _free.push_back(h.index);
    }

    /*
    void remove(uint32_t idx)
    {
        if (!isValid(idx))
            return;

        auto &slot = _slots[idx];
        slot.is_active = false;
        slot.generation++; // invalidate old handles
        _free.push_back(idx);
    }*/

    void clear()
    {
        _slots.clear();
        _free.clear();
    }

    // --- Access ---

    bool isValid(ShapeID h) const
    {
        return h.index < _slots.size() && _slots[h.index].is_active && h.generation == _slots[h.index].generation;
    }

    bool contains(ShapeID h) const { return isValid(h); }

    T *get(ShapeID h)
    {
        return isValid(h) ? &_slots[h.index].value : nullptr;
    }

    const T *get(ShapeID h) const
    {
        return isValid(h) ? &_slots[h.index].value : nullptr;
    }

    T &operator[](ShapeID h)
    {
        assert(isValid(h));
        return _slots[h.index].value;
    }
    const T &operator[](ShapeID h) const
    {
        assert(isValid(h));
        return _slots[h.index].value;
    }

    void append_active_ids(std::vector<ShapeID> &out_ids) const
    {
        out_ids.reserve(out_ids.size() + this->size());

        for (uint32_t i = 0; i < _slots.size(); ++i)
        {
            if (_slots[i].is_active)
            {
                out_ids.push_back(ShapeID{i, _slots[i].generation});
            }
        }
    }

public:
    template <bool IsConst>
    class IteratorImpl
    {
    private:
        using ContainerPtr = std::conditional_t<IsConst, const SparseSlotVector *, SparseSlotVector *>;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::conditional_t<IsConst, ConstEntry, Entry>;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type *;
        using reference = value_type; // Proxy value reference

        IteratorImpl() : _index(0), _container(nullptr) {}

        IteratorImpl(uint32_t index, ContainerPtr container)
            : _index(index), _container(container)
        {
            skip_inactive();
        }

        reference operator*() const
        {
            const auto &slot = _container->_slots[_index];
            if constexpr (IsConst)
            {
                return ConstEntry{ShapeID{_index, slot.generation}, slot.value};
            }
            else
            {
                // Note: slot.value is mutable if !IsConst
                return Entry{ShapeID{_index, slot.generation}, const_cast<T &>(slot.value)};
            }
        }

        IteratorImpl &operator++()
        {
            ++_index;
            skip_inactive();
            return *this;
        }

        IteratorImpl operator++(int)
        {
            auto copy = *this;
            ++(*this);
            return copy;
        }

        friend bool operator==(const IteratorImpl &a, const IteratorImpl &b)
        {
            return a._index == b._index && a._container == b._container;
        }

        friend bool operator!=(const IteratorImpl &a, const IteratorImpl &b)
        {
            return !(a == b);
        }

    private:
        void skip_inactive()
        {
            if (!_container)
                return;

            const auto size = static_cast<uint32_t>(_container->_slots.size());
            while (_index < size && !_container->_slots[_index].is_active)
            {
                ++_index;
            }
        }

        uint32_t _index = 0;
        ContainerPtr _container = nullptr;
    };

    // --- Iteration ---

    using iterator = IteratorImpl<false>;
    using const_iterator = IteratorImpl<true>;

    iterator begin()
    {
        return {0, this};
    }

    iterator end()
    {
        return {static_cast<uint32_t>(_slots.size()), this};
    }

    const_iterator begin() const
    {
        return {0, this};
    }

    const_iterator end() const
    {
        return {static_cast<uint32_t>(_slots.size()), this};
    }

private:
    struct Slot
    {
        T value;
        uint32_t generation = 0;
        bool is_active = false;
    };

    std::vector<Slot> _slots;
    std::vector<uint32_t> _free;
};


}