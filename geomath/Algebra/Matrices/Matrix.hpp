/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include "CommonMath.hpp"
#include "IMatrix.hpp"

#include <math.h>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <ostream>
#include <iomanip>

namespace Arns
{

namespace Math
{

class Matrix : public IMatrix
{
protected:
    std::vector<real_t> m_data;
    size_t m_rows;
    size_t m_columns;
    
public:

    Matrix() : m_data(), m_rows(0), m_columns(0) {}

    Matrix(size_t rows, size_t columns, real_t initValue = 0.f) 
        : m_data(rows * columns, initValue), m_rows(rows), m_columns(columns) {}

    Matrix(size_t rows, size_t columns, const std::vector<real_t> &data) 
        : m_rows(rows), m_columns(columns), m_data(data)
    {
        if (data.size() != rows * columns)
        {
            throw std::invalid_argument("Matrix data size does not match the number of rows and columns.");
        }
    }

    Matrix(std::initializer_list<std::initializer_list<real_t>> list)
    {
        m_rows = list.size();
        m_columns = m_rows > 0 ? list.begin()->size() : 0;

        m_data.reserve(m_rows * m_columns);
        for (const auto &row : list)
        {
            if (row.size() != m_columns)
            {
                throw std::invalid_argument("Rows must have uniform dimensions.");
            }
            m_data.insert(m_data.end(), row.begin(), row.end());
        }
    }

    Matrix(const Matrix &other) 
        : m_data(other.m_data), m_rows(other.m_rows), m_columns(other.m_columns) {}

    Matrix(const IMatrix &other)
        : m_rows(other.rows()), m_columns(other.columns()), m_data(other.rows() * other.columns())
    {
        for (size_t i = 0; i < m_rows; ++i)
        {
            for (size_t j = 0; j < m_columns; ++j)
            {
                m_data[i * m_columns + j] = other(i, j);
            }
        }
    }

    // IMatrix interface

    size_t rows()    const override { return m_rows;}

    size_t columns() const override { return m_columns;}

    real_t& operator()(size_t row, size_t column) override
    {
        return m_data[row * m_columns + column];
    }

    const real_t& operator()(size_t row, size_t column) const override
    {
        return m_data[row * m_columns + column];
    }

    //

    Matrix transpose() const
    {
        Matrix result(m_columns, m_rows);
        for (int i = 0; i < m_rows; i++)
            for (int j = 0; j < m_columns; j++)
                result(j, i) = m_data[i * m_columns + j];

        return result;
    }

    //

    friend std::ostream &operator<<(std::ostream &os, const Matrix &mat)
    {
        if (mat.rows() == 0 || mat.columns() == 0) {
            return os << "[] (Empty Matrix)";
        }

        for (size_t i = 0; i < mat.rows(); ++i) {
            os << "[ ";
            for (size_t j = 0; j < mat.columns(); ++j) {
                os << std::setw(8) << mat(i, j);
                if (j < mat.columns() - 1) {
                    os << " ";
                }
            }
            os << " ]";

            if (i < mat.rows() - 1) {
                os << "\n";
            }
        }
        return os;
    }
};

} // namespace Math

} // namespace Arns