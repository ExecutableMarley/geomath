/*
 * Copyright (c) Marley Arns
 * Licensed under the MIT License.
*/

#pragma once

#include <math.h>
#include <vector>
#include <stdexcept>

#include "CommonMath.hpp"
#include "IMatrix.hpp"

namespace Arns
{

namespace Math
{

class Matrix : public IMatrix
{
protected:
    std::vector<std::vector<real_t>> m_data;
    size_t m_rows;
    size_t m_columns;
public:

    Matrix() : m_data(), m_rows(0), m_columns(0) {}

    Matrix(size_t rows, size_t columns, real_t initValue = 0.f) : m_data(rows, std::vector<real_t>(columns, initValue)) 
    {
        this->m_rows = rows;
        this->m_columns = columns;
    }

    Matrix(size_t rows, size_t columns, const std::vector<real_t> &data) : m_data(rows, std::vector<real_t>(columns))
    {
        if (data.size() != rows * columns)
        {
            throw std::invalid_argument("Matrix data size does not match the number of rows and columns.");
        }

        this->m_rows = rows;
        this->m_columns = columns;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                m_data[i][j] = data[i * columns + j];
            }
        }
    }

    Matrix(const Matrix &other) : m_data(other.m_data), m_rows(other.m_rows), m_columns(other.m_columns) {}

    Matrix(const IMatrix &other)
    {
        m_rows = other.rows();
        m_columns = other.columns();
        m_data = std::vector<std::vector<real_t>>(m_rows, std::vector<real_t>(m_columns, 0));
        for (int i = 0; i < m_rows; i++)
            for (int j = 0; j < m_columns; j++)
                m_data[i][j] = other(i, j);
    }

    // IMatrix interface

    size_t rows()    const override { return m_rows;}

    size_t columns() const override { return m_columns;}

    real_t& operator()(size_t row, size_t column) override
    {
        return m_data[row][column];
    }

    const real_t& operator()(size_t row, size_t column) const override
    {
        return m_data[row][column];
    }

    //

    Matrix transpose() const
    {
        Matrix result(m_columns, m_rows);
        for (int i = 0; i < m_rows; i++)
            for (int j = 0; j < m_columns; j++)
                result(j, i) = m_data[i][j];

        return result;
    }
};

} // namespace Math

} // namespace Arns