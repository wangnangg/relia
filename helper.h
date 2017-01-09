#pragma once

#include "Matrix.h"
#include <iostream>
#include <iomanip>


template<typename M>
M create_matrix(uint_t dim, std::vector<double> val)
{}

template<>
RowSparseM create_matrix<RowSparseM>(uint_t dim, std::vector<double> val);

template<>
ColSparseM create_matrix<ColSparseM>(uint_t dim, std::vector<double> val);

inline double get_matrix_entry(const ColSparseM &m, uint_t row_ind, uint_t col_ind)
{
    const auto &col = m.get_raw_data()[col_ind];
    for (auto pair : col)
    {
        if (pair.index == row_ind)
        {
            return pair.val;
        }
    }
    return 0.0;
}

inline double get_matrix_entry(const RowSparseM &m, uint_t row_ind, uint_t col_ind)
{
    const auto &row = m.get_raw_data()[row_ind];
    for (auto pair : row)
    {
        if (pair.index == col_ind)
        {
            return pair.val;
        } else if(pair.index > col_ind)
        {
            return 0.0;
        }
    }
    return 0.0;
}

template<typename M>
std::string display(const M &m)
{
    std::stringstream ss;
    ss << ">>>>>>>>" << std::endl;
    ss << std::setprecision(2);
    for (uint_t i = 0; i < m.dim(); i++)
    {
        for (uint_t j = 0; j < m.dim(); j++)
        {
            ss << std::setw(4) << get_matrix_entry(m, i, j) << " ";
        }
        ss << std::endl;
    }
    ss << "<<<<<<<<" << std::endl;
    return ss.str();
}


template<typename M>
std::string display_raw(const M &m)
{
    std::stringstream ss;
    ss << "raw view>>>>>>>>" << std::endl;
    ss << std::setprecision(2);
    for (uint_t i = 0; i < m.dim(); i++)
    {
        for (auto e : m(i))
        {
            ss << "(" << e.index << ", " << e.val << ") ";
        }
        ss << std::endl;
    }
    ss << "raw view<<<<<<<<" << std::endl;
    return ss.str();
}

std::string display(const Vector &v);


