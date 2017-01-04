#pragma once
#include "Matrix.h"
#include <iostream>
#include <iomanip>


template <typename M>
M create_matrix(uint_t dim, std::vector<double> val) {}

template <>
RowSparseM create_matrix<RowSparseM>(uint_t dim, std::vector<double> val);

template <>
ColSparseM create_matrix<ColSparseM>(uint_t dim, std::vector<double> val);

inline double get_matrix_entry(const ColSparseM& m, uint_t row_ind, uint_t col_ind)
{
	const auto& col = m.get_raw_data()[col_ind];
	for (auto pair : col)
	{
		if (pair.index == row_ind)
		{
			return pair.val;
		}
	}
	return 0.0;
}

inline double get_matrix_entry(const RowSparseM& m, uint_t row_ind, uint_t col_ind)
{
	const auto& row = m.get_raw_data()[row_ind];
	for (auto pair : row)
	{
		if (pair.index == col_ind)
		{
			return pair.val;
		}
	}
	return 0.0;
}

template <typename M>
void display(const M& m)
{
	std::cout << ">>>>>>>>" << std::endl;
	std::cout << std::setprecision(2);
	for (uint_t i = 0; i < m.dim(); i++)
	{
		for (uint_t j = 0; j < m.dim(); j++)
		{
			std::cout << std::setw(4) << get_matrix_entry(m, i, j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "<<<<<<<<" << std::endl;
}


template <typename M>
void display_raw(const M& m)
{
	std::cout << "raw view>>>>>>>>" << std::endl;
	std::cout << std::setprecision(2);
	for (uint_t i = 0; i < m.dim(); i++)
	{
		for(auto e : m(i))
		{
			std::cout << "(" << e.index << ", " << e.val << ") ";
		}
		std::cout << std::endl;
	}
	std::cout << "raw view<<<<<<<<" << std::endl;
}

void display(const Vector& v);


