#pragma once

#include "Type.h"
#include <vector>

class Vector
{
	std::vector<double> data;
public:
	Vector(uint_t size, double val = 0.0) : data(size, val)
	{}
	Vector() : data(0)
	{}
	explicit Vector(const Vector &) = default;
	Vector(Vector &&) = default;


	Vector &operator=(Vector &&) = default;

	Vector &operator=(const Vector &) = default;

	Vector(std::vector<double> &&data) : data(std::move(data))
	{}

	void scale(double val)
	{
		for (double &e : data)
		{
			e = e * val;
		}
	}

	uint_t dim() const
	{
		return data.size();
	}

	void fill(double val)
	{
		std::fill(data.begin(), data.end(), val);
	}

	double &operator()(uint_t i)
	{
		return data[i];
	}

	double operator()(uint_t i) const
	{
		return data[i];
	}

	const std::vector<double> &get_raw_data() const
	{
		return data;
	}
};

struct SparseEntry
{
	uint_t index;
	double val;

	SparseEntry(uint_t index, double val) : index(index), val(val)
	{}

	SparseEntry() : index(0), val(0.0)
	{}
};


class SparseMatrix
{
	std::vector<std::vector<SparseEntry>> data; //non diag part

public:
	SparseMatrix(uint_t dim) : data(dim)
	{}

	SparseMatrix(SparseMatrix &&) = default;
	SparseMatrix(const SparseMatrix &) = delete;

	void clear(uint_t ind)
	{
		data[ind].clear();
	}

	void alloc_major(uint_t ind, uint_t size)
	{
		data[ind].reserve(size);
	}

	void add_entry(uint_t major_ind, uint_t minor_ind, double value)
	{
		data[major_ind].emplace_back(minor_ind, value);
	}

	void add_entry(uint_t major_ind, SparseEntry entry)
	{
		data[major_ind].push_back(entry);
	}

	const std::vector<std::vector<SparseEntry>> &get_raw_data() const
	{
		return data;
	}

	uint_t dim() const
	{
		return data.size();
	}

	const std::vector<SparseEntry> &operator()(uint_t i) const
	{
		return data[i];
	}

	void scale(double v)
	{
		for (auto& maj : data)
		{
			for (auto& e : maj)
			{
				e.val = e.val * v;
			}
		}
	}
	void assemble(uint_t i);
};

class ColSparseM : public SparseMatrix
{
public:
	ColSparseM(uint_t dim) : SparseMatrix(dim)
	{}
	ColSparseM(ColSparseM &&) = default;
	ColSparseM &operator=(ColSparseM &&m) = default;
};

class RowSparseM : public SparseMatrix
{
public:
	RowSparseM(uint_t dim) : SparseMatrix(dim)
	{}
	RowSparseM(RowSparseM &&) = default;
	RowSparseM &operator=(RowSparseM &&m) = default;
};

RowSparseM to_row_sparse(const ColSparseM &col_matrix);

ColSparseM to_col_sparse(const RowSparseM &row_matrix);


template <typename M>
M Identity(uint_t dim)
{
	M matrix(dim);
	for(uint_t i=0; i<dim; i++)
	{
		matrix.add_entry(i, i, 1.0);
	}
	return matrix;
}
/*
 * norm of vector
 */
double norm1(const Vector &x);
double norm2(const Vector &x);

/**
 * \brief compute v0 = alpha * v1
 * \param v0 
 * \param alpha 
 * \param v1 
 */
void vec_scalar(Vector& v0, double alpha, Vector& v1);

/*
 * m0 = alpha * m1;
 */
template <typename MType>
void mat_scale(MType& m0, double alpha, const MType& m1)
{
	for (uint_t i = 0; i < m0.dim(); i++)
	{
		m0.alloc_major(i, m1(i).size());
		for (auto v : m1(i))
		{
			v.val *= alpha;
			m0.add_entry(i, v);
		}
	}
}

/*
 * v0 == v1? within the error bound
 */
bool vec_eq(const Vector& v0, const Vector& v1, double error);

/*
 * m0 == m1? within the error bound
 */
template <typename MType>
bool mat_eq(const MType& m0, const MType& m1, double error)
{
	if (m0.dim() != m1.dim())
	{
		return false;
	}
	for (uint_t i = 0; i < m0.dim(); i++)
	{
		const auto& major0 = m0(i);
		const auto& major1 = m1(i);
		if (major0.size() != major1.size())
		{
			return false;
		}
		for (uint_t j = 0; j < major0.size(); j++)
		{
			const auto entry0 = major0[j];
			const auto entry1 = major1[j];
			if (entry0.index != entry1.index)
			{
				return false;
			}
			if (std::abs(entry0.val - entry1.val) > error)
			{
				return false;
			}
		}
	}
	return true;
}




/*
 * v0 = alpha * v1 + beta * v2
 */
void vec_add(Vector &v0, double alpha, const Vector &v1, double beta, const Vector &v2);

/*
 * v0 += alpha * v1;
 */
void vec_inc(Vector& v0, double alpha, const Vector& v1);



int_t next_index(const std::vector<SparseEntry> &row_a, uint_t &a_j,
	const std::vector<SparseEntry> &row_b, uint_t &b_j);
SparseEntry next_entry(const std::vector<SparseEntry> &row_a, uint_t &a_j,
	const std::vector<SparseEntry> &row_b, uint_t &b_j, double alpha, double beta);
/*
 * C = alpha*A + beta*B;
 */
template <typename MType>
void mat_add(MType &C, double alpha, const MType &A, double beta, const MType &B)
{
	uint_t dim = A.dim();
	for (uint_t i = 0; i < dim; i++)
	{
		uint_t c_counter;
		uint_t a_j;
		uint_t b_j;
		const auto &a_row = A(i);
		const auto &b_row = B(i);
		a_j = 0;
		b_j = 0;
		c_counter = 0;
		while (next_index(a_row, a_j, b_row, b_j) != -1)
		{
			c_counter += 1;
		}
		c_counter += b_row.size() - b_j + a_row.size() - a_j;
		C.clear(i);
		C.alloc_major(i, c_counter);
		a_j = 0;
		b_j = 0;
		for (uint_t j = 0; j < c_counter; j++)
		{
			C.add_entry(i, next_entry(a_row, a_j, b_row, b_j, alpha, beta));
		}
	}
}

/*
 * v0 = alpha * A * v1
 */
void matvec(Vector &v0, double alpha, const ColSparseM &A, const Vector &v1);
void matvec(Vector &v0, double alpha, const RowSparseM &A, const Vector &v1);



/*
 * split A into L, U and D
 */
void split_lud(RowSparseM &L, RowSparseM &U, RowSparseM &D, const RowSparseM &A);

/*
 * solve Lx=alpha*b where L is a lower triangular matrix.
 */
void backsolve(Vector &x, const RowSparseM &L, double alpha, const Vector &b);

/*
 * solve Lx=B*b
 */
void backsolve(Vector &x, const RowSparseM &L, const RowSparseM &B, const Vector &b);

/*
 * solve Lx=B*b + c
 */
void backsolve(Vector &x, const RowSparseM &L, const RowSparseM &B, const Vector &b, const Vector& c);






