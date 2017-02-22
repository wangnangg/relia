#include "Matrix.h"
#include <tuple>
#include <algorithm>
#include "helper.h"
#include <cmath>

void SparseMatrix::assemble(uint_t i)
{
	auto &major = data[i];
	std::sort(major.begin(), major.end(), [](SparseEntry a, SparseEntry b)
	{
		return a.index < b.index;
	});
	uint_t working_slot = 0;
	uint_t next_slot = 1;
	while (next_slot < major.size())
	{
		if (major[next_slot].index != major[working_slot].index)
		{
			working_slot += 1;
			major[working_slot] = major[next_slot];
		}
		else
		{
			major[working_slot].val += major[next_slot].val;
		}
		next_slot += 1;
	}
	major.resize(working_slot + 1);
}

RowSparseM to_row_sparse(const ColSparseM &col_matrix)
{
	uint_t dim = col_matrix.dim();
	std::vector<uint_t> counter(dim, 0);
	for (uint_t i = 0; i < dim; i++)
	{
		for (auto e : col_matrix(i))
		{
			counter[e.index] += 1;
		}
	}
	RowSparseM result(dim);
	for (uint_t i = 0; i < dim; i++)
	{
		result.alloc_major(i, counter[i]);
	}
	for (uint_t i = 0; i < dim; i++)
	{
		for (auto e : col_matrix(i))
		{
			result.add_entry(e.index, i, e.val);
		}
	}
	for (uint_t i = 0; i < dim; i++)
	{
		result.assemble(i);
	}
	return result;
}

ColSparseM to_col_sparse(const RowSparseM &row_matrix)
{
	uint_t dim = row_matrix.dim();
	std::vector<uint_t> counter(dim, 0);
	for (uint_t i = 0; i < dim; i++)
	{
		for (auto e : row_matrix(i))
		{
			counter[e.index] += 1;
		}
	}
	ColSparseM result(dim);
	for (uint_t i = 0; i < dim; i++)
	{
		result.alloc_major(i, counter[i]);
	}
	for (uint_t i = 0; i < dim; i++)
	{
		for (auto e : row_matrix(i))
		{
			result.add_entry(e.index, i, e.val);
		}
	}
	for (uint_t i = 0; i < dim; i++)
	{
		result.assemble(i);
	}
	return result;
}

void matvec(Vector& v0, double alpha, const ColSparseM& A, const Vector& v1)
{
	v0.fill(0.0);
	for (uint_t col_i = 0; col_i < v0.dim(); col_i++)
	{
		for (const auto pair : A(col_i))
		{
			v0(pair.index) += alpha * pair.val * v1(col_i);
		}
	}
}

void matvec(Vector& v0, double alpha, const RowSparseM& A, const Vector& v1)
{
	v0.fill(0.0);
	for (uint_t row_i = 0; row_i < v0.dim(); row_i++)
	{
		for (const auto pair : A(row_i))
		{
			v0(row_i) += alpha * pair.val * v1(pair.index);
		}
	}
}

void vec_add(Vector& v0, double alpha, const Vector& v1, double beta, const Vector& v2)
{
	for (uint_t i = 0; i < v0.dim(); i++)
	{
		v0(i) = alpha * v1(i) + beta * v2(i);
	}
}

void vec_inc(Vector& v0, double alpha, const Vector& v1)
{
	for (uint_t i = 0; i < v0.dim(); i++)
	{
		v0(i) += alpha * v1(i);
	}
}

int_t next_index(const std::vector<SparseEntry> &row_a, uint_t &a_j,
	const std::vector<SparseEntry> &row_b, uint_t &b_j)
{
	if (a_j >= row_a.size() || b_j >= row_b.size())
	{
		return -1;
	}
	int_t index;
	auto a_e = row_a[a_j];
	auto b_e = row_b[b_j];
	if (a_e.index < b_e.index)
	{
		index = a_e.index;
		a_j += 1;
	}
	else if (a_e.index > b_e.index)
	{
		index = b_e.index;
		b_j += 1;
	}
	else
	{
		index = a_e.index;
		a_j += 1;
		b_j += 1;
	}
	return index;
}

SparseEntry next_entry(const std::vector<SparseEntry> &row_a, uint_t &a_j,
	const std::vector<SparseEntry> &row_b, uint_t &b_j, double alpha, double beta)
{
	SparseEntry entry(0, 0);
	if (a_j >= row_a.size())
	{
		auto b_e = row_b[b_j];
		entry.index = b_e.index;
		entry.val = b_e.val * beta;
		b_j += 1;
		return entry;
	}
	if (b_j >= row_b.size())
	{
		auto a_e = row_a[a_j];
		entry.index = a_e.index;
		entry.val = a_e.val * alpha;
		a_j += 1;
		return entry;
	}
	auto a_e = row_a[a_j];
	auto b_e = row_b[b_j];
	if (a_e.index < b_e.index)
	{
		entry.index = a_e.index;
		entry.val = a_e.val * alpha;
		a_j += 1;
	}
	else if (a_e.index > b_e.index)
	{
		entry.index = b_e.index;
		entry.val = b_e.val * beta;
		b_j += 1;
	}
	else
	{
		entry.index = a_e.index;
		entry.val = a_e.val * alpha + b_e.val * beta;
		a_j += 1;
		b_j += 1;
	}
	return entry;
}


double norm1(const Vector &x)
{
	double v = 0;
	for (auto e : x.get_raw_data())
	{
		v += std::abs(e);
	}
	return v;
}

double norm2(const Vector &x)
{
	double v = 0;
	for (auto e : x.get_raw_data())
	{
		v += e * e;
	}
	return std::sqrt(v);
}

void vec_scalar(Vector& v0, double alpha, Vector& v1)
{
	for (uint_t i = 0; i < v0.dim(); i++)
	{
		v0(i) = alpha * v1(i);
	}
}

bool vec_eq(const Vector& v0, const Vector& v1, double error)
{
	if (v0.dim() != v1.dim())
	{
		return false;
	}
	for (uint_t i = 0; i < v0.dim(); i++)
	{
		if (std::abs(v0(i) - v1(i)) > error)
		{
			return false;
		}
	}
	return true;
}

void split_lud(RowSparseM &L, RowSparseM &U, RowSparseM &D, const RowSparseM &A)
{
	uint_t dim = A.dim();
	for (uint_t i = 0; i < dim; i++)
	{
		L.clear(i);
		U.clear(i);
		D.clear(i);
		uint_t l_counter = 0;
		uint_t u_counter = 0;
		uint_t diag_counter = 0;
		for (auto e : A(i))
		{
			if (e.index < i)
			{
				l_counter += 1;
			}
			else if (e.index > i)
			{
				u_counter += 1;
			}
			else
			{
				diag_counter += 1;
			}
		}
		L.alloc_major(i, l_counter);
		U.alloc_major(i, u_counter);
		D.alloc_major(i, diag_counter);
	}
	for (uint_t i = 0; i < dim; i++)
	{
		for (auto e : A(i))
		{
			if (e.index < i)
			{
				L.add_entry(i, e.index, e.val);
			}
			else if (e.index > i)
			{
				U.add_entry(i, e.index, e.val);
			}
			else
			{
				D.add_entry(i, e.index, e.val);
			}
		}
	}
}

void backsolve(Vector &x, const RowSparseM &L, double alpha, const Vector &b)
{
	uint_t dim = x.dim();
	for (uint_t i = 0; i < dim; i++)
	{
		const auto &l_row_i = L(i);
		double right_entry = alpha * b(i);
		uint_t j;
		for (j = 0; j < l_row_i.size() - 1; j++)
		{
			auto e = l_row_i[j];
			right_entry -= e.val * x(e.index);
		}
		x(i) = right_entry / l_row_i[j].val;
	}
}

void backsolve(Vector &x, const RowSparseM &L, const RowSparseM &B, const Vector &b)
{
	uint_t dim = x.dim();
	for (uint_t i = 0; i < dim; i++)
	{
		const auto &b_row_i = B(i);
		const auto &l_row_i = L(i);
		double right_entry = 0;
		for (auto e : b_row_i)
		{
			right_entry += e.val * b(e.index);
		}
		uint_t j;
		for (j = 0; j < l_row_i.size() - 1; j++)
		{
			auto e = l_row_i[j];
			right_entry -= e.val * x(e.index);
		}
		x(i) = right_entry / l_row_i[j].val;
	}
}

void backsolve(Vector& x, const RowSparseM& L, const RowSparseM& B, const Vector& b, const Vector& c)
{
	uint_t dim = x.dim();
	for (uint_t i = 0; i < dim; i++)
	{
		const auto &b_row_i = B(i);
		const auto &l_row_i = L(i);
		double right_entry = 0;
		for (auto e : b_row_i)
		{
			right_entry += e.val * b(e.index);
		}
		right_entry += c(i);
		uint_t j;
		for (j = 0; j < l_row_i.size() - 1; j++)
		{
			auto e = l_row_i[j];
			right_entry -= e.val * x(e.index);
		}
		x(i) = right_entry / l_row_i[j].val;
	}

}

