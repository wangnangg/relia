#pragma once
#include "Type.h"
#include <vector>

class Vector
{
private:
	std::vector<double> data;
	Vector(const Vector&) = default;
public:
	Vector(uint_t size) : data(size) {}
	Vector(Vector&&) = default;
	Vector& operator=(Vector&& v) = default;
	Vector(std::vector<double>&& data) :data(std::move(data)) {}
	void scale(double val)
	{
		for(double& e: data)
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
	double& operator()(uint_t i)
	{
		return data[i];
	}
	double operator()(uint_t i) const 
	{
		return data[i];
	}
	const std::vector<double>& get_raw_data() const { return data; }
};

struct SparseMatrixEntry
{
	uint_t index;
	double val;
	SparseMatrixEntry(uint_t index, double val):index(index), val(val){}
	SparseMatrixEntry():index(0), val(0.0){}
};
class SparseMatrix
{
private:
	std::vector<std::vector<SparseMatrixEntry>> data;
	SparseMatrix(const SparseMatrix&) = default;
public:	
	SparseMatrix(uint_t dim):data(dim) {}
	SparseMatrix(SparseMatrix&&) = default;
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
	void add_entry(uint_t major_ind, SparseMatrixEntry entry)
	{
		data[major_ind].push_back(entry);
	}
	const std::vector<std::vector<SparseMatrixEntry>>& get_raw_data() const
	{
		return data;
	}
	uint_t dim() const { return data.size(); }

	const std::vector<SparseMatrixEntry>& operator()(uint_t i) const { return data[i]; }

	void assemble(uint_t i);
};

class ColSparseM : public SparseMatrix
{
public:
	ColSparseM(uint_t dim) : SparseMatrix(dim) {}
	ColSparseM(ColSparseM&&) = default;
	ColSparseM& operator=(ColSparseM&& m) = default;
};

class RowSparseM : public SparseMatrix
{
public:
	RowSparseM(uint_t dim) : SparseMatrix(dim) {}
	RowSparseM(RowSparseM&&) = default;
	RowSparseM& operator=(RowSparseM&& m) = default;
};

RowSparseM to_row_sparse(const ColSparseM& col_matrix);

ColSparseM to_col_sparse(const RowSparseM& row_matrix);

/*
 * v0 = A * v1;
 */
void matvec(Vector& v0, const RowSparseM& A, const Vector& v1);
void matvec(Vector& v0, const ColSparseM& A, const Vector& v1);

/*
 * v0 = v1 - v2;
 */
void sub(Vector& v0, const Vector& v1, const Vector& v2);

/*
 * y = alpha * A * x + y;
 */
void matvec(Vector& y, double alpha, const ColSparseM& A, const Vector& x);
void matvec(Vector& y, double alpha, const RowSparseM& A, const Vector& x);

double norm1(const Vector& x);
double norm2(const Vector& x);


/*
 * C = alpha*A + beta*B;
 */
void add(RowSparseM& C, double alpha, const RowSparseM& A, double beta, const RowSparseM& B);

/*
 * split A into L, U and D
 */
void split_lud(RowSparseM& L, RowSparseM& U, RowSparseM& D, const RowSparseM& A);

/*
 * solve Lx=alpha*b where L is a lower triangular matrix.
 */
void backsolve(Vector& x, const RowSparseM& L, double alpha, const Vector& b);

/*
 * solve Lx=alpha*B*b
 */
void backsolve(Vector& x, const RowSparseM& L, double alpha, const RowSparseM& B, const Vector& b);



