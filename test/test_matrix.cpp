#include "gtest/gtest.h"
#include "Matrix.h"
#include "helper.h"
bool equal(const Vector& v0, const Vector& v1, double error = 0.0001)
{
	for(uint_t i=0; i <v0.dim(); i++)
	{
		if(std::abs(v0(i) - v1(i)) > error)
		{
			return false;
		}
	}
	return true;
}

TEST(Matrix, create)
{
	ColSparseM m0 = create_matrix<ColSparseM>(3, {
		0, 2.0 / 3.0, 0,
			1.0, 0, 0,
				0, 1.0 / 3.0, 1
	});
	display(m0);
	ASSERT_EQ(get_matrix_entry(m0, 0, 1), 2.0 / 3.0);
	ASSERT_EQ(get_matrix_entry(m0, 1, 0), 1.0);
	ASSERT_EQ(get_matrix_entry(m0, 2, 2), 1.0);

	Vector v0 = Vector({ 0,1,2 });
	display(v0);
	ASSERT_EQ(v0(0),  0.0);
	ASSERT_EQ(v0(1),  1.0);
	ASSERT_EQ(v0(2),  2.0);
}

TEST(Matrix, basic_op)
{
	ColSparseM m0 = create_matrix<ColSparseM>(3, {
		0, 2.0 / 3.0, 0,
			1.0, 0, 0,
				0, 1.0 / 3.0, 1
	});
	Vector v0 = Vector({ 1,0,0 });
	Vector v1(3);
	matvec(v1, m0, v0);
	display(v1);
	ASSERT_TRUE(equal(v1, Vector({ 0.0, 1.0, 0.0 })));
	sub(v0, v0, v1);
	display(v0);
	ASSERT_TRUE(equal(v0, Vector({1.0, -1.0, 0.0})));

}

TEST(Matrix, backsolve)
{
	Vector x(4);
	auto m = create_matrix<RowSparseM>(4, {
		1, 0, 0, 0,
		1, 2, 0, 0,
		3, 4, 5, 0,
		6, 7, 8, 9
	});
	display(m);
	Vector b({ 4, 3, 2, 1 });
	backsolve(x, m, 1.5, b);
	display(x);
	Vector ans({ 6, -0.75, -2.4, -1.116666667 });
	ASSERT_TRUE(equal(ans, x));
}

TEST(Matrix, backsolve2)
{
	Vector x(4);
	auto m = create_matrix<RowSparseM>(4, {
		1, 0, 0, 0,
		1, 2, 0, 0,
		3, 4, 5, 0,
		6, 7, 8, 9
	});
	display(m);
	auto mr = create_matrix<RowSparseM>(4, {
		1,2,3,4,
		5,6,7,8,
		9,10,11,12,
		13,14,15,16
	});
	display(mr);
	Vector b({ 4, 3, 2, 1 });
	backsolve(x, m, 1.5, mr, b);
	display(x);
	Vector ans({ 30, 30, -12, -9.3333333});
	ASSERT_TRUE(equal(ans, x));
}

TEST(Matrix, matrix_add)
{
	RowSparseM C(4);
	auto A = create_matrix<RowSparseM>(4, {
		1,2,3,0,
		0,0,0,0,
		0,0,1,2,
		1,2,3,4
	});
	auto B = create_matrix<RowSparseM>(4, {
		0,0,0,0,
		3,0,1,0,
		5,2,1,0,
		4,3,2,1
	});
	add(C, 1.0, A, 1.0, B);
	display(C);
	add(C, 2.0, A, 3.0, B);
	display(C);
}


TEST(Matrix, split_lud)
{
	RowSparseM L(4), U(4), D(4);
	auto A = create_matrix<RowSparseM>(4, {
		1,2,3,4,
		5,6,7,8,
		9,10,11,12,
		13,14,15,16
	});
	split_lud(L, U, D, A);
	display(A);
	display(L);
	display(U);
	display(D);
}

TEST(Matrix, assemble)
{
	RowSparseM m(4);
	m.add_entry(0, 2, 2.0);
	m.add_entry(0, 1, 0.5);
	m.add_entry(0, 3, 1.0);
	m.add_entry(0, 3, 1.0);
	m.add_entry(0, 1, 0.5);
	m.add_entry(0, 3, 1.0);
	m.add_entry(1, 3, 1.0);
	m.add_entry(1, 2, 1.0);
	m.add_entry(1, 1, 1.0);
	m.add_entry(1, 0, 1.0);
	display_raw(m);
	m.assemble(0);
	m.assemble(1);
	display_raw(m);
	
}

TEST(Matrix, convert)
{
	auto A = create_matrix<RowSparseM>(4, {
		1,2,3,4,
		5,6,7,8,
		9,10,11,12,
		13,14,15,16
	});
	auto B = to_col_sparse(A);
	auto A_ = to_row_sparse(B);
	display(A);
	display(B);
	display(A_);
	display_raw(A);
	display_raw(B);
	display_raw(A_);
}