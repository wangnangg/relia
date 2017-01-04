#include <gtest/gtest.h>
#include "Algebra.h"
#include "helper.h"

TEST(Algebra, power_method)
{
	ColSparseM m0 = create_matrix<ColSparseM>(3, {
			0, 2.0 / 3.0, 0,
			1.0, 0, 0,
			0, 1.0 / 3.0, 1
	});
	Vector v0 = Vector({1, 0, 0});
	IterStopCondition stop_condtion(100, 1e-6);
	power_method(m0, v0, stop_condtion);
	display(v0);
	std::cout << "used iteration:" << stop_condtion.get_used_iter() << std::endl;
	ASSERT_TRUE(stop_condtion.is_precision_reached());

	ColSparseM m1 = create_matrix<ColSparseM>(4, {
			0, 0, 0, 0,
			0.5, 0, 0, 0,
			0.5, 0.5, 0, 0,
			0, 0.5, 1, 1
	});
	Vector v1 = Vector({1, 0, 0, 0});
	IterStopCondition stop_condtion1(10, 1e-10);
	power_method(m1, v1, stop_condtion1);
	display(v1);
	std::cout << "used iteration:" << stop_condtion1.get_used_iter() << std::endl;
	ASSERT_TRUE(stop_condtion1.is_precision_reached());

}

TEST(Algebra, sor_method)
{
	auto Q = create_matrix<RowSparseM>(4, {
			-1, 2, 0, 0,
			1, -3, 2, 0,
			0, 1, -3, 2,
			0, 0, 1, -2
	});
	IterStopCondition stop_condtion(100, 1e-6);
	Vector x(4);
	x.fill(1.0);
	sor_method(x, Q, stop_condtion, 1.0);
	x.scale(1.0 / norm1(x));
	display(x);
	std::cout << "used iteration:" << stop_condtion.get_used_iter() << std::endl;
	ASSERT_TRUE(stop_condtion.is_precision_reached());
}