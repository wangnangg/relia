#pragma once

#include "Matrix.h"
#include "helper.h"
#include "Poisson.h"
#include <cmath>
class IterStopCondition
{
	const uint_t max_iteration;
	const double precision;
	double reached_precision;
	uint_t used_iteration_count;
	uint_t check_interval;
public:
	IterStopCondition() = delete;
	IterStopCondition(const IterStopCondition& other) = default;

	IterStopCondition(uint_t max_iter, double precision, uint_t check_interval) :
		max_iteration(max_iter),
		precision(precision),
		reached_precision(NAN),
		used_iteration_count(0),
		check_interval(check_interval)
	{}

	bool should_stop(uint_t iter_count, const Vector &res)
	{
		double res_norm = norm2(res);
		if (iter_count > max_iteration)
		{
			used_iteration_count = iter_count;
			reached_precision = res_norm;
			return true;
		}
		if (res_norm < precision)
		{
			used_iteration_count = iter_count;
			reached_precision = res_norm;
			return true;
		}
		return false;
	}

	bool is_precision_reached() const
	{
		return reached_precision < precision;
	}

	uint_t get_used_iter() const
	{
		return used_iteration_count;
	}

	uint_t get_check_interval() const
	{
		return check_interval;
	}
	uint_t get_max_iter() const
	{
		return max_iteration;
	}
	double get_target_precision() const
	{
		return precision;
	}

	double get_reached_precision() const
	{
		return reached_precision;
	}

	void assert_precision_reached() const
	{
		if (is_precision_reached())
		{
			return;
		}
		std::stringstream ss;
		ss << "Precision not reached." << std::endl;
		ss << "Used iteration:" << used_iteration_count << std::endl;
		ss << "Reached precision:" << reached_precision << std::endl;
		ss << "Target precision:" << precision << std::endl;
		throw Exception(ss.str());
	}
};


/*solve Px=x*/
template<typename M>
void power_method(const M &P, Vector &x, IterStopCondition &stop_condition)
{
	Vector x_next(x.dim());
	Vector res(x.dim());
	uint_t iter_count = 0;
	do
	{
		for (uint_t i = 0; i < stop_condition.get_check_interval(); i++)
		{
			matvec(x_next, 1.0, P, x);
			std::swap(x, x_next);
			iter_count += 1;
		}
		vec_add(res, 1.0, x, -1.0, x_next);
	} while (!stop_condition.should_stop(iter_count, res));
	stop_condition.assert_precision_reached();
}

/*
 * solve Qx=0
 */
void sor_method(Vector &x, const RowSparseM &Q, IterStopCondition &stop_condition, double omega);

/*
 * solve Qx=alpha * b
 */
void sor_method(Vector &x, const RowSparseM &Q, double alpha, Vector b, IterStopCondition &stop_condition, double omega);



/**
 * \brief compute x = v0 * e^(Qt)
 * \param x
 * \param v0 
 * \param P = Q / unif_rate + I
 * \param unif_rate 
 * \param t 
 * \param precision required precision
 * \param overflowed indicating if overflow happened
 */
void transient_prob_unif(Vector& x, Vector v0, const ColSparseM& P, double unif_rate, double t, double precision, bool& overflowed);
void transient_cum_unif(Vector& x, Vector v0, const ColSparseM& P, double unif_rate, double t, double precision);
