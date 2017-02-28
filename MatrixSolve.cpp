#include "Matrix.h"
#include <utility>
#include "MatrixSolve.h"
#include "helper.h"
#include "easylogging++.h"

void sor_method(Vector &x, const RowSparseM &Q, IterStopCondition &stop_condition, double omega)
{
    TIMED_SCOPE(funcTimer, "SOR Ax=0");
	uint_t dim = Q.dim();
	RowSparseM L(dim), U(dim), D(dim);
	split_lud(L, U, D, Q);
	RowSparseM L_w(dim), U_w(dim);
	mat_add(L_w, 1.0, D, omega, L);
	mat_add(U_w, -omega, U, 1 - omega, D);
	/*solve L_w . x_next = U_w x */
	Vector x_next(dim);
	uint_t iter_count = 0;
	Vector res(dim);
	do
	{
		for (uint_t i = 0; i < stop_condition.get_check_interval(); i++)
		{
			backsolve(x_next, L_w, U_w, x);
			std::swap(x, x_next);
			iter_count += 1;
		}
		vec_add(res, 1.0, x, -1.0, x_next);
	} while (!stop_condition.should_stop(iter_count, res));
	if( !stop_condition.is_precision_reached())
	{
		LOG(WARNING) << "Precisoin not reached." << " Target: " << stop_condition.get_target_precision() << " Reached: " << stop_condition.get_reached_precision();
	}
}
void sor_method(Vector &x, const RowSparseM &Q, double alpha, Vector b, IterStopCondition &stop_condition, double omega)
{
    TIMED_SCOPE(funcTimer, "SOR Ax=b");
	uint_t dim = Q.dim();
	RowSparseM L(dim), U(dim), D(dim);
	split_lud(L, U, D, Q);
	RowSparseM L_w(dim), U_w(dim);
	mat_add(L_w, 1.0, D, omega, L);
	mat_add(U_w, -omega, U, 1 - omega, D);
	b.scale(alpha * omega);
	/*solve L_w . x_next = U_w x + b */
	Vector x_next(dim);
	uint_t iter_count = 0;
	Vector res(dim);
	do
	{
		for (uint_t i = 0; i < stop_condition.get_check_interval(); i++)
		{
			backsolve(x_next, L_w, U_w, x, b);
			std::swap(x, x_next);
			iter_count += 1;
		}
		vec_add(res, 1.0, x, -1.0, x_next);
	} while (!stop_condition.should_stop(iter_count, res));
	if( !stop_condition.is_precision_reached())
	{
		LOG(WARNING) << "Precisoin not reached." << " Target: " << stop_condition.get_target_precision() << " Reached: " << stop_condition.get_reached_precision();
	}
}

void transient_prob_unif(Vector& x, Vector v0, const ColSparseM& P, double unif_rate, double t, double precision, bool& overflowed)
{
	if(unif_rate <= 0)
	{
		x = v0;
		overflowed = false;
		return;
	}
	uint_t ss_check_interval = 50;
	Vector v1(v0.dim());

	int_t left_p, right_p;
	bool result;
	fox_find_trunc_point(unif_rate * t, left_p, right_p, precision, result);
	overflowed = result;
	auto weight = fox_term_weight(unif_rate * t, left_p, right_p, precision, result);
	double weight_sum = fox_weight_sum(weight);
	overflowed = result | overflowed;

	Vector checkpoint(v0);
	Vector diff(v0.dim());
	uint_t iter_counter = 0;
	bool ss_reached = false;
	x.fill(0);
	for (int_t k = 0; k < left_p; k++)
	{
		matvec(v1, 1.0, P, v0);
		std::swap(v1, v0);
		iter_counter += 1;
		if (iter_counter % ss_check_interval == 0)
		{
			vec_add(diff, 1.0, v0, -1.0, checkpoint);
			if (norm2(diff) < precision)
			{
				ss_reached = true;
				break;
			}
			checkpoint = v0;
		}
	}
	double summed_term = 0;
	if (!ss_reached)
	{
		for (int_t k = left_p; k <= right_p; k++)
		{
			double p_term = weight[k - left_p] / weight_sum;
			summed_term += p_term;
			vec_inc(x, p_term, v0);

			matvec(v1, 1.0, P, v0);
			std::swap(v1, v0);
			iter_counter += 1;
			if (iter_counter % ss_check_interval == 0)
			{
				vec_add(diff, 1.0, v0, -1.0, checkpoint);
				if (norm2(diff) < precision)
				{
					ss_reached = true;
					break;
				}
				checkpoint = v0;
			}
		}
	}
	if (ss_reached)
	{
		vec_inc(x, 1 - summed_term, v0);
	}
}

void transient_cum_unif(Vector& x, Vector v0, const ColSparseM& P, double unif_rate, double t, double precision)
{
	if(unif_rate <= 0)
	{
		vec_scalar(x, t, v0);
		return;
	}
	uint_t ss_check_interval = 50;
	double qt = unif_rate * t;
	uint_t right_p = find_cum_rtp(qt, precision, t);

	double log_qt = log(qt);
	double tmpti = -qt;
	double pterm = exp(tmpti);
	double right_cum = 1.0 - pterm;
	double sum_right_cum = right_cum;
	x.fill(0);
	Vector v1(v0.dim());

	vec_inc(x, right_cum, v0);
	bool ss_reached = false;
	uint_t iter_counter = 0;
	Vector diff(v0.dim());
	Vector checkpoint(v0);
	for (uint_t i = 1; i <= right_p; i++)
	{
		matvec(v1, 1.0, P, v0);
		std::swap(v1, v0);

		tmpti += log_qt - log((double)(i));
		pterm = exp(tmpti);
		right_cum -= pterm;
		sum_right_cum += right_cum;

		vec_inc(x, right_cum, v0);
		iter_counter += 1;
		if (iter_counter % ss_check_interval == 0)
		{
			vec_add(diff, 1.0, v0, -1.0, checkpoint);
			if (norm2(diff) < (precision / t))
			{
				ss_reached = true;
				break;
			}
			checkpoint = v0;
		}
	}
	x.scale(1.0 / unif_rate);
	if (ss_reached)
	{
		vec_inc(x, t - sum_right_cum / unif_rate, v0);
	}
}
