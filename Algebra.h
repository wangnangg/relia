#pragma once

#include "Matrix.h"

class IterStopCondition
{
    uint_t max_iteration;
    double precision;
    double reached_precision = 0.0;
    uint_t used_iteration_count = 0;
    uint_t check_interval = 10;
public:
    IterStopCondition(uint_t max_iter, double precision) : max_iteration(max_iter), precision(precision)
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
    { return reached_precision < precision; }

    uint_t get_used_iter() const
    { return used_iteration_count; }

    uint_t get_check_interval() const
    { return check_interval; }
};

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
            matvec(x_next, P, x);
            std::swap(x, x_next);
            iter_count += 1;
        }
        sub(res, x, x_next);
    } while (!stop_condition.should_stop(iter_count, res));
}


/*
 * solve Qx=0
 */
void sor_method(Vector &x, const RowSparseM &Q, IterStopCondition &stop_condition, double omega);


