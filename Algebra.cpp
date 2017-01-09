#include "Matrix.h"
#include <utility>
#include "Algebra.h"
#include "logger.h"
#include "helper.h"

void sor_method(Vector &x, const RowSparseM &Q, IterStopCondition &stop_condition, double omega)
{
    LOG2("sor solving Qx=0:\n" << display(Q))
    LOG2("start with x = :" << display(x))
    uint_t dim = Q.dim();
    RowSparseM L(dim), U(dim), D(dim);
    split_lud(L, U, D, Q);
    RowSparseM L_w(dim), U_w(dim);
    add(L_w, 1.0, D, omega, L);
    add(U_w, -omega, U, 1 - omega, D);
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
        sub(res, x, x_next);
    } while (!stop_condition.should_stop(iter_count, res));
    LOG2("sor x = :" << display(x))
}
void sor_method(Vector &x, const RowSparseM &Q, double alpha, Vector b, IterStopCondition &stop_condition, double omega)
{
    LOG2("sor solving Qx=b\n" << display(Q) << "b = " << display(b))
    LOG2("start with x = " << display(x))
    uint_t dim = Q.dim();
    RowSparseM L(dim), U(dim), D(dim);
    split_lud(L, U, D, Q);
    RowSparseM L_w(dim), U_w(dim);
    add(L_w, 1.0, D, omega, L);
    add(U_w, -omega, U, 1 - omega, D);
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
        sub(res, x, x_next);
    } while (!stop_condition.should_stop(iter_count, res));
    LOG2("sor x = " << display(x))
}
