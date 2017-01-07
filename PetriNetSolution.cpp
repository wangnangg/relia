//
// Created by wangnan on 1/4/17.
//

#include "PetriNetSolution.h"
LOG_INIT;
std::pair<Vector, ChainMatrixMapper>
solve_ss_power(const MarkingChain<BasicChainElement> &chain, const MarkingChainInitState &chain_init,
               IterStopCondition &stop_condition)
{
    auto tuple = markingchain_to_Pmatrix(chain, chain_init);
    auto &PMatrix = std::get<0>(tuple);
    auto &sol = std::get<1>(tuple);
    auto &mapper = std::get<2>(tuple);
    power_method(PMatrix, sol, stop_condition);
    return std::make_pair(std::move(sol), std::move(mapper));
}

std::pair<Vector, ChainMatrixMapper>
solve_ss_sor(const MarkingChain<BasicChainElement> &chain, const MarkingChainInitState &chain_init,
             IterStopCondition &stop_condition, double omega)
{
    auto tuple = markingchain_to_Qmatrix(chain, chain_init);
    auto &QMatrix = std::get<0>(tuple);
    auto &sol = std::get<1>(tuple);
    auto &mapper = std::get<2>(tuple);
    sol.fill(1.0);
    auto QMat_row = to_row_sparse(QMatrix);
    sor_method(sol, QMat_row, stop_condition, omega);
    sol.scale(1.0 / norm1(sol));
    return std::make_pair(std::move(sol), std::move(mapper));
}

std::pair<Vector, ChainMatrixMapper>
solve_ss_auto(const MarkingChain<BasicChainElement> &chain, const MarkingChainInitState &chain_init,
              IterStopCondition &stop_condition)
{
    return solve_ss_sor(chain, chain_init, stop_condition, 1.0);
}
