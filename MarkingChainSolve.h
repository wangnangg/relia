//
// Created by wangnan on 1/11/17.
//

#ifndef RELIAPY_MARKINGCHAINSOLVE_H
#define RELIAPY_MARKINGCHAINSOLVE_H

struct SSChainSolution
{
    Vector prob;
    Vector stay_time;

    SSChainSolution(uint_t dim) : prob(dim), stay_time(dim)
    {}

};

SSChainSolution ss_divide_solve_marking_chain(const MarkingChain<SubchainElement> &chain,
                                              const std::vector<ElementProb> &chain_init_vec,
                                              const IterStopCondition &stop_condition);


#endif //RELIAPY_MARKINGCHAINSOLVE_H
