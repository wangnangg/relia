//
// Created by wangnan on 1/4/17.
//

#include "PetriNetSolution.h"
LOG_INIT;

std::vector<MarkingVal> solve_ss_auto(const PetriNet &petri_net, IterStopCondition &stop_condition)
{
    auto generate_result = generate_marking_chain<ChainElement>(petri_net, stop_condition);
    auto& chain  = generate_result.first;
    auto& chain_init = generate_result.second;
    Vector sol = solve_marking_chain(chain, chain_init, stop_condition);
    std::vector<MarkingVal> result;
    result.reserve(chain.size());
    for(uint_t i=0; i<chain.size(); i++)
    {
       result.emplace_back(&chain[i]->get_marking(), sol(i), 0.0);
    }
    return result;
}

std::vector<MarkingVal> solve_ss_power(const PetriNet &petri_net, IterStopCondition &stop_condition)
{
//    auto generate_result = generate_marking_chain<BasicChainElement>(petri_net, stop_condition);
//    auto &chain = generate_result.first;
//    auto &chain_init = generate_result.second;
//    auto Pmat = markingchain_to_Pmatrix(chain);
    return std::vector<MarkingVal>();
}

std::vector<MarkingVal> solve_ss_sor(const PetriNet &petri_net, IterStopCondition &stop_condition, double omega)
{
    return std::vector<MarkingVal>();
}
