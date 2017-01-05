//
// Created by wangnan on 1/4/17.
//

#ifndef RELIA_PETRINETSOLUTION_H
#define RELIA_PETRINETSOLUTION_H

#include "PetriNet.h"
#include "MarkingChain.h"
#include "Matrix.h"

class PetriNetSolution
{
public:
    typedef std::function<double(PetriNetContext *)> RewardFuncType;
private:
    std::vector<double> inst_reward;
    std::vector<double> cum_reward;
    std::vector<RewardFuncType> inst_reward_func;
    std::vector<RewardFuncType> cum_reward_func;
public:
    PetriNet petri_net;
    void *tag;

    uint_t add_inst_reward_func(RewardFuncType reward_func)
    {
        inst_reward_func.push_back(reward_func);
        return inst_reward_func.size() - 1;
    }

    uint_t add_cum_reward_func(RewardFuncType reward_func)
    {
        cum_reward_func.push_back(reward_func);
        return cum_reward_func.size() - 1;
    }

    double get_inst_reward(uint_t reward_index)
    {
        return inst_reward[reward_index];
    }

    double get_cum_reward(uint_t reward_index)
    {
        return cum_reward[reward_index];
    }

    IterStopCondition solve_steady_state()
    {
        auto result = generate_marking_chain<BasicChainElement>(petri_net);
        auto &chain = result.first;
        auto &chain_init = result.second;
        auto tuple = markingchain_to_Qmatrix(chain, chain_init);
        auto &QMatrix = std::get<0>(tuple);
        auto &sol = std::get<1>(tuple);
        auto &mapper = std::get<2>(tuple);
        IterStopCondition stopCondition(1000, 1e-10);
        //power_method(PMatrix, sol, stopCondition);
        sol.fill(1.0);
        auto QMat_row = to_row_sparse(QMatrix);
        sor_method(sol, QMat_row, stopCondition, 1.0);
        sol.scale(1.0 / norm1(sol));
        inst_reward.clear();
        for (uint_t i = 0; i < inst_reward_func.size(); i++)
        {
            inst_reward.push_back(eval_reward(chain, mapper, sol, inst_reward_func[i]));
        }
        return stopCondition;
    }

private:

    double eval_reward(const MarkingChain<BasicChainElement> &chain,
                       const ChainMatrixMapper &mapper,
                       const Vector &sol,
                       RewardFuncType func)
    {
        double reward = 0.0;
        for (uint_t i = 0; i < chain.size(); i++)
        {
            const auto &ele = chain[i];
            PetriNetContext context{&petri_net, &ele->get_marking()};
            double val = func(&context);
            uint_t chain_ind = ele->get_index();
            uint_t mat_ind = mapper.chain2mat(chain_ind);
            reward += sol(mat_ind) * val;
        }
        return reward;
    }

};


#endif //RELIA_PETRINETSOLUTION_H
