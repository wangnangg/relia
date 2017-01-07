//
// Created by wangnan on 1/4/17.
//

#ifndef RELIA_PETRINETSOLUTION_H
#define RELIA_PETRINETSOLUTION_H

#include "PetriNet.h"
#include "MarkingChain.h"
#include "Matrix.h"
#include "logger.h"

struct Option
{
    enum SSMethod
    {
        SS_Auto = 0,
        SS_SOR = 1,
        SS_Power = 2
    } steady_state_method;
    enum TSMethod
    {
        TS_Auto = 0
    } transient_state_method;
    double sor_omega;
    uint_t max_interation;
    double precision;
    uint_t van_chain_max_iter;
    double van_chain_precision;
};

std::pair<Vector, ChainMatrixMapper> solve_ss_power(const MarkingChain<BasicChainElement> &chain,
                                                    const MarkingChainInitState &chain_init,
                                                    IterStopCondition &stop_condition);;

std::pair<Vector, ChainMatrixMapper> solve_ss_sor(const MarkingChain<BasicChainElement> &chain,
                                                  const MarkingChainInitState &chain_init,
                                                  IterStopCondition &stop_condition,
                                                  double omega);;

std::pair<Vector, ChainMatrixMapper> solve_ss_auto(const MarkingChain<BasicChainElement> &chain,
                                                   const MarkingChainInitState &chain_init,
                                                   IterStopCondition &stop_condition);

class PetriNetSolution
{
public:
    typedef std::function<double(PetriNetContext *)> RewardFuncType;
private:
    std::vector<double> inst_reward;
    std::vector<double> cum_reward;
    std::vector<RewardFuncType> inst_reward_func;
    std::vector<RewardFuncType> cum_reward_func;
    Option _option;
    const Option &option;
public:
    PetriNet petri_net;
    void *tag;

    PetriNetSolution() : option(_option)
    {
        LOG1(__FUNCTION__);
        _option.steady_state_method = Option::SS_Auto;
        _option.transient_state_method = Option::TS_Auto;
        _option.max_interation = 10000;
        _option.precision = 1e-10;
        _option.sor_omega = 1.0;
        _option.van_chain_max_iter = 10000;
        _option.van_chain_precision = 1e-10;
    }

    ~PetriNetSolution()
    {
        LOG1(__FUNCTION__);
    }

    //option
    void set_ss_method(Option::SSMethod method)
    {
        LOG1(__FUNCTION__);
        _option.steady_state_method = method;
    }

    void set_ts_method(Option::TSMethod method)
    {
        LOG1(__FUNCTION__);
        _option.transient_state_method = method;
    }

    void set_sor_omega(double omega)
    {
        LOG1(__FUNCTION__);
        _option.sor_omega = omega;
    }

    void set_max_iter(uint_t iter)
    {
        LOG1(__FUNCTION__);
        _option.max_interation = iter;
    }

    void set_precision(double prec)
    {
        LOG1(__FUNCTION__);
        _option.precision = prec;
    }

    uint_t add_inst_reward_func(RewardFuncType reward_func)
    {
        LOG1(__FUNCTION__);
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
        LOG1(__FUNCTION__);
        return inst_reward[reward_index];
    }

    double get_cum_reward(uint_t reward_index)
    {
        LOG1(__FUNCTION__);
        return cum_reward[reward_index];
    }

    IterStopCondition solve_steady_state()
    {
        LOG1(__FUNCTION__);
        petri_net.finalize();
        IterStopCondition van_chain_stop_condition(option.van_chain_max_iter, option.van_chain_precision);
        auto result = generate_marking_chain<BasicChainElement>(petri_net, van_chain_stop_condition);
        auto &chain = result.first;
        auto &chain_init = result.second;
        IterStopCondition stop_condition(option.max_interation, option.precision);
        Vector sol;
        ChainMatrixMapper mapper;
        switch (option.steady_state_method)
        {
            case Option::SS_Power:
                std::tie(sol, mapper) = solve_ss_power(chain, chain_init, stop_condition);
                break;
            case Option::SS_SOR:
                std::tie(sol, mapper) = solve_ss_sor(chain, chain_init, stop_condition, option.sor_omega);
                break;
            default:
                std::tie(sol, mapper) = solve_ss_auto(chain, chain_init, stop_condition);
                break;
        }
        inst_reward.clear();
        for (uint_t i = 0; i < inst_reward_func.size(); i++)
        {
            inst_reward.push_back(eval_reward(chain, mapper, sol, inst_reward_func[i]));
        }
        return stop_condition;


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
