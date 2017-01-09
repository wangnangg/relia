//
// Created by wangnan on 1/4/17.
//

#include "PetriNetSolution.h"
LOG_INIT;

template <typename ElementType>
std::vector<MarkingVal> translate_solution(MarkingChain<ElementType> chain, const SSChainSolution& sol)
{
	std::vector<MarkingVal> result;
    result.reserve(chain.size());
    for(uint_t i=0; i<chain.size(); i++)
    {
       result.emplace_back(std::move(chain[i]->get_marking()), sol.prob(i), sol.stay_time(i));
    }
    return result;
}
std::vector<MarkingVal> solve_ss_divide(const PetriNet &petri_net, const IterStopCondition &stop_condition)
{
	auto generate_result = generate_marking_chain<ChainElement>(petri_net, stop_condition);
    auto& chain  = generate_result.first;
    auto& chain_init = generate_result.second;
    SSChainSolution sol = ss_divide_solve_marking_chain(chain, chain_init, stop_condition);
	return translate_solution(std::move(chain), sol);
}

void PetriNetSolution::solve_steady_state()
{
	LOG1(__FUNCTION__);
	petri_net.finalize();
	IterStopCondition stop_condition(option.max_interation, option.precision);
	std::vector<MarkingVal> result;
	switch (option.steady_state_method)
	{
	case Option::SS_Power:
		result = solve_ss_power(petri_net, stop_condition);
	case Option::SS_SOR:
		result = solve_ss_sor(petri_net, stop_condition, option.sor_omega);
	case Option::SS_Divide:
		result = solve_ss_divide(petri_net, stop_condition);
	default:
		result = solve_ss_auto(petri_net, stop_condition);
		break;
	}
	update_reward(result);
}

void PetriNetSolution::update_reward(const std::vector<MarkingVal>& result)
{
	clear_reward();
	for (const MarkingVal& m_val : result)
	{
		PetriNetContext context{&petri_net, &m_val.marking};
		for (uint_t i = 0; i < inst_reward_func.size(); i++)
		{
			double reward_rate = inst_reward_func[i](&context);
			inst_reward[i] += reward_rate * m_val.prob;
		}
		for (uint_t i = 0; i < cum_reward_func.size(); i++)
		{
			double reward_rate = cum_reward_func[i](&context);
			if(reward_rate != 0.0)
			{
				cum_reward[i] += reward_rate * m_val.stay_time;
			}
		}
	}
}

std::vector<MarkingVal> solve_ss_power(const PetriNet &petri_net, const IterStopCondition &stop_condition)
{
//    auto generate_result = generate_marking_chain<BasicChainElement>(petri_net, stop_condition);
//    auto &chain = generate_result.first;
//    auto &chain_init = generate_result.second;
//    auto Pmat = markingchain_to_Pmatrix(chain);
    return std::vector<MarkingVal>();
}

std::vector<MarkingVal> solve_ss_sor(const PetriNet &petri_net, const IterStopCondition &stop_condition, double omega)
{
    return std::vector<MarkingVal>();
}


std::vector<MarkingVal> solve_ss_auto(const PetriNet &petri_net, const IterStopCondition &stop_condition)
{
	return solve_ss_divide(petri_net, stop_condition);
}


