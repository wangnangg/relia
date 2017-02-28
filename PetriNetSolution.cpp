//
// Created by wangnan on 1/4/17.
//

#include "PetriNetSolution.h"
#include "MarkingChainSolve.h"
#include "easylogging++.h"

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
	auto generate_result = generate_marking_chain<SubchainElement>(petri_net, stop_condition);
    auto& chain  = generate_result.first;
    auto& chain_init = generate_result.second;
    SSChainSolution sol = ss_divide_solve_marking_chain(chain, chain_init, stop_condition);
	return translate_solution(std::move(chain), sol);
}

void PetriNetSolution::solve_steady_state()
{
	petri_net.finalize();
	IterStopCondition stop_condition(option.max_interation, option.precision, option.check_interval);
	std::vector<MarkingVal> result;
	switch (option.steady_state_method)
	{
	case Option::SS_Power:
		result = solve_ss_power(petri_net, stop_condition);
		break;
	case Option::SS_SOR:
		result = solve_ss_sor(petri_net, stop_condition, option.sor_omega);
		break;
	case Option::SS_Divide:
		result = solve_ss_divide(petri_net, stop_condition);
		break;
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

void PetriNetSolution::solve_transient_state(double time)
{
    petri_net.finalize();
    IterStopCondition stop_condition(option.max_interation, option.precision, option.check_interval);
    std::vector<MarkingVal> result;
    auto chain_pair = generate_marking_chain<SubchainElement>(petri_net, stop_condition);
    auto& chain = chain_pair.first;
    auto& chain_init = chain_pair.second;
    auto Q = markingchain_to_Qmatrix(chain);
    double unif_rate = max_out_rate_sum(chain);
    unif_Qmatrix(Q, unif_rate);
    SSChainSolution sol(Q.dim());
    auto init_vec = chain_init_to_vector(chain_init, chain.size());
    bool overflowed;
    transient_prob_unif(sol.prob, Vector(init_vec), Q, unif_rate, time, option.precision, overflowed);
    if(overflowed)
    {
        LOG(WARNING) << "unification overflowed.";
    }
    transient_cum_unif(sol.stay_time, Vector(init_vec), Q, unif_rate, time, option.precision);
    auto solution = translate_solution(std::move(chain), sol);
    update_reward(solution);
}

std::vector<MarkingVal> solve_ss_power(const PetriNet &petri_net, const IterStopCondition &stop_condition)
{
	auto generate_result = generate_marking_chain<BasicChainElement>(petri_net, stop_condition);
    auto &chain = generate_result.first;
    auto &chain_init = generate_result.second;
	auto Pmat = markingchain_to_Qmatrix(chain);
	unif_Qmatrix(Pmat, max_out_rate_sum(chain) * 1.1);

//    auto generate_result = generate_marking_chain<BasicChainElement>(petri_net, stop_condition);
//    auto &chain = generate_result.first;
//    auto &chain_init = generate_result.second;
//    auto Pmat = markingchain_to_Pmatrix(chain);
    return std::vector<MarkingVal>();
}

std::vector<MarkingVal> solve_ss_sor(const PetriNet &petri_net, const IterStopCondition &stop_condition, double omega)
{
	auto generate_result = generate_marking_chain<BasicChainElement>(petri_net, stop_condition);
    auto &chain = generate_result.first;
    auto &chain_init = generate_result.second;
	auto Qmat = markingchain_to_Qmatrix(chain);
	auto Qmat_row = to_row_sparse(Qmat);
	SSChainSolution solution(Qmat_row.dim());
	solution.prob.fill(1.0);
	IterStopCondition cond1(stop_condition);
	sor_method(solution.prob, Qmat_row, cond1, 1.0);
	solution.prob.scale(1.0 / norm1(solution.prob));
	auto b = chain_init_to_vector(chain_init, Qmat_row.dim());
	IterStopCondition cond2(stop_condition);
	solution.stay_time.fill(1.0);
	sor_method(solution.stay_time, Qmat_row, -1.0, std::move(b), cond2, 1.0);

	return translate_solution(std::move(chain), solution);
}


std::vector<MarkingVal> solve_ss_auto(const PetriNet &petri_net, const IterStopCondition &stop_condition)
{
	return solve_ss_divide(petri_net, stop_condition);
}


