//
// Created by wangnan on 1/11/17.
//

#include "MarkingChain.h"
#include "MarkingChainSplit.h"
#include "MarkingChainSolve.h"

#include <limits>

struct SSTangibleSubchainSolution
{
	Vector stay_time;
	std::vector<ElementProb> prob;
	SSTangibleSubchainSolution(uint_t dim): stay_time(dim), prob() {}
};

SSTangibleSubchainSolution ss_solve_tangible_subchain(const Subchain &subchain,
	Vector subchain_init,
	const IterStopCondition& stop_condition)
{
	SSTangibleSubchainSolution solution(subchain.home_element_count());
    if(subchain.home_element_count() == 1)
    {
        SubchainElement* home_ele = subchain[0];
        double init_prob = subchain_init(0);
        solution.prob.reserve(subchain.foreign_element_count());
        for(auto arc : home_ele->get_to_arc_list())
        {
            SubchainElement* ele_ptr = (SubchainElement*)arc.dest_ele;
            solution.prob.emplace_back(ele_ptr, arc.val / home_ele->get_out_sum() * init_prob);
        }
        solution.stay_time(0) = init_prob / home_ele->get_out_sum();
        return solution;
    }
	auto Pmat = subchain_to_Pmatrix(subchain);
	Vector prob_sol(Pmat.dim());
	double prob_sum = norm1(subchain_init);
	for (uint_t i = 0; i < subchain_init.dim(); i++)
	{
		prob_sol(i) = subchain_init(i);
	}
	IterStopCondition stop_cond = stop_condition;
	power_method(Pmat, prob_sol, stop_cond);
	double sol_sum = 0;
	for (uint_t i = subchain.home_element_count(); i < Pmat.dim(); i++)
	{
		prob_sol(i) = std::abs(prob_sol(i));
		sol_sum += prob_sol(i);
	}
	double factor = prob_sum / sol_sum;
	solution.prob.reserve(subchain.foreign_element_count());
	for (uint_t i = subchain.home_element_count(); i < Pmat.dim(); i++)
	{
		solution.prob.emplace_back(subchain[i], prob_sol(i) * factor);
	}

	auto Qmat = subchain_to_Qmatrix(subchain);
	auto Qmat_row = to_row_sparse(Qmat);
	solution.stay_time.fill(0.0);

	IterStopCondition sor_stop_cond(stop_condition);
	sor_method(solution.stay_time, Qmat_row, -1.0, std::move(subchain_init), sor_stop_cond, 1.0);
	return solution;
}

Vector ss_solve_absorbing_subchain(const Subchain &subchain,
	double prob_sum,
	IterStopCondition &stop_condition)
{
    if(subchain.home_element_count() == 1)
    {
        return Vector(1, prob_sum);
    }
	auto Qmat = subchain_to_Qmatrix(subchain);
	auto Qmat_row = to_row_sparse(Qmat);
	Vector sol(Qmat.dim(), 1.0);
	sor_method(sol, Qmat_row, stop_condition, 1.0);
	double norm = norm1(sol);
	double factor = prob_sum / norm;
	for (uint_t i = 0; i < sol.dim(); i++)
	{
		sol(i) = std::abs(sol(i)) * factor;
	}
	return sol;
}

inline void feed_init_prob(std::vector<Vector> &subchain_init_prob_vec,
	const std::vector<ElementProb> &input_prob_list)
{
	for (const ElementProb& init_prob : input_prob_list)
	{
		SubchainElement *ele = (SubchainElement *)init_prob.ele;
		Vector &v = subchain_init_prob_vec[ele->get_subchain_index()];
		v(ele->get_subindex()) += init_prob.prob;
	}
}


SSChainSolution ss_divide_solve_marking_chain(const MarkingChain<SubchainElement> &chain,
	const std::vector<ElementProb> &chain_init_vec,
	const IterStopCondition& stop_condition)
{
	TIMED_SCOPE(funcTimer, "Divide method");
	SSChainSolution solution(chain.size());
	std::vector<uint_t> start_ind;
	start_ind.reserve(chain_init_vec.size());
	for (auto prob_pair : chain_init_vec)
	{
		auto ele_ptr = (SubchainElement*)prob_pair.ele;
		start_ind.push_back(ele_ptr->get_index());
	}
	auto subchain_list = split_to_subchains(chain, start_ind);
	std::vector<Vector> subchain_init_prob_vec;
	subchain_init_prob_vec.reserve(subchain_list.size());
	for (uint_t i = 0; i < subchain_list.size(); i++)
	{
		subchain_init_prob_vec.emplace_back(subchain_list[i].home_element_count());
	}
	feed_init_prob(subchain_init_prob_vec, chain_init_vec);
	//eval topology sort order
	for (int_t i = subchain_list.size() - 1; i >= 0; i--)
	{
		Vector &init_vec = subchain_init_prob_vec[i];
		if (norm1(init_vec) == 0.0) //unreachable subchain
		{
			continue;
		}
		const Subchain &subchain = subchain_list[i];
		if (subchain.is_absorbing())
		{
			IterStopCondition stop_condition_inst = stop_condition;
			Vector sub_sol = ss_solve_absorbing_subchain(subchain, norm1(init_vec), stop_condition_inst);
			for (uint_t local_ind = 0; local_ind < subchain.home_element_count(); local_ind++)
			{
				uint_t global_index = subchain[local_ind]->get_index();
				double prob = sub_sol(local_ind);
				solution.prob(global_index) = prob;
				solution.stay_time(global_index) = std::numeric_limits<double>::infinity();
			}
		}
		else
		{
			IterStopCondition stop_condition_inst = stop_condition;
			auto sub_sol = ss_solve_tangible_subchain(subchain, std::move(init_vec), stop_condition_inst);
			for (uint_t local_ind = 0; local_ind < subchain.home_element_count(); local_ind++)
			{
				uint_t global_index = subchain[local_ind]->get_index();
				double stay_time = sub_sol.stay_time(local_ind);
				solution.stay_time(global_index) = stay_time;
				//prob = 0
			}
			feed_init_prob(subchain_init_prob_vec, sub_sol.prob);
		}
	}
	return solution;
}



