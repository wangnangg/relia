#include "MarkingChain.h"
#include "MarkingChainSplit.h"
#include "helper.h"
#include "logger.h"
#include <limits>

ColSparseM subchain_to_Pmatrix(const Subchain &subchain, ChainIndexMapper *mapper = nullptr)
{
	uint_t dim = subchain.home_element_count() + subchain.foreign_element_count();
	ColSparseM Pmat(dim);
	for (uint_t i = 0; i < subchain.home_element_count(); i++)
	{
		auto ele_ptr = subchain[i];
		uint_t mat_ind = mapper == nullptr ? i : mapper->from_chain(i);
		Pmat.alloc_major(mat_ind, ele_ptr->get_out_degree());
		for (auto arc : ele_ptr->get_to_arc_list())
		{
			if (arc.val != 0.0)
			{
				auto dest_ptr = static_cast<ChainElement *>(arc.dest_ele);
				uint_t dest_mat_ind;
				if (subchain.owns(dest_ptr))
				{
					dest_mat_ind = mapper == nullptr ?
						dest_ptr->get_subindex() : mapper->from_chain(dest_ptr->get_subindex());
				}
				else
				{
					dest_mat_ind = subchain.find_foreign_element_index(dest_ptr);
				}
				Pmat.add_entry(mat_ind, dest_mat_ind, arc.val / ele_ptr->get_out_sum());
			}
		}
		Pmat.assemble(mat_ind);
	}
	for (uint_t i = subchain.home_element_count(); i < dim; i++)
	{
		Pmat.alloc_major(i, 1);
		Pmat.add_entry(i, i, 1.0);
	}
	return Pmat;
}


ColSparseM subchain_to_Qmatrix(const Subchain &subchain, ChainIndexMapper *mapper = nullptr)
{
	uint_t dim = subchain.home_element_count();
	ColSparseM Qmat(dim);
	for (uint_t i = 0; i < dim; i++)
	{
		auto ele_ptr = subchain[i];
		uint_t mat_ind = mapper == nullptr ? i : mapper->from_chain(i);
		Qmat.alloc_major(mat_ind, ele_ptr->get_out_degree() + 1);
		for (auto arc : ele_ptr->get_to_arc_list())
		{
			if (arc.val != 0.0)
			{
				auto dest_ptr = static_cast<ChainElement *>(arc.dest_ele);
				if (subchain.owns(dest_ptr))
				{
					uint_t dest_mat_ind = mapper == nullptr ?
						dest_ptr->get_subindex() : mapper->from_chain(dest_ptr->get_subindex());
					Qmat.add_entry(mat_ind, dest_mat_ind, arc.val);
				}
			}
		}
		Qmat.add_entry(mat_ind, mat_ind, -ele_ptr->get_out_sum());
		Qmat.assemble(mat_ind);
	}
	return Qmat;
}
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
	auto Pmat = subchain_to_Pmatrix(subchain);
	LOG2("solve tangible subchain:\n" << display(Pmat));
	Vector prob_sol(Pmat.dim());
	double prob_sum = norm1(subchain_init);
	for (uint_t i = 0; i < subchain_init.dim(); i++)
	{
		prob_sol(i) = subchain_init(i);
	}
	IterStopCondition stop_cond = stop_condition;
	power_method(Pmat, prob_sol, stop_cond);
	double sol_sum = 0;
	LOG2("prob = " << display(prob_sol));
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
	
	stop_cond = stop_condition;
	sor_method(solution.stay_time, Qmat_row, -1.0, std::move(subchain_init), stop_cond, 1.0);
	LOG2("stay_time = " << display(solution.stay_time));
	return solution;
}

Vector ss_solve_absorbing_subchain(const Subchain &subchain,
	double prob_sum,
	IterStopCondition &stop_condition)
{
	auto Qmat = subchain_to_Qmatrix(subchain);
	LOG2("solve absorbing subchain:\n" << display(Qmat));
	auto Qmat_row = to_row_sparse(Qmat);
	Vector sol(Qmat.dim(), 1.0);
	sor_method(sol, Qmat_row, stop_condition, 1.0);
	double norm = norm1(sol);
	double factor = prob_sum / norm;
	for (uint_t i = 0; i < sol.dim(); i++)
	{
		sol(i) = std::abs(sol(i)) * factor;
	}
	LOG2("normed result(" << prob_sum << "):" << display(sol));
	return sol;
}

inline void feed_init_prob(std::vector<Vector> &subchain_init_prob_vec,
	const std::vector<ElementProb> &input_prob_list)
{
	for (const ElementProb& init_prob : input_prob_list)
	{
		ChainElement *ele = (ChainElement *)init_prob.ele;
		Vector &v = subchain_init_prob_vec[ele->get_subchain_index()];
		v(ele->get_subindex()) += init_prob.prob;
	}
}


SSChainSolution ss_divide_solve_marking_chain(const MarkingChain<ChainElement> &chain,
	const std::vector<ElementProb> &chain_init_vec,
	const IterStopCondition& stop_condition)
{
	SSChainSolution solution(chain.size());
	LOG2("solve marking chain:\n" << display(markingchain_to_Qmatrix(chain)));
	std::vector<uint_t> start_ind;
	start_ind.reserve(chain_init_vec.size());
	for (auto prob_pair : chain_init_vec)
	{
		auto ele_ptr = (ChainElement*)prob_pair.ele;
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
			Vector sub_sol = subchain.home_element_count() == 1 ?
				std::move(init_vec) :
				ss_solve_absorbing_subchain(subchain, norm1(init_vec), stop_condition_inst);
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
	LOG2("solution prob = " << display(solution.prob));
	LOG2("solution stay time = " << display(solution.stay_time));
	return solution;
}

Vector markingchain_init_to_vector(std::vector<ElementProb> init, uint_t dim)
{
	Vector v(dim);
	for (auto e : init)
	{
		auto ele = (BasicChainElement *)e.ele;
		v(ele->get_index()) = e.prob;
	}
	return v;
}

Subchain::Subchain(uint_t index, const MarkingChain<ChainElement> &chain, const std::vector<uint_t> &ele_ind_list) :
	myindex(index),
	element_list(0),
	foreign_element_list(0)
{
	element_list.reserve(ele_ind_list.size());
	for (auto ind : ele_ind_list)
	{
		ChainElement *ele = chain[ind];
		element_list.push_back(ele);
		ele->set_subchain_index(index);
		ele->set_subindex(element_list.size() - 1);
	}
	for (auto ele_ptr : element_list)
	{
		for (auto arc : ele_ptr->get_to_arc_list())
		{
			auto dest_ele = (ChainElement *)arc.dest_ele;
			if (dest_ele->get_subchain_index() != index)
			{//foreign element
				if (foreign_map.find(dest_ele) == foreign_map.end())
				{
					// new foreign element
					foreign_map[dest_ele] = foreign_element_list.size();
					foreign_element_list.push_back(dest_ele);
				}
			}
		}
	}
}
