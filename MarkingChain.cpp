#include "MarkingChain.h"
#include "MarkingChainSplit.h"
#include "helper.h"
#include "logger.h"
#include <limits>

ColSparseM subchain_to_Pmatrix(const Subchain &subchain, ChainIndexMapper *mapper)
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
				auto dest_ptr = static_cast<SubchainElement *>(arc.dest_ele);
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


ColSparseM subchain_to_Qmatrix(const Subchain &subchain, ChainIndexMapper *mapper)
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
				auto dest_ptr = static_cast<SubchainElement *>(arc.dest_ele);
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



Subchain::Subchain(uint_t index, const MarkingChain<SubchainElement> &chain, const std::vector<uint_t> &ele_ind_list) :
	myindex(index),
	element_list(0),
	foreign_element_list(0)
{
	element_list.reserve(ele_ind_list.size());
	for (auto ind : ele_ind_list)
	{
		SubchainElement *ele = chain[ind];
		element_list.push_back(ele);
		ele->set_subchain_index(index);
		ele->set_subindex(element_list.size() - 1);
	}
	for (auto ele_ptr : element_list)
	{
		for (auto arc : ele_ptr->get_to_arc_list())
		{
			auto dest_ele = (SubchainElement *)arc.dest_ele;
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