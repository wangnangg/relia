#pragma once
#include "PetriNet.h"
#include <vector>
#include "Type.h"
#include <map>
#include <memory>
#include "Matrix.h"
#include "Algebra.h"
struct ChainArc
{
	void* dest_ele;
	double val;
};

class BasicChainElement
{
	Marking marking;
	std::vector<ChainArc> to_arc_list = {};
	double out_sum = 0;
	uint_t index;
public:
	BasicChainElement(const BasicChainElement&) = delete;
	BasicChainElement(Marking&& marking, uint_t index) :marking(std::move(marking)), index(index){}
	void add_arc(ChainArc arc)
	{
		to_arc_list.push_back(arc);
		out_sum += arc.val;
	}
	double get_out_sum() const
	{
		return out_sum;
	}
	uint_t get_out_degree() const
	{
		return to_arc_list.size();
	}
	const Marking& get_marking() const
	{
		return marking;
	}
	const std::vector<ChainArc>& get_to_arc_list() const
	{
		return to_arc_list;
	}
	uint_t get_index() const { return index; }
};


template <typename ElementType>
class MarkingChain
{
	std::vector<std::unique_ptr<ElementType>> element_list;
public:
	MarkingChain() = default;
	MarkingChain(MarkingChain&&) = default;
	MarkingChain(const MarkingChain&) = delete;
	MarkingChain& operator=(MarkingChain&&) = default;

	ElementType* operator[](uint_t index) const
	{
		return element_list[index].get();
	}
	uint_t size() const
	{
		return element_list.size();
	}
	ElementType* search(const Marking& marking) const
	{
		for (const auto& ele_ptr : element_list)
		{
			if (ele_ptr->get_marking() == marking)
			{
				return ele_ptr.get();
			}
		}
		return nullptr;	
	}
	ElementType* add_marking(Marking&& marking)
	{
		element_list.emplace_back(std::make_unique<ElementType>(std::move(marking), element_list.size()));
		return element_list[element_list.size() - 1].get();
	}

	bool contains(ElementType* ele) const
	{
		uint_t index = ele->get_index();
		if(index >= element_list.size())
		{
			return false;
		}
		return element_list[ele->get_index()].get() == ele;
	}
};


struct ElementInitProb
{
	void* ele;
	double prob;
};

typedef std::vector<ElementInitProb> MarkingChainInitState;


template <typename T>
class PointerIndexMapper
{
	std::vector<T*> mat2chain;
	std::map<T*, uint_t> chain2mat;
public:
	PointerIndexMapper() = default;
	void set_map(T* ele_id, uint_t mat_ind)
	{
		chain2mat[ele_id] = mat_ind;
		if (mat_ind >= mat2chain.size())
		{
			mat2chain.resize(mat_ind + 1, nullptr);
		}
		mat2chain[mat_ind] = ele_id;
	}
	uint_t pointer_to_index(T* ele_ptr) const
	{
		return chain2mat.at(ele_ptr);
	}
	T* index_to_pointer(uint_t mat_ind) const
	{
		return mat2chain[mat_ind];
	}
	bool has_pointer(T* ele_ptr) const
	{
		return chain2mat.find(ele_ptr) != chain2mat.end();
	}
	bool has_index(uint_t mat_ind) const
	{
		if (mat_ind >= mat2chain.size())
		{
			return false;
		}
		return mat2chain[mat_ind] != nullptr;
	}
	uint_t size() const
	{
		return mat2chain.size();
	}
};


class ChainMatrixMapper
{
	std::vector<uint_t> chain_to_mat;
	std::vector<uint_t> mat_to_chain;
public:
	ChainMatrixMapper(uint_t size):chain_to_mat(size), mat_to_chain(size) {}
	void set_map(uint_t chain_ind, uint_t matrix_ind)
	{
		chain_to_mat[chain_ind] = matrix_ind;
		mat_to_chain[matrix_ind] = chain_ind;
	}
	uint_t chain2mat(uint_t chain_ind) const
	{
		return chain_to_mat[chain_ind];
	}
	uint_t mat2chain(uint_t mat_ind) const
	{
		return mat_to_chain[mat_ind];
	}
	uint_t size() const
	{
		return chain_to_mat.size();
	}
};

template <typename ElementType>
ElementType* search_element(const std::vector<const MarkingChain<ElementType>*>& chain_list, const Marking& mk)
{
	ElementType* result = nullptr;
	for (const MarkingChain<ElementType>* chain : chain_list)
	{
		result = chain->search(mk);
		if (result != nullptr)
		{
			break;
		}
	}
	return result;
}

template <typename ElementType>
void explore_vanishing_marking(const PetriNet& petri_net,
	const std::vector<const MarkingChain<ElementType>*>& chain_list,
	MarkingChain<ElementType>& tan_container_chain,
	MarkingChain<ElementType>& van_container_chain,
	ElementType* explore_ele)
{
	auto next_mk_list = petri_net.next_markings(explore_ele->get_marking());
	for (auto& mk_pair : next_mk_list)
	{
		ChainArc arc;
		Marking& mk = mk_pair.first;
		arc.val = mk_pair.second;
		arc.dest_ele = search_element(chain_list, mk);
		if (arc.dest_ele == nullptr)
		{
			if (mk.type == Marking::Vanishing)
			{
				arc.dest_ele = van_container_chain.add_marking(std::move(mk));
			}
			else
			{
				arc.dest_ele = tan_container_chain.add_marking(std::move(mk));
			}
		}
		explore_ele->add_arc(arc);
	}
}

template <typename ElementType>
MarkingChain<ElementType> recursive_explore_vanishing_marking(const PetriNet& petri_net,
	std::vector<const MarkingChain<ElementType>*> chain_list,
	MarkingChain<ElementType>& tan_container_chain,
	Marking&& van_marking)
{
	MarkingChain<ElementType> van_container_chain;
	chain_list.push_back(&van_container_chain);
	van_container_chain.add_marking(std::move(van_marking));
	for (uint_t current_index = 0; current_index < van_container_chain.size(); current_index++)
	{
		explore_vanishing_marking(petri_net, chain_list, tan_container_chain, van_container_chain, van_container_chain[current_index]);
	}
	return van_container_chain;
}



extern IterStopCondition van_chain_solve_stop_condition;

template <typename ElementType>
std::vector<ChainArc> collapse_vanishing_chain(const MarkingChain<ElementType>& van_container_chain)
{
	PointerIndexMapper<ElementType> mapper = index_vanishing_chain(van_container_chain);
	uint_t dim = mapper.size();
	uint_t van_count = van_container_chain.size();
	ColSparseM P(dim);
	for (uint_t i = 0; i < van_count; i++)
	{
		auto van_ele = mapper.index_to_pointer(i);
		P.alloc_major(i, van_ele->get_out_degree());
		for (ChainArc arc : van_ele->get_to_arc_list())
		{
			uint_t matrix_ind = mapper.pointer_to_index(static_cast<ElementType*>(arc.dest_ele));
			double prob = arc.val / van_ele->get_out_sum();
			P.add_entry(i, matrix_ind, prob);
		}
		P.assemble(i);
	}
	for (uint_t i = van_count; i < dim; i++)
	{
		P.alloc_major(i, 1);
		P.add_entry(i, i, 1.0);
	}
	Vector x_sol(dim);
	x_sol(mapper.pointer_to_index(van_container_chain[0])) = 1.0;
	power_method(P, x_sol, van_chain_solve_stop_condition);

	std::vector<ChainArc> out_arc_list;
	out_arc_list.reserve(dim - van_count);
	for (uint_t i = van_count; i < dim; i++)
	{
		auto dest_chain_ele = mapper.index_to_pointer(i);
		out_arc_list.push_back(ChainArc{ dest_chain_ele, x_sol(i) });
	}
	return out_arc_list;
}


template <typename ElementType>
PointerIndexMapper<ElementType> index_vanishing_chain(
	const MarkingChain<ElementType>& van_chain
)
{
	PointerIndexMapper<ElementType> mapper;
	uint_t van_size = van_chain.size();
	uint_t tan_counter = 0;
	for (uint_t i = 0; i < van_size; i++)
	{
		ElementType* ele = van_chain[i];
		mapper.set_map(ele, i);
		for (ChainArc arc : ele->get_to_arc_list())
		{
			if (static_cast<ElementType*>(arc.dest_ele)->get_marking().type != Marking::Vanishing)
			{
				if (!mapper.has_pointer(static_cast<ElementType*>(arc.dest_ele)))
				{
					mapper.set_map(static_cast<ElementType*>(arc.dest_ele), van_size + tan_counter);
					tan_counter += 1;
				}
			}
		}
	}

	return mapper;

}




template <typename ElementType>
void explore_tangible_marking(const PetriNet& petri_net,
	const std::vector<const MarkingChain<ElementType>*>& chain_list,
	MarkingChain<ElementType>& tan_container_chain,
	ElementType* tan_ele)
{
	auto next_mk_list = petri_net.next_markings(tan_ele->get_marking());
	for (auto& mk_pair : next_mk_list)
	{
		Marking& mk = mk_pair.first;
		double val = mk_pair.second;
		if (mk.type == Marking::Vanishing)
		{
			auto van_chain = recursive_explore_vanishing_marking(petri_net, chain_list, tan_container_chain, std::move(mk));
			auto arc_list = collapse_vanishing_chain(van_chain);
			for (auto arc : arc_list)
			{
				arc.val = arc.val * val;
				tan_ele->add_arc(arc);
			}
		}
		else {
			ElementType* ele = search_element(chain_list, mk);
			if (ele == nullptr)
			{
				ele = tan_container_chain.add_marking(std::move(mk));
			}
			ChainArc arc{ ele, val };
			tan_ele->add_arc(arc);
		}
	}

}

template <typename ElementType>
std::pair<MarkingChain<ElementType>, MarkingChainInitState> generate_marking_chain(const PetriNet& petri_net)
{
	MarkingChain<ElementType> tan_container_chain;
	MarkingChainInitState init_state;
	std::vector<const MarkingChain<ElementType>*> chain_list{ &tan_container_chain };
	if (petri_net.get_init_marking().type == Marking::Vanishing)
	{
		auto van_chain = recursive_explore_vanishing_marking(petri_net, chain_list, tan_container_chain, petri_net.get_init_marking().clone());
		auto arc_list = collapse_vanishing_chain(van_chain);
		for (auto arc : arc_list)
		{
			init_state.push_back({ arc.dest_ele, arc.val });
		}
	}
	else
	{
		tan_container_chain.add_marking(petri_net.get_init_marking().clone());
		init_state.push_back({ tan_container_chain[0], 1.0 });
	}
	for (uint_t current_index = 0; current_index < tan_container_chain.size(); current_index++)
	{
		explore_tangible_marking(petri_net, chain_list, tan_container_chain, tan_container_chain[current_index]);
	}
	return std::make_pair<MarkingChain<ElementType>, MarkingChainInitState>(std::move(tan_container_chain), std::move(init_state));
}


ChainMatrixMapper index_tangible_chain(const MarkingChain<BasicChainElement>& chain);


ColSparseM markingchain_to_Qmatrix(const MarkingChain<BasicChainElement>& chain);
