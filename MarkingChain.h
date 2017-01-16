#pragma once

#include "PetriNet.h"
#include <vector>
#include "Type.h"
#include <map>
#include <memory>
#include "Matrix.h"
#include "Algebra.h"
#include <cmath>

struct ChainArc
{
    void *dest_ele;
    double val;
};

class BasicChainElement
{
    Marking marking;
    std::vector<ChainArc> to_arc_list = {};
    double out_sum = 0;
    uint_t index;
public:
    BasicChainElement(const BasicChainElement &) = delete;

    BasicChainElement(Marking &&marking, uint_t index) : marking(std::move(marking)), index(index)
    {}

    void add_arc(ChainArc arc)
    {
        //self arc is not allowed
        if (arc.dest_ele != this)
        {
            to_arc_list.push_back(arc);
            out_sum += arc.val;
        }
    }

    double get_out_sum() const
    {
        return out_sum;
    }

    uint_t get_out_degree() const
    {
        return to_arc_list.size();
    }

    const Marking &get_marking() const
    {
        return marking;
    }

    Marking &get_marking()
    {
        return marking;
    }

    const std::vector<ChainArc> &get_to_arc_list() const
    {
        return to_arc_list;
    }

    uint_t get_index() const
    {
        return index;
    }
};

class SubchainElement : public BasicChainElement
{
    int_t subindex;
    int_t subchain_index;
public:
    SubchainElement(Marking &&marking, uint_t index) : BasicChainElement(std::move(marking), index),
                                                    subindex(-1),
                                                    subchain_index(-1)
    {}

    void set_subindex(int_t ind)
    {
        subindex = ind;
    }

    void set_subchain_index(int_t ind)
    {
        subchain_index = ind;
    }


    uint_t get_subindex() const
    {
        return subindex;
    }

    uint_t get_subchain_index() const
    {
        return subchain_index;
    }

};

struct mk_ptr_less
{
    bool operator()(const Marking *lhs, const Marking *rhs) const
    {
        for (uint_t i = 0; i < lhs->token_list.size(); i++)
        {
            if (lhs->token_list[i] < rhs->token_list[i])
            {
                return true;
            } else if (lhs->token_list[i] > rhs->token_list[i])
            {
                return false;
            }
        }
        return false;
    }
};

template<typename ElementType>
class MarkingChain
{
    std::vector<std::unique_ptr<ElementType>> element_list;
    std::map<const Marking *, ElementType *, mk_ptr_less> mk_map;
public:
    MarkingChain() = default;

    MarkingChain(MarkingChain &&) = default;

    MarkingChain(const MarkingChain &) = delete;

    MarkingChain &operator=(MarkingChain &&) = default;

    ElementType *operator[](uint_t index) const
    {
        return element_list[index].get();
    }

    uint_t size() const
    {
        return element_list.size();
    }

    //TODO: faster search
    ElementType *search(const Marking &marking) const
    {
        auto ele = mk_map.find(&marking);
        if (ele == mk_map.end())
        {
            return nullptr;
        } else
        {
            return ele->second;
        }
    }

    ElementType *add_marking(Marking &&marking)
    {
        element_list.emplace_back(std::make_unique<ElementType>(std::move(marking), element_list.size()));
        mk_map[&(element_list.back().get()->get_marking())] = element_list.back().get();
        return element_list.back().get();
    }

    bool contains(ElementType *ele) const
    {
        uint_t index = ele->get_index();
        if (index >= element_list.size())
        {
            return false;
        }
        return element_list[ele->get_index()].get() == ele;
    }

    void clear()
    {
        element_list.clear();
        mk_map.clear();
    }

};

class Subchain
{
    uint_t myindex;
    std::vector<SubchainElement *> element_list;
    std::vector<SubchainElement *> foreign_element_list;
    std::map<SubchainElement *, uint_t> foreign_map;
public:
    Subchain(uint_t index, const MarkingChain<SubchainElement> &chain, const std::vector<uint_t> &ele_ind_list);

    bool is_absorbing() const
    {
        return foreign_element_list.size() == 0;
    }

    bool owns(SubchainElement *ele) const
    {
        return ele->get_subchain_index() == myindex;
    }

    uint_t home_element_count() const
    {
        return element_list.size();
    }

    uint_t foreign_element_count() const
    {
        return foreign_element_list.size();
    }

    SubchainElement *operator[](uint_t ind) const
    {
        if (ind < element_list.size())
        {
            return element_list[ind];
        }
        return foreign_element_list[ind - element_list.size()];
    }

    uint_t find_foreign_element_index(SubchainElement *ele) const
    {
        return foreign_map.at(ele) + element_list.size();
    }
};

struct ElementProb
{
    void *ele;
    double prob;

    ElementProb(void *ele, double prob) : ele(ele), prob(prob)
    {}
};

typedef std::vector<ElementProb> MarkingChainSparseState;


template<typename T>
class PointerIndexMapper
{
    std::vector<T *> mat2chain;
    std::map<T *, uint_t> chain2mat;
public:
    PointerIndexMapper() = default;

    void set_map(T *ele_id, uint_t mat_ind)
    {
        chain2mat[ele_id] = mat_ind;
        if (mat_ind >= mat2chain.size())
        {
            mat2chain.resize(mat_ind + 1, nullptr);
        }
        mat2chain[mat_ind] = ele_id;
    }

    uint_t pointer_to_index(T *ele_ptr) const
    {
        return chain2mat.at(ele_ptr);
    }

    T *index_to_pointer(uint_t mat_ind) const
    {
        return mat2chain[mat_ind];
    }

    bool has_pointer(T *ele_ptr) const
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


class ChainIndexMapper
{
    std::vector<uint_t> chain_to_index;
    std::vector<uint_t> index_to_chain;
public:
    ChainIndexMapper() : chain_to_index(0), index_to_chain(0)
    {}

    ChainIndexMapper(const ChainIndexMapper &) = delete;

    ChainIndexMapper(uint_t size) : chain_to_index(size), index_to_chain(size)
    {}

    ChainIndexMapper(ChainIndexMapper &&) = default;

    ChainIndexMapper &operator=(ChainIndexMapper &&) = default;

    ChainIndexMapper &operator=(const ChainIndexMapper &) = delete;

    void set_map(uint_t chain_ind, uint_t other_ind)
    {
        chain_to_index[chain_ind] = other_ind;
        index_to_chain[other_ind] = chain_ind;
    }

    uint_t from_chain(uint_t chain_ind) const
    {
        return chain_to_index[chain_ind];
    }

    uint_t to_chain(uint_t other_ind) const
    {
        return index_to_chain[other_ind];
    }

    uint_t size() const
    {
        return chain_to_index.size();
    }
};

template<typename ElementType>
ElementType *search_element(const std::vector<const MarkingChain<ElementType> *> &chain_list, const Marking &mk)
{
    ElementType *result = nullptr;
    for (const MarkingChain<ElementType> *chain : chain_list)
    {
        result = chain->search(mk);
        if (result != nullptr)
        {
            break;
        }
    }
    return result;
}

template<typename ElementType>
void explore_vanishing_marking(const PetriNet &petri_net,
                               const std::vector<const MarkingChain<ElementType> *> &chain_list,
                               MarkingChain<ElementType> &tan_container_chain,
                               MarkingChain<ElementType> &van_container_chain,
                               ElementType *explore_ele)
{
    auto next_mk_list = petri_net.next_markings(explore_ele->get_marking());
    for (auto &mk_pair : next_mk_list)
    {
        ChainArc arc;
        Marking &mk = mk_pair.first;
        arc.val = mk_pair.second;
        arc.dest_ele = search_element(chain_list, mk);
        if (arc.dest_ele == nullptr)
        {
            if (mk.type == Marking::Vanishing)
            {
                arc.dest_ele = van_container_chain.add_marking(std::move(mk));
            } else
            {
                arc.dest_ele = tan_container_chain.add_marking(std::move(mk));
            }
        }
        explore_ele->add_arc(arc);
    }
}

template<typename ElementType>
MarkingChain<ElementType> recursive_explore_vanishing_marking(const PetriNet &petri_net,
                                                              std::vector<const MarkingChain<ElementType> *> chain_list,
                                                              MarkingChain<ElementType> &tan_container_chain,
                                                              Marking &&van_marking)
{
    MarkingChain<ElementType> van_container_chain;
    chain_list.push_back(&van_container_chain);
    van_container_chain.add_marking(std::move(van_marking));
    for (uint_t current_index = 0; current_index < van_container_chain.size(); current_index++)
    {
        explore_vanishing_marking(petri_net, chain_list, tan_container_chain, van_container_chain,
                                  van_container_chain[current_index]);
    }
    return van_container_chain;
}


template<typename ElementType>
std::vector<ChainArc> collapse_vanishing_chain(const MarkingChain<ElementType> &van_container_chain,
                                               IterStopCondition &van_chain_solve_stop_condition)
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
            uint_t matrix_ind = mapper.pointer_to_index(static_cast<ElementType *>(arc.dest_ele));
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
        out_arc_list.push_back(ChainArc{dest_chain_ele, x_sol(i)});
    }
    return out_arc_list;
}


template<typename ElementType>
PointerIndexMapper<ElementType> index_vanishing_chain(
        const MarkingChain<ElementType> &van_chain
)
{
    PointerIndexMapper<ElementType> mapper;
    uint_t van_size = van_chain.size();
    uint_t tan_counter = 0;
    for (uint_t i = 0; i < van_size; i++)
    {
        ElementType *ele = van_chain[i];
        mapper.set_map(ele, i);
        for (ChainArc arc : ele->get_to_arc_list())
        {
            if (static_cast<ElementType *>(arc.dest_ele)->get_marking().type != Marking::Vanishing)
            {
                if (!mapper.has_pointer(static_cast<ElementType *>(arc.dest_ele)))
                {
                    mapper.set_map(static_cast<ElementType *>(arc.dest_ele), van_size + tan_counter);
                    tan_counter += 1;
                }
            }
        }
    }

    return mapper;

}


template<typename ElementType>
void explore_tangible_marking(const PetriNet &petri_net,
                              const std::vector<const MarkingChain<ElementType> *> &chain_list,
                              MarkingChain<ElementType> &tan_container_chain,
                              ElementType *tan_ele,
                              const IterStopCondition van_chain_stop_condition)
{
    auto next_mk_list = petri_net.next_markings(tan_ele->get_marking());
    for (auto &mk_pair : next_mk_list)
    {
        Marking &mk = mk_pair.first;
        double val = mk_pair.second;
        if (mk.type == Marking::Vanishing)
        {
            auto van_chain = recursive_explore_vanishing_marking(petri_net, chain_list, tan_container_chain,
                                                                 std::move(mk));
            IterStopCondition stop_cond = van_chain_stop_condition;
            auto arc_list = collapse_vanishing_chain(van_chain, stop_cond);
            for (auto arc : arc_list)
            {
                arc.val = arc.val * val;
                tan_ele->add_arc(arc);
            }
        } else
        {
            ElementType *ele = search_element(chain_list, mk);
            if (ele == nullptr)
            {
                ele = tan_container_chain.add_marking(std::move(mk));
            }
            ChainArc arc{ele, val};
            tan_ele->add_arc(arc);
        }
    }

}

template<typename ElementType>
std::pair<MarkingChain<ElementType>, MarkingChainSparseState> generate_marking_chain(const PetriNet &petri_net,
                                                                                     const IterStopCondition van_chain_stop_condition)
{
    MarkingChain<ElementType> tan_container_chain;
    MarkingChainSparseState init_state;
    std::vector<const MarkingChain<ElementType> *> chain_list{&tan_container_chain};
    if (petri_net.get_init_marking().type == Marking::Vanishing)
    {
        auto van_chain = recursive_explore_vanishing_marking(petri_net, chain_list, tan_container_chain,
                                                             petri_net.get_init_marking().clone());
        IterStopCondition stop_cond = van_chain_stop_condition;
        auto arc_list = collapse_vanishing_chain(van_chain, stop_cond);
        for (auto arc : arc_list)
        {
            init_state.push_back({arc.dest_ele, arc.val});
        }
    } else
    {
        tan_container_chain.add_marking(petri_net.get_init_marking().clone());
        init_state.push_back({tan_container_chain[0], 1.0});
    }
    for (uint_t current_index = 0; current_index < tan_container_chain.size(); current_index++)
    {
        explore_tangible_marking(petri_net, chain_list, tan_container_chain, tan_container_chain[current_index],
                                 van_chain_stop_condition);
    }

    LOG1("reachability graph state count:" << tan_container_chain.size());
    return std::make_pair<MarkingChain<ElementType>, MarkingChainSparseState>(std::move(tan_container_chain),
                                                                              std::move(init_state));
}

template<typename ElementType>
ColSparseM markingchain_to_Qmatrix(
        const MarkingChain<ElementType> &chain,
        ChainIndexMapper *mapper = nullptr)
{
    uint_t dim = chain.size();
    ColSparseM Qmat(dim);
    for (uint_t i = 0; i < dim; i++)
    {
        auto ele_ptr = chain[i];
        uint_t mat_ind = mapper == nullptr ? i : mapper->from_chain(i);
        if (ele_ptr->get_out_degree() != 0)
        {
            Qmat.alloc_major(mat_ind, ele_ptr->get_out_degree() + 1);
            for (auto arc : ele_ptr->get_to_arc_list())
            {
                if (arc.val != 0.0)
                {
                    auto dest_ptr = static_cast<ElementType *>(arc.dest_ele);
                    uint_t dest_mat_ind = mapper == nullptr ?
                                          dest_ptr->get_index() : mapper->from_chain(dest_ptr->get_index());
                    Qmat.add_entry(mat_ind, dest_mat_ind, arc.val);
                }
            }
            Qmat.add_entry(mat_ind, mat_ind, -ele_ptr->get_out_sum());
            Qmat.assemble(mat_ind);
        }
    }
    return Qmat;
}


template<typename ElementType>
ColSparseM markingchain_to_Pmatrix(const MarkingChain<ElementType> &chain, ChainIndexMapper *mapper = nullptr)
{
    uint_t dim = chain.size();
    ColSparseM Pmat(dim);
    for (uint_t i = 0; i < dim; i++)
    {
        auto ele_ptr = chain[i];
        uint_t mat_ind = mapper == nullptr ? i : mapper->from_chain(i);
        if (ele_ptr->get_out_degree() == 0)
        {
            Pmat.alloc_major(mat_ind, 1);
            Pmat.add_entry(mat_ind, mat_ind, 1.0);
        } else
        {
            Pmat.alloc_major(mat_ind, ele_ptr->get_out_degree());
            for (auto arc : ele_ptr->get_to_arc_list())
            {
                auto dest_ptr = static_cast<ElementType *>(arc.dest_ele);
                uint_t dest_mat_ind = mapper == nullptr ?
                                      dest_ptr->get_index() : mapper->from_chain(dest_ptr->get_index());
                if (arc.val != 0.0)
                {
                    Pmat.add_entry(mat_ind, dest_mat_ind, arc.val / ele_ptr->get_out_sum());
                }
            }
        }

        Pmat.assemble(mat_ind);
    }

    return Pmat;
}


ColSparseM subchain_to_Qmatrix(const Subchain &subchain, ChainIndexMapper *mapper = nullptr);
ColSparseM subchain_to_Pmatrix(const Subchain &subchain, ChainIndexMapper *mapper = nullptr);

