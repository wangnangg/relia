#include "MarkingChain.h"
#include "MarkingChainSplit.h"
#include "helper.h"
#include "logger.h"


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
                } else
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
                uint_t dest_mat_ind = mapper == nullptr ?
                                      dest_ptr->get_subindex() : mapper->from_chain(dest_ptr->get_subindex());
                Qmat.add_entry(mat_ind, dest_mat_ind, arc.val);
            }
        }
        Qmat.add_entry(mat_ind, mat_ind, -ele_ptr->get_out_sum());
        Qmat.assemble(mat_ind);
    }
    return Qmat;
}


std::vector<ElementProb> solve_tangible_subchain(const Subchain &subchain,
                                                 Vector subchain_init,
                                                 IterStopCondition& stop_condition)
{
    auto Pmat = subchain_to_Pmatrix(subchain);
    Vector sol(Pmat.dim());
    double prob_sum = norm1(subchain_init);
    for(uint_t i=0; i<subchain_init.dim(); i++)
    {
        sol(i) = subchain_init(i);
    }
    power_method(Pmat, sol, stop_condition);
    double sol_sum = 0;
    for(uint_t i=subchain.foreign_element_count(); i<Pmat.dim(); i++)
    {
        sol(i) = std::abs(sol(i));
        sol_sum += sol(i);
    }
    double factor = prob_sum / sol_sum;
    std::vector<ElementProb> result;
    result.reserve(subchain.foreign_element_count());
    for(uint_t i=subchain.foreign_element_count(); i<Pmat.dim(); i++)
    {
        result.emplace_back(subchain[i], sol(i) * factor);
    }
    return result;
}


Vector solve_absorbing_subchain(const Subchain &subchain,
                                double prob_sum,
                                IterStopCondition &stop_condition)
{
    auto Qmat = subchain_to_Qmatrix(subchain);
    LOG2("solve absorbing subchain:\n" << display(Qmat));
    auto Qmat_row = to_row_sparse(Qmat);
    Vector sol(Qmat.dim(), 1.0);
    display(Qmat_row);
    sor_method(sol, Qmat_row, stop_condition, 1.0);
    double norm = norm1(sol);
    double factor = prob_sum / norm;
    for (uint_t i = 0; i < sol.dim(); i++)
    {
        sol(i) = std::abs(sol(i)) * factor;
    }
    LOG2("result:"<<display(sol));
    return sol;
}

inline void feed_init_prob(std::vector<Vector> &subchain_init_prob_vec,
                           const std::vector<ElementProb> &input_prob_list)
{
    for (const ElementProb& init_prob : input_prob_list)
    {
        ChainElement *ele = (ChainElement *) init_prob.ele;
        Vector &v = subchain_init_prob_vec[ele->get_subchain_index()];
        v(ele->get_subindex()) = init_prob.prob;
    }
}


Vector solve_marking_chain(const MarkingChain<ChainElement> &chain,
                           const std::vector<ElementProb> &chain_init_vec,
                           const IterStopCondition& stop_condition)
{
    Vector sol(chain.size());
    std::vector<uint_t> start_ind;
    start_ind.reserve(chain_init_vec.size());
    for(auto prob_pair : chain_init_vec)
    {
        auto ele_ptr = (ChainElement*)prob_pair.ele;
        start_ind.push_back(ele_ptr->get_index());
    }
    auto subchain_list = split_to_subchains(chain, start_ind);
    std::vector<Vector> subchain_init_prob_vec;
    subchain_init_prob_vec.reserve(subchain_list.size());
    for (uint_t i = 0; i < subchain_list.size(); i++)
    {
        subchain_init_prob_vec.emplace_back(subchain_list.size());
    }
    feed_init_prob(subchain_init_prob_vec, chain_init_vec);
    for (uint_t i = 0; i < subchain_list.size(); i++)
    {
        Vector &init_vec = subchain_init_prob_vec[i];
        const Subchain &subchain = subchain_list[i];
        if (subchain.is_absorbing())
        {
            IterStopCondition stop_condition_inst = stop_condition;
            Vector sub_sol = solve_absorbing_subchain(subchain, norm1(init_vec), stop_condition_inst);
            for (uint_t local_ind = 0; local_ind < subchain.home_element_count(); local_ind++)
            {
                uint_t global_index = subchain[local_ind]->get_index();
                double prob = sub_sol(local_ind);
                sol(global_index) = prob;
            }
        } else
        {
            IterStopCondition stop_condition_inst = stop_condition;
            auto dest_prob_list = solve_tangible_subchain(subchain, std::move(init_vec), stop_condition_inst);
            feed_init_prob(subchain_init_prob_vec, dest_prob_list);
        }
    }
    return sol;
}

Vector markingchain_init_to_vector(std::vector<ElementProb> init, uint_t dim)
{
    Vector v(dim);
    for (auto e: init)
    {
        auto ele = (BasicChainElement *) e.ele;
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
            auto dest_ele = (ChainElement *) arc.dest_ele;
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
