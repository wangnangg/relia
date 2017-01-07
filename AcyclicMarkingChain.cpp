#include "AcyclicMarkingChain.h"

class AcyclicChainElement : public BasicChainElement
{
    double tau_sum;
public:
    AcyclicChainElement(Marking &&mk, uint_t index) : BasicChainElement(std::move(mk), index), tau_sum(0)
    {}

    void accu_tau(double val)
    { tau_sum += val; }

    double compute_tau() const
    { return tau_sum / get_out_sum(); }
};


typedef MarkingChain<AcyclicChainElement> AcyclicMarkingChain;

std::vector<uint_t> topology_sort(const AcyclicMarkingChain &chain)
{
    std::vector<uint_t> in_degree_counter(chain.size(), 0);
    for (uint_t i = 0; i < chain.size(); i++)
    {
        for (ChainArc arc : chain[i]->get_to_arc_list())
        {
            auto dest_ele = static_cast<AcyclicChainElement *>(arc.dest_ele);
            if (chain.contains(dest_ele))
            {
                in_degree_counter[dest_ele->get_index()] += 1;
            }
        }
    }

    std::vector<uint_t> result;
    result.reserve(chain.size());
    for (uint_t i = 0; i < chain.size(); i++)
    {
        if (in_degree_counter[i] == 0)
        {
            result.push_back(i);
        }
    }

    for (uint_t i = 0; i < result.size(); i++)
    {
        auto ele = chain[result[i]];
        for (ChainArc arc : ele->get_to_arc_list())
        {
            auto dest_ele = static_cast<AcyclicChainElement *>(arc.dest_ele);
            if (chain.contains(dest_ele))
            {
                uint_t index = dest_ele->get_index();
                in_degree_counter[index] -= 1;
                if (in_degree_counter[index] == 0)
                {
                    result.push_back(index);
                }
            }
        }
    }
    if (result.size() != chain.size())
    {
        throw -1;
    }
    return result;
}

double eval_level_tau(const AcyclicMarkingChain &chain, const std::vector<uint_t> &eval_order)
{
    double level_tau = 0;
    for (uint_t i = 0; i < eval_order.size(); i++)
    {
        AcyclicChainElement *ele = chain[eval_order[i]];
        if (ele->get_marking().type == Marking::Absorbing)
        {
            continue;
        }
        double ele_tau = ele->compute_tau();
        level_tau += ele_tau;
        for (auto arc : ele->get_to_arc_list())
        {
            AcyclicChainElement *dest_ele = static_cast<AcyclicChainElement *>(arc.dest_ele);
            dest_ele->accu_tau(arc.val * ele_tau);
        }
    }
    return level_tau;
}

MarkingChain<AcyclicChainElement>
generate_next_level_chain(const PetriNet &petri_net, const MarkingChain<AcyclicChainElement> &current_level,
                          const IterStopCondition van_chain_stop_condition)
{
    MarkingChain<AcyclicChainElement> next_level;
    std::vector<const MarkingChain<AcyclicChainElement> *> chain_list{&current_level, &next_level};
    for (uint_t i = 0; i < current_level.size(); i++)
    {
        explore_tangible_marking(petri_net, chain_list, next_level, current_level[i], van_chain_stop_condition);
    }
    return next_level;
}

double compute_acyclic_mtta(const PetriNet &petri_net, const IterStopCondition van_chain_stop_condition)
{
    MarkingChain<AcyclicChainElement> current_level;
    std::vector<const MarkingChain<AcyclicChainElement> *> chain_list{&current_level};
    if (petri_net.get_init_marking().type == Marking::Vanishing)
    {
        auto van_chain = recursive_explore_vanishing_marking(petri_net, chain_list, current_level,
                                                             petri_net.get_init_marking().clone());
        IterStopCondition stop_cond = van_chain_stop_condition;
        auto arc_list = collapse_vanishing_chain(van_chain, stop_cond);
        for (auto arc : arc_list)
        {
            auto ele_ptr = static_cast<AcyclicChainElement *>(arc.dest_ele);
            ele_ptr->accu_tau(arc.val);
        }
    } else
    {
        current_level.add_marking(petri_net.get_init_marking().clone());
        current_level[0]->accu_tau(1.0);
    }

    double tau = 0;
    MarkingChain<AcyclicChainElement> next_level;
    do
    {
        next_level = generate_next_level_chain(petri_net, current_level, van_chain_stop_condition);
        auto eval_order = topology_sort(current_level);
        tau += eval_level_tau(current_level, eval_order);
        std::swap(current_level, next_level);
    } while (current_level.size() != 0);
    return tau;
}
