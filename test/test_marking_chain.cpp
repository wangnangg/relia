#include <gtest/gtest.h>
#include <MarkingChainSplit.h>
#include "petri_net_collection.h"
#include "helper.h"

IterStopCondition stop_cond;

TEST(MarkingChain, generate)
{
    display(generate_marking_chain<BasicChainElement>(trivial_petri_net(), stop_cond));
    display(generate_marking_chain<BasicChainElement>(molloys_petri_net(), stop_cond));
    display(generate_marking_chain<BasicChainElement>(acyclic_imme_petri_net(), stop_cond));
    display(generate_marking_chain<BasicChainElement>(acyclic_imme_petri_net2(), stop_cond));
    display(generate_marking_chain<BasicChainElement>(cyclic_imme_petri_net(), stop_cond));
    display(generate_marking_chain<BasicChainElement>(cyclic_imme_petri_net2(), stop_cond));
}

typedef PetriNet(*petri_net_func)();

petri_net_func petri_nets[] = {trivial_petri_net, molloys_petri_net, acyclic_imme_petri_net,
                               acyclic_imme_petri_net2,
                               cyclic_imme_petri_net,
                               cyclic_imme_petri_net2,
                                acyclic_petri_net,
							   mixed_class_petri_net,
	mixed_class_petri_net2
};

TEST(MarkingChain, solve)
{
    for (auto f : petri_nets)
    {
        std::cout << "start of solve:" << std::endl;
        auto chain_pair = generate_marking_chain<ChainElement>(f(), stop_cond);
        auto& chain = chain_pair.first;
        auto& chain_init = chain_pair.second;
        IterStopCondition stop_condition(1000, 1e-6);
        Vector sol = solve_marking_chain(chain, chain_init, stop_condition);
		std::cout << "solution:" << display(sol) << std::endl;
    }
}

TEST(MarkingChain, split)
{
    for (auto f : petri_nets)
    {
        auto chain_pair = generate_marking_chain<ChainElement>(f(), stop_cond);
        auto& chain = chain_pair.first;
        auto& chain_init = chain_pair.second;
        std::vector<uint_t> start_ind;
        start_ind.reserve(chain_init.size());
        for(auto prob_pair : chain_init)
        {
            auto ele_ptr = (ChainElement*)prob_pair.ele;
            start_ind.push_back(ele_ptr->get_index());
        }
        auto subchain_list = split_to_subchains(chain, start_ind);
        uint_t size = 0;
        for(auto& subchain : subchain_list)
        {
            size += subchain.home_element_count();
        }
        ASSERT_EQ(size, chain.size());
    }
}

