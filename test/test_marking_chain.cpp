#include <gtest/gtest.h>
#include "petri_net_collection.h"
#include "helper.h"

TEST(MarkingChain, generate)
{
    display(generate_marking_chain<BasicChainElement>(trivial_petri_net()));
    display(generate_marking_chain<BasicChainElement>(molloys_petri_net()));
    display(generate_marking_chain<BasicChainElement>(acyclic_imme_petri_net()));
    display(generate_marking_chain<BasicChainElement>(acyclic_imme_petri_net2()));
    display(generate_marking_chain<BasicChainElement>(cyclic_imme_petri_net()));
    display(generate_marking_chain<BasicChainElement>(cyclic_imme_petri_net2()));
}

typedef PetriNet(*petri_net_func)();

petri_net_func petri_nets[] = {trivial_petri_net, molloys_petri_net, acyclic_petri_net};

TEST(MarkingChain, solve)
{
    for (auto f : petri_nets)
    {
        std::cout << "start of solve:" << std::endl;
        auto chain_pair = generate_marking_chain<BasicChainElement>(f());
        auto Qmat = markingchain_to_Qmatrix(chain_pair.first);
        auto Qmat_row = to_row_sparse(Qmat);
        Vector x(Qmat_row.dim());
        x.fill(1.0);
        IterStopCondition stop_condition(100, 1e-6);
        sor_method(x, Qmat_row, stop_condition, 1.0);
        x.scale(1.0 / norm1(x));
        display(Qmat);
        display(x);
        std::cout << "used iter:" << stop_condition.get_used_iter() << std::endl << " " << std::endl;
    }

}

