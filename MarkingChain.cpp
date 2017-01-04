#include "MarkingChain.h"


IterStopCondition van_chain_solve_stop_condition(1000, 1e-9);

ChainMatrixMapper index_tangible_chain(const MarkingChain<BasicChainElement> &chain)
{
    ChainMatrixMapper mapper(chain.size());
    for (uint_t i = 0; i < chain.size(); i++)
    {
        mapper.set_map(i, i);
    }
    return mapper;
}

ColSparseM markingchain_to_Qmatrix(const MarkingChain<BasicChainElement> &chain)
{
    auto mapper = index_tangible_chain(chain);
    uint_t dim = chain.size();
    ColSparseM Qmat(dim);
    for (uint_t i = 0; i < dim; i++)
    {
        auto ele_ptr = chain[i];
        uint_t mat_ind = mapper.chain2mat(i);
        Qmat.alloc_major(mat_ind, ele_ptr->get_out_degree() + 1);
        for (auto arc : ele_ptr->get_to_arc_list())
        {
            auto dest_ptr = static_cast<BasicChainElement *>(arc.dest_ele);
            uint_t dest_mat_ind = mapper.chain2mat(dest_ptr->get_index());
            if (mat_ind != dest_mat_ind && arc.val != 0.0)
            {
                Qmat.add_entry(mat_ind, dest_mat_ind, arc.val);
            }
        }
        if (ele_ptr->get_out_sum() != 0.0)
        {
            Qmat.add_entry(mat_ind, mat_ind, -ele_ptr->get_out_sum());
        }
        Qmat.assemble(mat_ind);
    }
    return Qmat;
}

