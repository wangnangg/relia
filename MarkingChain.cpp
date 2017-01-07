#include "MarkingChain.h"


ChainMatrixMapper index_tangible_chain(const MarkingChain<BasicChainElement> &chain)
{
    ChainMatrixMapper mapper(chain.size());
    for (uint_t i = 0; i < chain.size(); i++)
    {
        mapper.set_map(i, i);
    }
    return mapper;
}

std::tuple<ColSparseM, Vector, ChainMatrixMapper> markingchain_to_Qmatrix(
        const MarkingChain<BasicChainElement> &chain,
        const MarkingChainInitState &chain_init)
{
    auto mapper = index_tangible_chain(chain);
    uint_t dim = chain.size();
    ColSparseM Qmat(dim);
    for (uint_t i = 0; i < dim; i++)
    {
        auto ele_ptr = chain[i];
        uint_t mat_ind = mapper.chain2mat(i);
        if (ele_ptr->get_out_degree() != 0)
        {
            Qmat.alloc_major(mat_ind, ele_ptr->get_out_degree() + 1);
            for (auto arc : ele_ptr->get_to_arc_list())
            {
                auto dest_ptr = static_cast<BasicChainElement *>(arc.dest_ele);
                uint_t dest_mat_ind = mapper.chain2mat(dest_ptr->get_index());
                //assume no self loop
                if (arc.val != 0.0)
                {
                    Qmat.add_entry(mat_ind, dest_mat_ind, arc.val);
                }
            }
            Qmat.add_entry(mat_ind, mat_ind, -ele_ptr->get_out_sum());
            Qmat.assemble(mat_ind);
        }
    }

    Vector init(dim);
    for (auto ele_init : chain_init)
    {
        auto dest_ptr = (BasicChainElement *) ele_init.ele;
        uint_t chain_index = dest_ptr->get_index();
        uint_t mat_ind = mapper.chain2mat(chain_index);
        init(mat_ind) = ele_init.prob;
    }

    return std::make_tuple(std::move(Qmat), std::move(init), std::move(mapper));
}

std::tuple<ColSparseM, Vector, ChainMatrixMapper> markingchain_to_Pmatrix(
        const MarkingChain<BasicChainElement> &chain,
        const MarkingChainInitState &chain_init)
{
    auto mapper = index_tangible_chain(chain);
    uint_t dim = chain.size();
    ColSparseM Pmat(dim);
    for (uint_t i = 0; i < dim; i++)
    {
        auto ele_ptr = chain[i];
        uint_t mat_ind = mapper.chain2mat(i);
        if (ele_ptr->get_out_degree() == 0)
        {
            Pmat.alloc_major(mat_ind, 1);
            Pmat.add_entry(mat_ind, mat_ind, 1.0);
        } else
        {
            Pmat.alloc_major(mat_ind, ele_ptr->get_out_degree());
            for (auto arc : ele_ptr->get_to_arc_list())
            {
                auto dest_ptr = static_cast<BasicChainElement *>(arc.dest_ele);
                uint_t dest_mat_ind = mapper.chain2mat(dest_ptr->get_index());
                //assume no self loop
                if (arc.val != 0.0)
                {
                    Pmat.add_entry(mat_ind, dest_mat_ind, arc.val / ele_ptr->get_out_sum());
                }
            }
        }

        Pmat.assemble(mat_ind);
    }

    Vector init(dim);
    for (auto ele_init : chain_init)
    {
        auto dest_ptr = (BasicChainElement *) ele_init.ele;
        uint_t chain_index = dest_ptr->get_index();
        uint_t mat_ind = mapper.chain2mat(chain_index);
        init(mat_ind) = ele_init.prob;
    }

    return std::make_tuple(std::move(Pmat), std::move(init), std::move(mapper));
}

