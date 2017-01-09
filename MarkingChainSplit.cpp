//
// Created by wangnan on 1/8/17.
//
#include "Type.h"
#include "MarkingChain.h"
#include "MarkingChainSplit.h"

struct SCCVertex
{
    uint_t index;
    uint_t lowlink;
    bool onstack;
    SCCVertex() : index(0), lowlink(0), onstack(false)
    {}
    SCCVertex(const SCCVertex&) = delete;
};

//TODO: A non-recursive version.
void recursive_visit(const MarkingChain<ChainElement> &chain,
                     std::vector<SCCVertex>& vertex_info,
                     uint_t v_ind,
                     uint_t& visit_index,
                     std::vector<uint_t>& scc_stack,
                     std::vector<Subchain>& subchain_list)
{
    auto& v = vertex_info[v_ind];
    v.index = visit_index;
    v.lowlink = visit_index;
    visit_index += 1;
    scc_stack.push_back(v_ind);
    v.onstack = true;
    for(auto arc : chain[v_ind]->get_to_arc_list())
    {
        auto dest_ele = (ChainElement *) arc.dest_ele;
        uint_t dest_ind = dest_ele->get_index();
        auto& w = vertex_info[dest_ind];
        if(w.index == 0)
        {
            recursive_visit(chain, vertex_info, dest_ind, visit_index, scc_stack, subchain_list);
            v.lowlink = std::min(v.lowlink, w.lowlink);
        } else if (w.onstack)
        {
            v.lowlink = std::min(v.lowlink, w.index);
        }
    }

    if(v.lowlink == v.index)
    {
        std::vector<uint_t> ele_ind_list;
        uint_t top_ele;
        do{
            top_ele = scc_stack.back();
            scc_stack.pop_back();
            vertex_info[top_ele].onstack = false;
            ele_ind_list.push_back(top_ele);
        } while(top_ele != v_ind);
        uint_t subchain_index = subchain_list.size();
        subchain_list.emplace_back(subchain_index, chain, ele_ind_list);
    }

}

std::vector<Subchain> split_to_subchains(const MarkingChain<ChainElement> &chain, std::vector<uint_t> start_ind)
{
    std::vector<Subchain> subchain_list;
    std::vector<SCCVertex> vertex_info(chain.size());
    std::vector<uint_t> scc_stack;
    uint_t visit_index = 1;
    for(auto ind : start_ind)
    {
        if(vertex_info[ind].index == 0) //not visit yet
        {
            recursive_visit(chain, vertex_info, ind, visit_index, scc_stack, subchain_list);
        }
    }
    return subchain_list;
}

