//
// Created by wangnan on 1/10/17.
//

#ifndef RELIAPY_EXPORT_H
#define RELIAPY_EXPORT_H

#include "PetriNet.h"
#include "MarkingChain.h"
#include "ExportType.h"

//TODO: Export marking chain


template<typename T>
void unfold_var(const ConstOrVar<T> &in, double &out)
{
    if (in.get_type() == ConstOrVar<T>::Const)
    {
        out = in.get_val(nullptr);
    } else
    {
        out = NAN;
    }
}

template<typename T>
Edge construct_edge(uint_t src, uint_t dest, uint_t type, ConstOrVar<T> param)
{
    double param_val;
    unfold_var(param, param_val);
    return Edge(src, dest, type, param_val);
}

template<typename T>
Node construct_node(uint_t index, uint_t type, ConstOrVar<T> param)
{
    double param_val;
    unfold_var(param, param_val);
    return Node(index, type, param_val);
}

Graph export_petri_net(const PetriNet &petri_net);

template<typename ElementType>
Graph export_marking_chain(const MarkingChain<ElementType> &chain)
{
    Graph g;
    g.node_list.reserve(chain.size());
    g.edge_list.reserve(chain.size());
    for (uint_t i = 0; i < chain.size(); i++)
    {
        auto ele_ptr = chain[i];
        g.node_list.emplace_back(i, 3, NAN);
        for(auto arc : ele_ptr->get_to_arc_list())
        {
            ElementType* dest_ele_ptr = (ElementType*)arc.dest_ele;
            g.edge_list.emplace_back(i, dest_ele_ptr->get_index(), 3, arc.val);
        }
    }
    return g;
}


#endif //RELIAPY_EXPORT_H
