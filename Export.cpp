//
// Created by wangnan on 1/10/17.
//
#include "Export.h"
#include <utility>
#include <vector>

template <typename T>
void unfold_var(const ConstOrVar<T>& in, T& out, bool& is_var)
{
    out = T();
    is_var = true;
    if(in.get_type() == ConstOrVar<T>::Const)
    {
       out = in.get_val(nullptr);
        is_var = false;
    }
}

Edge construct_edge(uint_t src, uint_t dest, uint_t type, ConstOrVar<uint_t> multi)
{
    unsigned multi_val;
    bool is_multi_var;
    unfold_var(multi, multi_val, is_multi_var);
    return Edge(src, dest, type, multi_val, is_multi_var);
}

Node construct_node(uint_t index, uint_t type, ConstOrVar<double> param)
{
    double param_val;
    bool is_param_var;
    unfold_var(param, param_val, is_param_var);
    return Node(index, type, param_val, is_param_var);
}

std::pair<std::vector<Node>, std::vector<Edge>> export_petri_net(const PetriNet &petri_net)
{
    std::vector<Node> node_list;
    node_list.reserve(petri_net.trans_count() + petri_net.place_count());
    std::vector<Edge> edge_list;
    edge_list.reserve(petri_net.trans_count() + petri_net.place_count());
    for (uint_t i = 0; i < petri_net.place_count(); i++)
    {
        node_list.emplace_back(Node(i, 0, 0, false));
    }
    for (const auto &tr : petri_net.trans_list)
    {
        uint_t type = tr.type == Imme ? 1 : 2;
        node_list.push_back(construct_node(tr.index, type, tr.param));
        for (const auto &arc : tr.in_arc_list)
        {
            edge_list.push_back(construct_edge(arc.p_index, tr.index, 0, arc.multi));
        }
        for (const auto &arc : tr.out_arc_list)
        {
            edge_list.push_back(construct_edge(tr.index, arc.p_index, 1, arc.multi));
        }
        for (const auto &arc : tr.inhibitor_arc_list)
        {
            edge_list.push_back(construct_edge(arc.p_index, tr.index, 2, arc.multi));
        }
    }
    return std::make_pair(std::move(node_list), std::move(edge_list));
};


