//
// Created by wangnan on 1/10/17.
//
#include "Export.h"
#include <utility>
#include <vector>
#include <limits>

Graph export_petri_net(const PetriNet &petri_net)
{
    Graph g;
    auto& node_list = g.node_list;
    node_list.reserve(petri_net.trans_count() + petri_net.place_count());
    auto& edge_list = g.edge_list;
    edge_list.reserve(petri_net.trans_count() + petri_net.place_count());
    for (uint_t i = 0; i < petri_net.place_count(); i++)
    {
        node_list.emplace_back(i, 0, NAN);
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
    return g;
};


