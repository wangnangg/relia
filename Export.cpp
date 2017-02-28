//
// Created by wangnan on 1/10/17.
//
#include "Export.h"
#include <utility>
#include <vector>
#include <limits>
#include "json.hpp"

json trans_type2json(TransType type)
{
    switch (type)
    {
        case TransType::Imme:
            return json("imme");
        default:
            return json("exp");
    }
}

template <typename T>
json constvar2json(const ConstOrVar<T>& x)
{
    switch(x.get_type())
    {
        case  ConstOrVar<T>::Const:
            return json(x.get_val(nullptr));
        default:
            return json("var");
    }
}

json arc2json(const Arc& arc)
{
    json j = json({});
    j["multi"] = constvar2json(arc.multi);
    j["p_index"] = arc.p_index;
    return j;
}


json transition2json(const Transition& t)
{
    json j = json({});
    j["type"] = trans_type2json(t.type);
    j["guard"] = constvar2json(t.guard);
    j["index"] = t.index;
    j["param"] = constvar2json(t.param);
    j["prio"] = t.priority;
    j["in_arc"] = json::array();
    for(const auto& arc: t.in_arc_list)
    {
        j["in_arc"].push_back(arc2json(arc));
    }
    j["out_arc"] = json::array();
    for(const auto& arc: t.out_arc_list)
    {
        j["out_arc"].push_back(arc2json(arc));
    }
    j["inh_arc"] = json::array();
    for(const auto& arc: t.inhibitor_arc_list)
    {
        j["inh_arc"].push_back(arc2json(arc));
    }
    return j;
}

json export_petri_net2json(const PetriNet &petri_net)
{
    json j = json({});

    j["place_count"] = petri_net.place_count();

    json trans_list = json::array();
    for(uint_t i=0; i<petri_net.trans_count(); i++)
    {
        trans_list.push_back(transition2json(petri_net.get_transition(i)));
    }
    j["transition"] = trans_list;
    return j;
}

json export_marking2json(const Marking &mk)
{
    json j = json::array();
    for(uint_t n : mk.token_list)
    {
        j.push_back(n);
    }
    return j;
}


Marking import_marking_from_json(const json& j)
{
    Marking m(j.size());
    for(uint_t i=0; i<j.size(); i++)
    {
        m[i] = j.at(i);
    }
    return m;
};


