//
// Created by wangnan on 1/10/17.
//

#ifndef RELIAPY_EXPORT_H
#define RELIAPY_EXPORT_H

#include "PetriNet.h"
#include "MarkingChain.h"
#include "json.hpp"
//TODO: Export marking chain



using json = nlohmann::json;
json export_petri_net2json(const PetriNet &petri_net);
json export_marking2json(const Marking &m);
Marking import_marking_from_json(const json& j);


#endif //RELIAPY_EXPORT_H
