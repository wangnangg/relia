//
// Created by wangnan on 1/10/17.
//

#ifndef RELIAPY_EXPORT_H
#define RELIAPY_EXPORT_H

#include "PetriNet.h"
#include "ExportType.h"
std::pair<std::vector<Node>, std::vector<Edge>> export_petri_net(const PetriNet &petri_net);


#endif //RELIAPY_EXPORT_H
