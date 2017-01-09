#pragma once

#include "MarkingChain.h"
#include <iostream>

void display(std::pair<MarkingChain<BasicChainElement>, MarkingChainSparseState> pair);

PetriNet trivial_petri_net();


PetriNet molloys_petri_net();

void add_trans(PetriNet &builder, TransType type, double val,
               std::vector<uint_t> in_place_list, std::vector<uint_t> out_place_list);

PetriNet acyclic_imme_petri_net();

PetriNet acyclic_imme_petri_net2();

PetriNet cyclic_imme_petri_net();

PetriNet cyclic_imme_petri_net2();


PetriNet acyclic_trivial_petri_net();


PetriNet acyclic_petri_net();


PetriNet securityCPS_petri_net();


PetriNet mixed_class_petri_net();
PetriNet mixed_class_petri_net2();
