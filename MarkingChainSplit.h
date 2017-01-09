//
// Created by wangnan on 1/8/17.
//

#ifndef RELIAPY_MARKINGCHAINSPLIT_H
#define RELIAPY_MARKINGCHAINSPLIT_H
#include "MarkingChain.h"
std::vector<Subchain> split_to_subchains(const MarkingChain<ChainElement> &chain, std::vector<uint_t> start_ind);
#endif //RELIAPY_MARKINGCHAINSPLIT_H
