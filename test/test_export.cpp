//
// Created by wangnan on 17-2-23.
//

#include <gtest/gtest.h>
#include "petri_net_collection.h"
#include "Export.h"

TEST(Export, petri_net_trivial)
{
    auto pn = trivial_petri_net();
    auto pn_json = export_petri_net2json(pn);
    std::cout << std::setw(4) << pn_json << std::endl;
}

TEST(Export, petri_net)
{
    for(auto& pn : petri_nets)
    {
        auto pn_json = export_petri_net2json(pn());
        std::cout << std::setw(4) << pn_json << std::endl;
    }
}


