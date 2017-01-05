//
// Created by wangnan on 1/4/17.
//
#include <gtest/gtest.h>
#include "PetriNetSolution.h"
#include "petri_net_collection.h"

TEST(Solution, reward_func)
{
    PetriNetSolution solution;
    solution.petri_net = trivial_petri_net();
    uint_t index = solution.add_inst_reward_func(static_cast<PetriNetSolution::RewardFuncType>(
                                                         [](PetriNetContext *context)
                                                         {
                                                             return 1.0;
                                                         }));
    auto stop = solution.solve_steady_state();
    ASSERT_TRUE(stop.is_precision_reached());
    double val = solution.get_inst_reward(index);
    ASSERT_DOUBLE_EQ(val, 1.0);

    stop = solution.solve_steady_state();
    ASSERT_TRUE(stop.is_precision_reached());
    val = solution.get_inst_reward(index);
    ASSERT_DOUBLE_EQ(val, 1.0);

}
