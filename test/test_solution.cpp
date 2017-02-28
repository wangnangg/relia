//
// Created by wangnan on 1/4/17.
//
#include <gtest/gtest.h>
#include "PetriNetSolution.h"
#include "petri_net_collection.h"
#include "AcyclicMarkingChain.h"

double mtta_func(PetriNetContext* context)
{
	if(context->petri_net->has_enabled_trans(context))
	{
		return 1.0;
	} else //absorbing
	{
		return 0.0;
	}
}


TEST(Solution, reward_func)
{
    PetriNetSolution solution;
    uint_t prob_index = solution.add_inst_reward_func(static_cast<PetriNetSolution::RewardFuncType>(
                                                         [](PetriNetContext *context)
                                                         {
                                                             return 1.0;
               
                                                         }));
	uint_t mtta_index = solution.add_cum_reward_func(mtta_func);
	for(auto f : petri_nets)
	{
		solution.petri_net = f();
		solution.set_ss_method(Option::SS_Divide);
		solution.solve_steady_state();
		double val = solution.get_inst_reward(prob_index);
		ASSERT_NEAR(val, 1.0, 1e-10);
		double mtta = solution.get_cum_reward(mtta_index);
		std::cout << "divide mtta:" << mtta << std::endl;
	}

}

TEST(Solution, CPS)
{
	PetriNetSolution solution;
    uint_t mtta_index = solution.add_cum_reward_func(mtta_func);
	uint_t n= 100;
    uint_t m = n / 4;
    LOG(INFO) << "number of nodes: " << n;
	solution.petri_net = securityCPS_petri_net(n, m);
	solution.set_ss_method(Option::SS_SOR);
	solution.solve_steady_state();
	double mtta1 = solution.get_cum_reward(mtta_index);
	std::cout << "security mtta1:" << mtta1 << std::endl;

	solution.set_ss_method(Option::SS_Divide);
	solution.solve_steady_state();
	double mtta2 = solution.get_cum_reward(mtta_index);
	std::cout << "security mtta2:" << mtta2 << std::endl;

	IterStopCondition stop_condition(10000, 1e-10, 10);
	double mtta_acyclic = compute_acyclic_mtta(securityCPS_petri_net(n, m), stop_condition);
	std::cout << "security mtta_acyclic:" << mtta_acyclic << std::endl;
    ASSERT_NEAR(mtta1, mtta_acyclic, 1e-9);
    ASSERT_NEAR(mtta2, mtta_acyclic, 1e-9);

}