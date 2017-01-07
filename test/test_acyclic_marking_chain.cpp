#include<gtest/gtest.h>
#include "AcyclicMarkingChain.h"
#include "petri_net_collection.h"

extern IterStopCondition stop_cond;

TEST(AcyclicMarkingChain, mtta)
{
    auto petri_net = acyclic_trivial_petri_net();
    ASSERT_EQ(compute_acyclic_mtta(petri_net, stop_cond), 10.0);

    petri_net = acyclic_petri_net();
    ASSERT_EQ(compute_acyclic_mtta(petri_net, stop_cond), 3.0);

    petri_net = securityCPS_petri_net();
    ASSERT_DOUBLE_EQ(compute_acyclic_mtta(petri_net, stop_cond), 0.77894902109188757);

}


