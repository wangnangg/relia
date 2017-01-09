#include<gtest/gtest.h>
#include "PetriNet.h"

TEST(petri_net, two_place)
{
    PetriNet builder;
    auto t1 = builder.add_transition(Exp, true, 0.1, 1);
    auto t2 = builder.add_transition(Exp, true, 0.2, 1);
    builder.add_arc(In, t1, 0, 1);
    builder.add_arc(Out, t1, 1, 1);
    builder.add_arc(In, t2, 1, 1);
    builder.add_arc(Out, t2, 0, 1);
    builder.set_init_token(0, 1);
    builder.finalize();
    PetriNet petri_net(std::move(builder));
    auto vec = petri_net.next_markings(petri_net.get_init_marking());
    ASSERT_EQ(vec.size(), 1);
    ASSERT_EQ(vec[0].second, 0.1);
    ASSERT_EQ(vec[0].first.token_list[0], 0);
    ASSERT_EQ(vec[0].first.token_list[1], 1);
    ASSERT_EQ(vec[0].first.f_enabled_trans_ind, 1);
    ASSERT_EQ(vec[0].first.type, Marking::Tangible);

    auto vec2 = petri_net.next_markings(vec[0].first);
    ASSERT_EQ(vec2.size(), 1);
    ASSERT_EQ(vec2[0].second, 0.2);
    auto &mk = vec2[0].first;
    ASSERT_EQ(mk.token_list[0], 1);
    ASSERT_EQ(mk.token_list[1], 0);
    ASSERT_EQ(mk.f_enabled_trans_ind, 0);
    ASSERT_EQ(mk.type, Marking::Tangible);
}


TEST(petri_net, sort_transition)
{
    PetriNet builder(2);
    auto t1 = builder.add_transition(Exp, true, 0.1, 2);
    auto t2 = builder.add_transition(Imme, true, 0.2, 4);
    auto t3 = builder.add_transition(Exp, true, 0.3, 1);
    auto t4 = builder.add_transition(Imme, true, 0.4, 3);
    builder.finalize();
    PetriNet petri_net(std::move(builder));
    PetriNetContext context{&petri_net, &petri_net.get_init_marking()};
    ASSERT_EQ(petri_net.trans_list[0].type, Imme);
    ASSERT_EQ(petri_net.trans_list[1].type, Imme);
    ASSERT_EQ(petri_net.trans_list[2].type, Exp);
    ASSERT_EQ(petri_net.trans_list[3].type, Exp);
    ASSERT_EQ(petri_net.trans_list[0].param.get_val(&context), 0.4);
    ASSERT_EQ(petri_net.trans_list[1].param.get_val(&context), 0.2);
    ASSERT_EQ(petri_net.trans_list[2].param.get_val(&context), 0.3);
    ASSERT_EQ(petri_net.trans_list[3].param.get_val(&context), 0.1);

}
