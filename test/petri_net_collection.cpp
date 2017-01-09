#include "petri_net_collection.h"

void display(std::pair<MarkingChain<BasicChainElement>, MarkingChainSparseState> pair)
{
    const auto &init_state = pair.second;
    std::cout << "init:";
    for (const auto &init_prob : init_state)
    {
        std::cout << "(id:" << init_prob.ele << ", " << "prob:" << init_prob.prob << ") ";
    }
    std::cout << std::endl;
    const auto &tan_chain = pair.first;
    for (uint_t i = 0; i < tan_chain.size(); i++)
    {
        std::cout << "{" << std::endl;
        std::cout << "id:" << tan_chain[i] << std::endl;
        std::cout << "marking: ";
        for (auto m : tan_chain[i]->get_marking().token_list)
        {
            std::cout << m << ", ";
        }
        std::cout << std::endl;
        std::cout << "arc:[" << std::endl;
        for (ChainArc arc : tan_chain[i]->get_to_arc_list())
        {
            std::cout << "\t(id:" << arc.dest_ele << ", val:" << arc.val << ")" << std::endl;
        }
        std::cout << "]}" << std::endl << std::endl;
    }
}

PetriNet trivial_petri_net()
{
    PetriNet builder(2);
    auto t1 = builder.add_transition(TransType::Exp, true, 0.1, 1);
    auto t2 = builder.add_transition(TransType::Exp, true, 0.2, 1);
    builder.add_arc(ArcType::In, t1, 0, 1);
    builder.add_arc(ArcType::Out, t1, 1, 1);
    builder.add_arc(ArcType::In, t2, 1, 1);
    builder.add_arc(ArcType::Out, t2, 0, 1);
    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

PetriNet molloys_petri_net()
{
    PetriNet builder(5);
    uint_t t0 = builder.add_transition(TransType::Exp, true, 1.0, 1);
    uint_t t1 = builder.add_transition(TransType::Exp, true, 3.0, 1);
    uint_t t2 = builder.add_transition(TransType::Exp, true, 7.0, 1);
    uint_t t3 = builder.add_transition(TransType::Exp, true, 9.0, 1);
    uint_t t4 = builder.add_transition(TransType::Exp, true, 5.0, 1);

    builder.add_arc(In, t0, 0, 1);
    builder.add_arc(Out, t0, 1, 1);
    builder.add_arc(Out, t0, 2, 1);
    builder.add_arc(In, t1, 1, 1);
    builder.add_arc(Out, t1, 3, 1);
    builder.add_arc(In, t2, 2, 1);
    builder.add_arc(Out, t2, 4, 1);
    builder.add_arc(In, t3, 3, 1);
    builder.add_arc(Out, t3, 1, 1);
    builder.add_arc(In, t4, 3, 1);
    builder.add_arc(In, t4, 4, 1);
    builder.add_arc(Out, t4, 0, 1);

    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

void add_trans(PetriNet &builder, TransType type, double val, std::vector<uint_t> in_place_list,
               std::vector<uint_t> out_place_list)
{
    auto t = builder.add_transition(type, true, val, 1);
    for (auto p : in_place_list)
    {
        builder.add_arc(ArcType::In, t, p, 1);
    }
    for (auto p : out_place_list)
    {
        builder.add_arc(ArcType::Out, t, p, 1);
    }
}

PetriNet acyclic_imme_petri_net()
{
    PetriNet builder(7);
    add_trans(builder, Imme, 1.0, {0}, {1});
    add_trans(builder, Imme, 0.6, {1}, {2});
    add_trans(builder, Imme, 0.4, {1}, {3});
    add_trans(builder, Imme, 0.2, {3}, {4});
    add_trans(builder, Imme, 0.3, {3}, {5});
    add_trans(builder, Imme, 0.5, {3}, {6});
    add_trans(builder, Imme, 1.0, {2}, {4});
    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

PetriNet acyclic_imme_petri_net2()
{
    PetriNet builder(8);
    add_trans(builder, Exp, 1.0, {7}, {0});
    add_trans(builder, Imme, 1.0, {0}, {1});
    add_trans(builder, Imme, 0.6, {1}, {2});
    add_trans(builder, Imme, 0.4, {1}, {3});
    add_trans(builder, Imme, 0.2, {3}, {4});
    add_trans(builder, Imme, 0.3, {3}, {5});
    add_trans(builder, Imme, 0.5, {3}, {6});
    add_trans(builder, Imme, 1.0, {2}, {4});
    builder.set_init_token(7, 1);
    builder.finalize();
    return builder;
}

PetriNet cyclic_imme_petri_net()
{
    PetriNet builder(7);
    add_trans(builder, Exp, 1.0, {0}, {1});
    add_trans(builder, Imme, 0.5, {1}, {4});
    add_trans(builder, Imme, 0.5, {1}, {2});
    add_trans(builder, Imme, 0.5, {2}, {5});
    add_trans(builder, Imme, 0.5, {2}, {3});
    add_trans(builder, Imme, 0.5, {3}, {6});
    add_trans(builder, Imme, 0.5, {3}, {1});
    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

PetriNet cyclic_imme_petri_net2()
{
    PetriNet builder(7);
    add_trans(builder, Exp, 1.0, {0}, {1});
    add_trans(builder, Imme, 0.5, {1}, {4});
    add_trans(builder, Imme, 0.5, {1}, {2});
    add_trans(builder, Imme, 0.5, {2}, {5});
    add_trans(builder, Imme, 0.5, {2}, {3});
    add_trans(builder, Imme, 0.5, {3}, {6});
    add_trans(builder, Imme, 0.5, {3}, {1});
    add_trans(builder, Imme, 0.5, {3}, {2});
    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

PetriNet acyclic_trivial_petri_net()
{
    PetriNet builder(2);
    auto t1 = builder.add_transition(Exp, true, 0.1, 1);
    builder.add_arc(In, t1, 0, 1);
    builder.add_arc(Out, t1, 1, 1);
    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

PetriNet acyclic_petri_net()
{
    PetriNet builder(7);
    add_trans(builder, Exp, 1.0, {0}, {1});
    add_trans(builder, Exp, 0.6, {1}, {2});
    add_trans(builder, Exp, 0.4, {1}, {3});
    add_trans(builder, Exp, 0.2, {3}, {4});
    add_trans(builder, Exp, 0.3, {3}, {5});
    add_trans(builder, Exp, 0.5, {3}, {6});
    add_trans(builder, Exp, 1.0, {2}, {4});
    builder.set_init_token(0, 1);
    builder.finalize();
    return builder;
}

uint_t numG(PetriNetContext *context)
{
    return context->marking->token_list[0];
}

uint_t numB(PetriNetContext *context)
{
    return context->marking->token_list[1];
}

uint_t numE(PetriNetContext *context)
{
    return context->marking->token_list[2];
}

uint_t numF(PetriNetContext *context)
{
    return context->marking->token_list[3];
}

uint_t numD(PetriNetContext *context)
{
    return context->marking->token_list[4];
}

uint_t numHalt(PetriNetContext *context)
{
    return context->marking->token_list[5];
}

static uint_t num_nodes = 20;
static uint_t m = 4;

PetriNet securityCPS_petri_net()
{
    PetriNet builder;
    builder.set_init_token(0, num_nodes);
    builder.set_init_token(4, 1);
    auto tge = builder.add_transition(Exp, true, static_cast<ConstOrVar<double>::CallBack>([](PetriNetContext *context)
    {
        return 0.5 * numG(context);
    }), 0);
    auto tgb = builder.add_transition(Exp, true, static_cast<ConstOrVar<double>::CallBack>([](PetriNetContext *context)
    {
        return 1.0 * numG(context);
    }), 0);
    auto tbe = builder.add_transition(Exp, true,
                                      static_cast<ConstOrVar<double>::CallBack>([](PetriNetContext *context) -> double
                                      {
                                          return 5.0 * numB(context);
                                      }), 0);
    auto tbf = builder.add_transition(Exp, true, static_cast<ConstOrVar<double>::CallBack>([](PetriNetContext *context)
    {
        return 0.0208333333333 * numB(context);
    }), 0);
    auto td = builder.add_transition(Exp, true, 0.0057142095238, 0);
    auto tflush = builder.add_transition(Imme,
                                         static_cast<ConstOrVar<bool>::CallBack>([](PetriNetContext *context) -> bool
                                         {
                                             if (((numG(context) + numB(context) < m) || (numF(context) >= 1) ||
                                                  (numD(context) == 0) ||
                                                  (numG(context) <= 2 * numB(context))) && (numHalt(context) == 0))
                                             {
                                                 return true;
                                             } else
                                             {
                                                 return false;
                                             }
                                         }), 1.0, 0);

    builder.add_arc(In, tgb, 0, 1);
    builder.add_arc(In, tge, 0, 1);
    builder.add_arc(In, tbe, 1, 1);
    builder.add_arc(In, tbf, 1, 1);
    builder.add_arc(In, td, 4, 1);

    builder.add_arc(In, tflush, 0, numG);
    builder.add_arc(In, tflush, 1, numB);
    builder.add_arc(In, tflush, 2, numE);
    builder.add_arc(In, tflush, 3, numF);
    builder.add_arc(In, tflush, 4, numD);

    builder.add_arc(Out, tgb, 1, 1);
    builder.add_arc(Out, tge, 2, 1);
    builder.add_arc(Out, tbe, 2, 1);
    builder.add_arc(Out, tbf, 3, 1);
    builder.add_arc(Out, tflush, 5, 1);
    builder.finalize();
    return builder;
}
