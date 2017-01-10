#include "PetriNet.h"
#include<algorithm>

bool Transition::is_enabled(PetriNetContext *context) const
{
    if (!guard.get_val(context))
    {
        return false;
    }
    for (Arc arc : inhibitor_arc_list)
    {
        auto p_token = context->marking->token_list[arc.p_index];
        if (p_token >= arc.multi.get_val(context))
        {
            return false;
        }
    }
    for (Arc arc : in_arc_list)
    {
        auto p_token = context->marking->token_list[arc.p_index];
        if (p_token < arc.multi.get_val(context))
        {
            return false;
        }
    }
    return true;
}

std::pair<Marking, double> Transition::fire(PetriNetContext *context) const
{

    auto result = std::make_pair(context->marking->clone(), param.get_val(context));
    Marking &mk = result.first;
    mk.type = Marking::Unknown;
    for (Arc arc : in_arc_list)
    {
        mk.token_list[arc.p_index] -= arc.multi.get_val(context);
    }
    for (Arc arc : out_arc_list)
    {
        mk.token_list[arc.p_index] += arc.multi.get_val(context);
    }
    return result;
}

uint_t PetriNet::add_transition(TransType type, ConstOrVar<bool> guard, ConstOrVar<double> param, uint_t priority)
{
    finalized = false;
    if (type == TransType::Imme)
    {
        imme_trans_count += 1;
    }
    uint_t index = trans_list.size();
    trans_list.push_back(Transition(type, guard, param, priority));
    return index;
}

void PetriNet::add_arc(ArcType type, uint_t trans_index, uint_t place_index, ConstOrVar<uint_t> multi)
{
    finalized = false;
    Transition &trans = trans_list[trans_index];
    switch (type)
    {
        case ArcType::In:
            trans.in_arc_list.push_back(Arc(place_index, multi));
            break;
        case ArcType::Out:
            trans.out_arc_list.push_back(Arc(place_index, multi));
            break;
        default:
            trans.inhibitor_arc_list.push_back(Arc(place_index, multi));
            break;
    }
    if (place_index >= init_marking.token_list.size())
    {
        init_marking.token_list.resize(place_index + 1, 0);
    }
}

void PetriNet::set_init_token(uint_t place_index, uint_t token)
{
    finalized = false;
    if (place_index >= init_marking.token_list.size())
    {
        init_marking.token_list.resize(place_index + 1, 0);
    }
    init_marking.token_list[place_index] = token;
}


void PetriNet::finalize()
{
    if (finalized)
    {
        return;
    }
    std::sort(trans_list.begin(), trans_list.end(),
              [](const Transition &t1, const Transition &t2)
              {
                  if (t1.type == TransType::Imme && t2.type != TransType::Imme)
                  {
                      return true;
                  } else if (t1.type != TransType::Imme && t2.type == TransType::Imme)
                  {
                      return false;
                  } else
                  {
                      return t1.priority < t2.priority;
                  }
              }
    );
    set_marking_type(init_marking);
    finalized = true;
}

std::vector<std::pair<Marking, double>> PetriNet::next_markings(const Marking &marking) const
{
    if (!finalized)
    {
        throw 0;
    }
    std::vector<std::pair<Marking, double>> result;
    PetriNetContext context{this, &marking};
    if(marking.type == Marking::Absorbing)
    {
        return result;
    }
    const Transition *f_enabled_trans = &trans_list[marking.f_enabled_trans_ind];
    uint_t prio = f_enabled_trans->priority;
    uint_t last_ind = f_enabled_trans->type == TransType::Imme ? imme_trans_count : trans_list.size();

    for (uint_t i = marking.f_enabled_trans_ind;
         i < last_ind && trans_list[i].priority == prio;
         i++)
    {
        const Transition *t = &trans_list[i];
        if (t->is_enabled(&context))
        {
            auto pair = t->fire(&context);
            set_marking_type(pair.first);
            result.push_back(std::move(pair));
        }
    }
    return result;
}

void PetriNet::set_marking_type(Marking &marking) const
{
    PetriNetContext context{this, &marking};
    if(halt_func && halt_func(&context))
    {
        marking.type = Marking::Absorbing;
        return;
    }
    for (uint_t i = 0; i < trans_list.size(); i++)
    {
        const Transition *t = &trans_list[i];
        if (t->is_enabled(&context))
        {
            switch (t->type)
            {
                case TransType::Imme:
                    marking.type = Marking::Vanishing;
                    break;
                default:
                    marking.type = Marking::Tangible;
                    break;
            }
            marking.f_enabled_trans_ind = i;
            return;
        }
    }
    if (marking.type == Marking::Unknown)
    {
        marking.type = Marking::Absorbing;
    }
}

