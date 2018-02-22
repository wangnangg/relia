#ifndef RELIA_PETRINETSOLUTION_H
#define RELIA_PETRINETSOLUTION_H

#include "PetriNet.h"
#include "MarkingChain.h"
#include "Matrix.h"
#include "AcyclicMarkingChain.h"
#include "easylogging++.h"
struct Option
{
    enum SSMethod
    {
        SS_Auto = 0,
        SS_SOR = 1,
        SS_Power = 2,
        SS_Divide = 3
    } steady_state_method;
    enum TSMethod
    {
        TS_Auto = 0
    } transient_state_method;
    double sor_omega;
    uint_t max_interation;
    uint_t check_interval;
    double precision;

    Option(const Option &) = delete;

    Option() = default;
};

struct MarkingVal
{
    Marking marking;
    double prob;
    double stay_time;

    MarkingVal(Marking marking, double prob, double stay_time) : marking(std::move(marking)), prob(prob),
                                                                 stay_time(stay_time)
    {}
};

std::vector<MarkingVal> solve_ss_power(const PetriNet &petri_net,
                                       const IterStopCondition &stop_condition);

std::vector<MarkingVal> solve_ss_sor(const PetriNet &petri_net,
                                     const IterStopCondition &stop_condition,
                                     double omega);

std::vector<MarkingVal> solve_ss_auto(const PetriNet &petri_net,
                                      const IterStopCondition &stop_condition);

std::vector<MarkingVal> solve_ss_divide(const PetriNet &petri_net, const IterStopCondition &stop_condition);


class PetriNetSolution
{
public:
    typedef std::function<double(PetriNetContext *)> RewardFuncType;
private:
    std::vector<double> inst_reward;
    std::vector<double> cum_reward;
    std::vector<RewardFuncType> inst_reward_func;
    std::vector<RewardFuncType> cum_reward_func;
    Option _option;
    const Option &option;
    std::map<std::string, uint_t> chain_generate_time;
    std::pair<MarkingChain<SubchainElement>, MarkingChainSparseState> mkchain_subchain;
    std::pair<MarkingChain<BasicChainElement>, MarkingChainSparseState> mkchain_basic;
public:
    PetriNet petri_net;
    void *tag;

    PetriNetSolution() : option(_option), tag(nullptr)
    {
        _option.steady_state_method = Option::SS_Auto;
        _option.transient_state_method = Option::TS_Auto;
        _option.max_interation = 100000;
        _option.precision = 1e-6;
        _option.sor_omega = 1.0;
        _option.check_interval = 10;
        el::Configurations defaultConf;
        defaultConf.setToDefault();
        defaultConf.setGlobally(
                el::ConfigurationType::Enabled, "false");
        el::Loggers::reconfigureAllLoggers(defaultConf);
        chain_generate_time["subchain"] = 0;
        chain_generate_time["basic"] = 0;
    }

    void update_chain(const std::string& name)
    {
        if(chain_generate_time[name] == petri_net.get_last_change_time())
        {
            return;
        }
        LOG(INFO) << "updating chain:" << name;
        chain_generate_time[name] = petri_net.get_last_change_time();
        IterStopCondition stop_condition(option.max_interation, option.precision, option.check_interval);
        if(name == "subchain")
        {
            mkchain_subchain = generate_marking_chain<SubchainElement>(petri_net, stop_condition);
        } else if (name == "basic")
        {
            mkchain_basic = generate_marking_chain<BasicChainElement>(petri_net, stop_condition);
        }

    }



    //option
    void set_ss_method(Option::SSMethod method)
    {
        _option.steady_state_method = method;
    }

    void set_ts_method(Option::TSMethod method)
    {
        _option.transient_state_method = method;
    }

    void set_sor_omega(double omega)
    {
        _option.sor_omega = omega;
    }

    void set_max_iter(uint_t iter)
    {
        _option.max_interation = iter;
    }

    void set_precision(double prec)
    {
        _option.precision = prec;
    }

    void set_halt_condition(std::function<bool(PetriNetContext *)> func)
    {
        petri_net.set_halt_condition(func);
    }

    void config_logger(const char *filename)
    {
        el::Configurations conf(filename);
        el::Loggers::reconfigureAllLoggers(conf);
    }


    template <typename ReturnType>
    void check_callback_function(std::function<ReturnType(PetriNetContext *)> func)
    {
        return;
        PetriNetContext context = {&petri_net, &petri_net.init_marking};
        func(&context);
    }

    uint_t add_inst_reward_func(RewardFuncType reward_func)
    {
        check_callback_function(reward_func);
        inst_reward_func.push_back(reward_func);
        inst_reward.push_back(0.0);
        return inst_reward_func.size() - 1;
    }

    uint_t add_cum_reward_func(RewardFuncType reward_func)
    {
        check_callback_function(reward_func);
        cum_reward_func.push_back(reward_func);
        cum_reward.push_back(0.0);
        return cum_reward_func.size() - 1;
    }


    double get_inst_reward(uint_t reward_index)
    {
        return inst_reward[reward_index];
    }

    double get_cum_reward(uint_t reward_index)
    {
        return cum_reward[reward_index];
    }

    void solve_steady_state();

    void solve_transient_state(double time);

    template<typename T>
    std::pair<MarkingChain<T>, MarkingChainSparseState> gen_marking_chain()
    {
        petri_net.finalize();
        IterStopCondition stop_condition(option.max_interation, option.precision, option.check_interval);
        auto chain_pair = generate_marking_chain<T>(petri_net, stop_condition);
        return chain_pair;
    }

    double get_acyclic_mtta()
    {
        IterStopCondition stop_condition(option.max_interation, option.precision, option.check_interval);
        return compute_acyclic_mtta(petri_net, stop_condition);
    }

    void bind_marking_modifier(const std::function<Marking (PetriNetContext *)>& modifier)
    {
        LOG(INFO) << "binding modifier";
        petri_net.bind_modifier(modifier);
    }

private:
    void update_reward(const std::vector<MarkingVal> &result);

    void clear_reward()
    {
        std::fill(inst_reward.begin(), inst_reward.end(), 0.0);
        std::fill(cum_reward.begin(), cum_reward.end(), 0.0);
    }
};

#endif //RELIA_PETRINETSOLUTION_H
