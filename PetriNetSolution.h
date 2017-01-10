#ifndef RELIA_PETRINETSOLUTION_H
#define RELIA_PETRINETSOLUTION_H

#include "PetriNet.h"
#include "MarkingChain.h"
#include "Matrix.h"
#include "logger.h"

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
    double precision;
    uint_t van_chain_max_iter;
    double van_chain_precision;
};
struct MarkingVal
{
    Marking marking;
    double prob;
    double stay_time;
    MarkingVal(Marking marking, double prob, double stay_time):marking(std::move(marking)),prob(prob), stay_time(stay_time)
    {}
};
std::vector<MarkingVal> solve_ss_power(const PetriNet& petri_net,
                                                    const IterStopCondition &stop_condition);;

std::vector<MarkingVal> solve_ss_sor(const PetriNet& petri_net,
                                                  const IterStopCondition &stop_condition,
                                                  double omega);;

std::vector<MarkingVal> solve_ss_auto(const PetriNet& petri_net,
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
public:
    PetriNet petri_net;
    void *tag;

    PetriNetSolution() : option(_option), tag(nullptr)
    {
        LOG1(__FUNCTION__);
        _option.steady_state_method = Option::SS_Auto;
        _option.transient_state_method = Option::TS_Auto;
        _option.max_interation = 10000;
        _option.precision = 1e-10;
        _option.sor_omega = 1.0;
        _option.van_chain_max_iter = 10000;
        _option.van_chain_precision = 1e-10;
    }

    ~PetriNetSolution()
    {
        LOG1(__FUNCTION__);
    }

    //option
    void set_ss_method(Option::SSMethod method)
    {
        LOG1(__FUNCTION__);
        _option.steady_state_method = method;
    }

    void set_ts_method(Option::TSMethod method)
    {
        LOG1(__FUNCTION__);
        _option.transient_state_method = method;
    }

    void set_sor_omega(double omega)
    {
        LOG1(__FUNCTION__);
        _option.sor_omega = omega;
    }

    void set_max_iter(uint_t iter)
    {
        LOG1(__FUNCTION__);
        _option.max_interation = iter;
    }

    void set_precision(double prec)
    {
        LOG1(__FUNCTION__);
        _option.precision = prec;
    }

    void set_halt_condition(std::function<bool(PetriNetContext*)> func)
    {
        petri_net.set_halt_condition(func);
    }

    uint_t add_inst_reward_func(RewardFuncType reward_func)
    {
        LOG1(__FUNCTION__);
        inst_reward_func.push_back(reward_func);
		inst_reward.push_back(0.0);
        return inst_reward_func.size() - 1;
    }

    uint_t add_cum_reward_func(RewardFuncType reward_func)
    {
        cum_reward_func.push_back(reward_func);
		cum_reward.push_back(0.0);
        return cum_reward_func.size() - 1;
    }


    double get_inst_reward(uint_t reward_index)
    {
        LOG1(__FUNCTION__);
        return inst_reward[reward_index];
    }

    double get_cum_reward(uint_t reward_index)
    {
        LOG1(__FUNCTION__);
        return cum_reward[reward_index];
    }

	void solve_steady_state();

private:
	void update_reward(const std::vector<MarkingVal>& result);

	void clear_reward()
	{
		std::fill(inst_reward.begin(), inst_reward.end(), 0.0);
		std::fill(cum_reward.begin(), cum_reward.end(), 0.0);
	}
};

#endif //RELIA_PETRINETSOLUTION_H
