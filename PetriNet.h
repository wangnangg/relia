#pragma once
#include <vector>
#include <utility>
#include <string>
#include "Type.h"


class PetriNet;
struct Marking;
class Transition;

struct PetriNetContext
{
	const PetriNet* petri_net;
	const Marking* marking;
};



template <typename T>
class ConstOrVar
{
public:
	typedef T(*CallBack)(void* context);
private:
	enum Type
	{
		Const,
		Var
	} type;
	union
	{
		T val;
		CallBack func;
	} store;

public:
	ConstOrVar() = default;
	ConstOrVar(T val) :type(Const)
	{
		store.val = val;
	}
	ConstOrVar(CallBack call_back) :type(Var)
	{
		store.func = call_back;
	}
	T get_val(PetriNetContext* context) const
	{
		switch (type)
		{
		case Const:
			return store.val;
		default:
			return store.func(context);
		}
	}
	ConstOrVar& operator=(T val) 
	{
		type = Type::Const;
		store.val = val;
		return *this;
	}
	ConstOrVar& operator=(CallBack call_back) 
	{
		type = Type::Var;
		store.func = call_back;
		return *this;
	}


};


struct Marking
{
	enum Type
	{
		Unknown,
		Vanishing,
		Tangible,
		Absorbing
	};


	Marking(const Marking& marking) = default;
	std::vector<uint_t> token_list;
	uint_t f_enabled_trans_ind;
	Type type;
public:
	Marking() = default;
	Marking(uint_t token_count) : token_list(token_count), f_enabled_trans_ind(0), type(Unknown) {}
	Marking(Marking&& marking) = default;
	Marking& operator=(Marking &&) = default;
	bool operator==(const Marking& other) const { return other.token_list == token_list; };
	Marking clone() const
	{
		return Marking(*this);
	}
};




class Arc
{
public:
	Arc() = default;
	Arc(uint_t p_index, ConstOrVar<uint_t> multi) : p_index(p_index), multi(multi) {}
	uint_t p_index;
	ConstOrVar<uint_t> multi;
};

class Transition
{
public:
	std::vector<Arc> in_arc_list;
	std::vector<Arc> out_arc_list;
	std::vector<Arc> inhibitor_arc_list;
	ConstOrVar<double> param;
	ConstOrVar<bool> guard;
	uint_t priority;
	TransType type;
public:
	Transition() = default;
	Transition(TransType type, ConstOrVar<bool> guard, ConstOrVar<double> param, uint_t priority) :
		in_arc_list(0),
		out_arc_list(0),
		inhibitor_arc_list(0),
		param(param),
		guard(guard),
		priority(priority),
		type(type)
	{}
	bool is_enabled(PetriNetContext *context) const;
	std::pair<Marking, double> fire(PetriNetContext *context) const;
};



class PetriNet
{
	bool finalized = false;
public:
	Marking init_marking;
	std::vector<Transition> trans_list;
	uint_t imme_trans_count;
public:
	PetriNet() = default;
	explicit PetriNet(uint_t place_count) :init_marking(place_count), trans_list(0), imme_trans_count(0) {}
	PetriNet(PetriNet &&) = default;
	PetriNet(const PetriNet &) = delete;

	//builder method:
	uint_t add_transition(TransType type, ConstOrVar<bool> guard, ConstOrVar<double> param, uint_t priority);
	void add_arc(ArcType type, uint_t trans_index, uint_t place_index, ConstOrVar<uint_t> multi);
	void set_init_token(uint_t place_index, uint_t token);
	void finalize();
	//lookup method:
	std::vector<std::pair<Marking, double>> next_markings(const Marking& marking) const;
	const Marking& get_init_marking() const
	{
		return init_marking;
	}

	PetriNet& operator=(PetriNet &&) = default;
private:
	void set_marking_type(Marking& marking) const;
};

