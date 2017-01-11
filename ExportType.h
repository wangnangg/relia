//
// Created by wangnan on 1/10/17.
//

#ifndef RELIAPY_EXPORTTYPE_H
#define RELIAPY_EXPORTTYPE_H
struct Node
{
    unsigned index;
    unsigned type; //0 for place ,1 for imm trans, 2 for exp trans
    double param;
    bool is_param_var;
    Node() = default;
    Node(unsigned index, unsigned type, double param, bool is_param_var):
            index(index),
            type(type),
            param(param),
            is_param_var(is_param_var){}
};

struct Edge
{
    unsigned src;
    unsigned dest;
    unsigned type; //0 for in, 1 for out, 2 for inhibitor
    unsigned multi;
    bool is_multi_var;
    Edge() = default;
    Edge(unsigned src, unsigned dest, unsigned type, unsigned multi, bool is_multi_var):
            src(src), dest(dest), type(type), multi(multi), is_multi_var(is_multi_var) {}
};
#endif //RELIAPY_EXPORTTYPE_H
