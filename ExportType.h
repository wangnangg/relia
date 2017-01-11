//
// Created by wangnan on 1/10/17.
//

#ifndef RELIAPY_EXPORTTYPE_H
#define RELIAPY_EXPORTTYPE_H

class Node
{
public:
    unsigned index;
    unsigned type; //0 for place ,1 for imm trans, 2 for exp trans, 3 for chain_element
    double param;
    Node() = default;
    Node(unsigned index, unsigned type, double param):
            index(index),
            type(type),
            param(param) {}
};

class Edge
{
public:
    unsigned src;
    unsigned dest;
    unsigned type; //0 for in, 1 for out, 2 for inhibitor, 3 for chain_arc
    double param; //multiplicity or rate
    Edge() = default;
    Edge(unsigned src, unsigned dest, unsigned type, double param):
            src(src), dest(dest), type(type), param(param) {}
};

class Graph
{
public:
    std::vector<Node> node_list;
    std::vector<Edge> edge_list;
    Graph() = default;
};
#endif //RELIAPY_EXPORTTYPE_H
