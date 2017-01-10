#pragma once

#include <exception>
#include <string>
typedef unsigned uint_t;
typedef int int_t;

enum ArcType
{
    In = 0,
    Out = 1,
    Inhibitor = 2
};

enum TransType
{
    Imme = 0,
    Exp = 1
};



class Exception : public std::exception
{
    std::string msg;
public:
    Exception(std::string msg):msg(std::move(msg)) {}
    virtual const char* what() const noexcept
    {
       return msg.c_str();
    }
};

