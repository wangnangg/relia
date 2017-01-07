//
// Created by wangnan on 1/6/17.
//

#ifndef RELIAPY_LOGGER_H
#define RELIAPY_LOGGER_H
#include <iostream>
#include <fstream>

extern std::ofstream __log_file;

#ifdef LOG_LEVEL_1
#define LOG1(s) __log_file << s << std::endl;
#define LOG2(s)
#define LOG_INIT std::ofstream __log_file("relia.log", std::ios::out | std::ios::trunc);
#elif defined(LOG_LEVEL_2)
#define LOG1(s) __log_file << s << std::endl;
#define LOG2(s) __log_file << s << std::endl;
#define LOG_INIT std::ofstream __log_file("relia.log", std::ios::out | std::ios::trunc);
#else
#define LOG1(s)
#define LOG2(s)
#define LOG_INIT
#endif

#endif //RELIAPY_LOGGER_H
