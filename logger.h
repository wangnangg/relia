//
// Created by wangnan on 1/6/17.
//

#ifndef RELIAPY_LOGGER_H
#define RELIAPY_LOGGER_H

#include <iostream>
#include <fstream>
extern std::ofstream __log_file;

//TODO: multi-level log

#ifdef _DEBUG
#define LOG_LEVEL_2
#endif
#ifdef LOG_LEVEL_1
#ifdef USE_STDOUT
#define LOG_FILE std::cout
#define LOG_INIT
#else
#define LOG_FILE __log_file
#define LOG_INIT std::ofstream __log_file("relia.log", std::ios::out | std::ios::trunc)
#endif
#define LOG1(s) LOG_FILE << s << std::endl
#define LOG2(s)
#elif defined(LOG_LEVEL_2)
#ifdef USE_STDOUT
#define LOG_FILE std::cout
#define LOG_INIT
#else
#define LOG_FILE __log_file
#define LOG_INIT std::ofstream __log_file("relia.log", std::ios::out | std::ios::trunc)
#endif
#define LOG1(s) LOG_FILE << s << std::endl
#define LOG2(s) LOG_FILE << s << std::endl
#else
#define LOG1(s)
#define LOG2(s)
#define LOG_INIT
#endif

#endif //RELIAPY_LOGGER_H
