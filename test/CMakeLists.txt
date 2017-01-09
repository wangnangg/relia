cmake_minimum_required(VERSION 3.5)
project(unit_test)
add_subdirectory(googletest)

set(CMAKE_CXX_STANDARD 14)
include_directories(../)
add_definitions(-DLOG_LEVEL_2)
add_definitions(-DUSE_STDOUT)
set(SOURCE_FILES
        petri_net_collection.cpp
        petri_net_collection.h
        test_petri_net.cpp
        test_acyclic_marking_chain.cpp
        test_algebra.cpp
        test_marking_chain.cpp
        test_matrix.cpp
        test_solution.cpp)

add_executable(unit_test ${SOURCE_FILES})
target_link_libraries(unit_test relia gtest gtest_main)
