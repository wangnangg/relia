cmake_minimum_required(VERSION 3.5)
project(relia)
add_subdirectory(test)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
add_definitions(-DLOG_LEVEL_2)
#add_definitions(-DUSE_STDOUT)
set(SOURCE_FILES
        helper.h
        helper.cpp
        AcyclicMarkingChain.cpp
        AcyclicMarkingChain.h
        Algebra.cpp
        Algebra.h
        MarkingChain.cpp
        MarkingChain.h
        Matrix.cpp
        Matrix.h
        PetriNet.cpp
        PetriNet.h
        Type.h
        PetriNetSolution.cpp
        PetriNetSolution.h logger.h MarkingChainSplit.cpp MarkingChainSplit.h)
add_library(relia ${SOURCE_FILES})
