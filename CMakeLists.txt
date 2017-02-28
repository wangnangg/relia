cmake_minimum_required(VERSION 3.5)
project(relia)
add_subdirectory(test)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror  -Wno-unused-variable")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
add_definitions(-DELPP_FEATURE_PERFORMANCE_TRACKING -DAUTO_INITIALIZE_EASYLOGGINGPP)
set(SOURCE_FILES
        helper.h
        helper.cpp
        AcyclicMarkingChain.cpp
        AcyclicMarkingChain.h
        MarkingChain.cpp
        MarkingChain.h
        Matrix.cpp
        Matrix.h
        PetriNet.cpp
        PetriNet.h
        Type.h
        PetriNetSolution.cpp
        PetriNetSolution.h
        MarkingChainSplit.cpp
        MarkingChainSplit.h
        Export.h
        Export.cpp
        MarkingChainSolve.cpp
        MarkingChainSolve.h
        easylogging++.cc easylogging++.h
        MatrixSolve.cpp MatrixSolve.h
        Poisson.cpp Poisson.h
        json.hpp)
add_library(relia ${SOURCE_FILES})
