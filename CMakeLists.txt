cmake_minimum_required(VERSION 3.16)
project(ContigConsensus)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(ContigConsensus
        Contig.cpp
        Contig.h
        main.cpp
        Match.cpp
        Match.h)
