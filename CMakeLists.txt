cmake_minimum_required (VERSION 2.8.11)
project (flow_cutter_pace20)
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -DNDEBUG -march=native -mtune=native -ffast-math")
add_executable (flow_cutter_pace20 src/include_all.cpp)
