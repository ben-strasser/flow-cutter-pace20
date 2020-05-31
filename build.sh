#!/bin/sh
g++ -Wall -std=c++11 -O3 -DNDEBUG -march=native -mtune=native -ffast-math src/include_all.cpp -o flow_cutter_pace20
g++ -Wall -std=c++11 -O3 -DNDEBUG -march=native -mtune=native -ffast-math -DPARALLELIZE -fopenmp src/include_all.cpp -o flow_cutter_parallel_pace20
 
