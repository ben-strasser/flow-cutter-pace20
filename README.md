# FlowCutter PACE 2020 Submission

This repository contains the FlowCutter code submitted to the [PACE 2020](https://pacechallenge.org/2020/td/) heuristic tree depth challenge. 

If you are running a Unix-like system, then getting started is very simple. Just clone the repository and build the programs, as follows:

```bash
git clone https://github.com/ben-strasser/flow-cutter-pace20.git
cd flow-cutter-pace20
./build.sh
```

There are no dependencies beyond a recent GCC. GCC 7.5.0 is recent enough. Older versions might work. Clang should also work but has not been tested by us. Building the code under Windows as Windows-native executable is probably possible but some code level modifications are likely necessary. We have not tested the Windows Subsystem for Linux, however, we expect our program to run under it without problems.

After executing the build script, the root directory of the repository should contain the two binary files `flow_cutter_pace20` and `flow_cutter_parallel_pace20`. `flow_cutter_pace20` is the program that is submitted to PACE2020. It is a sequential program that does not make use of multiple cores. `flow_cutter_parallel_pace20` has the exact same commandline but uses as many processor cores as available.

The program can be used as follows:

```bash
./flow_cutter_pace20 < in.gr > out.tree
```

The program runs until a SIGINT or SIGTERM is received. 

The program supports a few additional commandline options. Use `--help` to get a documentation.

## Publications

The following publications are related to this submission:

* Graph Bisection with Pareto Optimization.\
  Michael Hamann and Ben Strasser.\
  ACM Journal of Experimental Algorithmics

* Graph Bisection with Pareto-Optimization.\
  Michael Hamann and Ben Strasser.\
  Proceedings of the 18th Meeting on Algorithm Engineering and Experiments (ALENEX'16).

* Correspondence between Multilevel Graph Partitions and Tree Decompositions.\
  Michael Hamann and Ben Strasser.\
  MDPI Algorithms.

* Computing Tree Decompositions with FlowCutter - PACE 2017 Submission.\
  Ben Strasser.\
  ArXiv.
  
* Customizable Contraction Hierarchies.\
  Julian Dibbelt, Ben Strasser, and Dorothea Wagner.\
  ACM Journal of Experimental Algorithmics.

* Customizable Contraction Hierarchies.\
  Julian Dibbelt, Ben Strasser, and Dorothea Wagner.\
  Symposium on Experimental Algorithms (SEA).

