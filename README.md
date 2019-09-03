# Parallel (Shared Memory and External Memory) Wavelet Tree and Wavelet Matrix Construction [![Build Status](https://travis-ci.org/kurpicz/pwm.svg?branch=master)](https://travis-ci.org/kurpicz/pwm)

The wavelet tree and wavelet matrix are compact full text indices that can answer *access*, *rank*, and *select* queries (among others) for a text in time logarithmic in the size of the alphabet.
This project contains multiple wavelet tree and wavelet matrix construction algorithms, all implemented in C++.

## Sequential and Parallel Shared Memory Algorithms
We implemented different simple but very fast sequential and parallel wavelet tree (WT) and wavelet matrix (WM) construction algorithms.
The are based on two ideas, namely *prefix counting* (*pc*) and *prefix sorting* (*ps*).

Using these ideas, we have implemented:
1. a sequential version of the algorithms (*wx_pc* and *wx_ps*),
2. a parallel version of the algorithms (*wx_ppc* and *wx_pps*),
3. a parallel version of the algorithms using *domain decomposition* (*wx_dd_pc* and *wx_dd_ps*), and
4. a version of *pc* that scans the text just twice and fills all levels but the first level at once (indicated by an *_ss*-suffix).

In addition, there are also a naive construction algorithms (*wx_naive*).
Replace _wx_ with either _wt_ or _wm_ for the corresponding WT- or WM-construction algorithm.

A detailed description and benchmarks of the implemented algorithms can be found [here](https://doi.org/10.1137/1.9781611975055.2) ([arXiv preprint](https://arxiv.org/abs/1702.07578)).

    @inproceedings{Fischer2018WT,
      author    = {Johannes Fischer and Florian Kurpicz and Marvin L{\"{o}}bel},
      title     = {Simple, Fast and Lightweight Parallel Wavelet Tree Construction},
      booktitle = {Proceedings of the 20th Workshop on Algorithm Engineering and Experiments ({ALENEX})},
      pages     = {9--20},
      publisher = {{SIAM}},
      year      = {2018}
    }

## Sequential and Parallel (Semi-)External Memory Algorithms
We also implemented (semi-)external memory WT and wavelet WM algorithms.
In the semi-external memory model, we keep all data that requires random access in main memory, whereas all other properties/limitations of the [external memory mode](https://en.wikipedia.org/wiki/External_memory_algorithm) are still effective.

We have implemented the following algorithms:
1. a sequential semi-external version of *wx_pc*,
2. two different semi-external versions of *wx_ps*, where one works in-place,
3. a parallel semi-external version of *wx_ppc*,
4. a sequential fully external version of *wx_ps*, and
5. a parallel fully external version that is similar to *wx_dd_pc*.

## How to get it?
First clone this repository, then build all executables.
```sh
git clone https://github.com/kurpicz/pwm.git
cd pwm
mkdir build
cd build
cmake ..
make
```
This will create an executable `src/benchmark` which allows you to run all our algorithms.
To this end, we provide a simple command line interface.
Constructing WTs and WMs for a text with (all) our algorithms can be done by running `./src/benchmark -f <path_to_file>`.
To get a list of the other available options run `./src/benchmark --help`.

If you want to test the algorithms simply type `make check` in the build directory.
Then, we create WTs and WMs with all our algorithms, use them to reconstruct the text, and tell you if something went wrong.

### Building Other Wavelet Tree Construction Algorithms for Benchmarking
We want to make the comparison of our algorithms with the state of the art as easy as possible.
To this end, we include all wavelet tree construction algorithms that have publicly available code in this [project](/tree/master/external/wavelet_construction).
Since some have dependencies, not all are build by default.
We offer the following options:
1. `PWM_BUILD_COMPETITORS` Build competitors. Default: *enabled*.
2. `PWM_ENABLE_CILK_ALGORITHM` Build parallel competitors that require Cilk (only available in gcc < 8.0.0). Default: *enabled*.
3. `PWM_BUILD_SDSL_WT` Build competitors that depend on the [SDSL](https://github.com/simongog/sdsl-lite), which has to be installed locally. Default: *disabled*.