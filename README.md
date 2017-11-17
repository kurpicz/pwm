# Parallel Wavelet Tree and Wavelet Matrix Construction

[![Build Status](https://travis-ci.org/kurpicz/pwm.svg?branch=master)](https://travis-ci.org/kurpicz/pwm)

## What is it?
We implemented different simple but very fast sequential and parallel wavelet tree (WT) and wavelet matrix (WM) construction algorithms.
The are based on two ideas, namely *prefix counting* (*pc*) and *prefix sorting* (*ps*).

Using these ideas, we have implemented:
1. a sequential version of the algorithms (*wx_ps* and *wx_ps*),
2. a parallel version of the algorithms (*wx_pps* and *wx_pps*) and
3. a parallel version of the algorithms using *domain decomposition* (*wx_dd_ps* and *wx_dd_ps*).
4. a version of *pc* that scans the text just twice and fills all levels but the first level at once (indicated by an *_ss*-suffix).

In addition, there are also a naive construction algorithms (*wx_naive*).
Replace _wx_ with either _wt_ or _wm_ for the corresponding WT- or WM-construction algorithm.

A detailed description and benchmarks of the implemented algorithms can be found in this [arXiv preprint](https://arxiv.org/abs/1702.07578).

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

## Contributors
- Florian Kurpicz (author)
- Johannes Fischer
- Marvin Löbel
