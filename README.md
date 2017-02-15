# Parallel Wavelet Tree and Wavelet Matrix Construction
## What is it?
We implemented different simple but very fast sequential and parallel wavelet matrix and wavelet tree construction algorithms.
The two algorithms are named *pc* and *ps* (short for *prefix counting* and *prefix sorting*).

## How to get it?
First clone this repository, then build all executables.
```sh
git clone https://github.com/kurpicz/pwm.git
cd pwm
git submodule update --init
make
```
Now there are many different executables in `benchmark/bin/wm` and `benchmark/bin/wm` for wavelet matrices and wavelet trees, resp.
For each of the two algorithms and two data structures (matrix and tree), there are three different executables:

1. construct_*[wm|wt]_name*
2. check_*[wm|wt]_name*
3. memory_*[wm|wt]_name*

Here, *construct* builds the data structure 5 (default) times and returns the median of the running times,  *check* builds the data structure naively and checks for correctness and *memory* returns the memory peak of the algorithm during construction.

## Contributors
- Florian Kurpicz (author)
- Johannes Fischer
