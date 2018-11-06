# Huffman-shaped Wavelet Trees and Wavelet Matrices

Here, we extended our algorithms such that their output is a Huffman-shaped wavelet tree (HWT) or Huffman-shaped wavelet matrix (HWM).
To this end, we compute the Huffman codes for all symbols in our input and compute the bit vectors for the codes.
Hence, the size of the levels of the tree/matrix will be decreasing.
The first level of a HWT/HWM contains the exactly as many bits as the first level of a WT/WM.

All (non-naive) algorithms for Huffman-shaped wavelet trees and matrices are implemented in two settings: the _destructive_ and the _non-destructive_ construction of the wavelet tree or matrix.
We indicate the mode by appending a suffix _d_ to _destructive_ construction algorithms.

We have implemented:
1. TODO

In addition to the algorithms listed above, we have also implemented naive construction algorithms (*wx_naive*).

Replace _wx_ with either _hwt_ or _hwm_ for the corresponding HWT- or HWM-construction algorithm.
