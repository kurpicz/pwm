// This code is part of the paper "Parallel Wavelet Tree Construction"
// Copyright (c) 2014 Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>

#include "../getmemory.hpp"
#include "../parallel.hpp"
#include "../IO.hpp"
#include "WT.hpp"
#include "../sequence.hpp"

#include <tlx/cmdline_parser.hpp>

using namespace std;
using namespace benchIO;

//naive linear time rank
intT naiveRank(long* bitmap, long start, long finish, int bit) {
  intT rank = 0; 
  for(long i=start;i<finish;i++) {
    long w = bitmap[i/64];
    if((1 & (w >> (i % 64))) == bit) rank++;
  }
  return rank;
}

//access queries---use for testing only; for performance, need to
//implement a fast rank structure
intT query(pair<WTnode*,long*> R, long index, long sigma,
           [[maybe_unused]] long n) {
  [[maybe_unused]] long levels = max(1,utils::log2Up(sigma));
  WTnode* nodes = R.first;
  long* bitmap = R.second;
  WTnode curr = nodes[0]; 
  long min = 0, max = 1 << utils::log2Up(sigma);

  //loop until no children nodes left
  while(max-min > 2 && (curr.leftChild != UINT_T_MAX || curr.rightChild != UINT_T_MAX)) {
    long offset = curr.bitmapPtr;
    long w = bitmap[(index+offset)/64];
    int bit = (w >> ((index+offset) % 64)) & 1;
    if(bit) {
      //subtract number of zeros
      index -= naiveRank(bitmap, offset, offset+index, 0);
      min += (long) 1 << (utils::log2Up(max-min)-1);
      curr = nodes[curr.rightChild];
    }
    else {
      index -= naiveRank(bitmap, offset, offset+index, 1);
      max = min + ((long) 1 << (utils::log2Up(max-min)-1)); 
      curr = nodes[curr.leftChild];
    }
  }
  long offset = curr.bitmapPtr;
  long w = bitmap[(index+offset)/64];
  if(1 & (w >> ((index+offset) % 64))) return min+1;
  else return min;
} 

void timeWT(symbol* s, long n, int rounds,char* inFile, char* outFile,
            int check) {
  //sigma
  long k = 1 + sequence::reduce(s, n, utils::maxF<symbol>());
  int* A = newA(int,k+1);
  parallel_for(long i=0;i<k+1;i++) A[i] = 0;
  parallel_for(long i=0;i<n;i++) {
    if(!A[s[i]]) A[s[i]] = 1;
  }
  long sigma = sequence::plusScan(A,A,k+1);
  //cout << "n = " << n << endl;
  //cout << "sigma = " << sigma << endl;
  int* reverseMap = newA(int,sigma);
  parallel_for(long i=0;i<k;i++) 
    if(A[i] != A[i+1]) reverseMap[A[i]] = i;

  parallel_for(long i=0;i<n;i++) s[i] = A[s[i]];
  free(A);

  std::cout << "RESULT algo=wt_sdsl_pc ";

  // Version for sdsl
  typedef sdsl::int_vector<sizeof(symbol)*8> input_type;
  input_type input(n, 0);
  for (int64_t i = 0; i < n; i++) {
	input[i] = s[i];
  }
#ifdef INT
  typedef sdsl::wt_int<sdsl::bit_vector,
      sdsl::rank_support_scan<>,
      sdsl::select_support_scan<>,
      sdsl::select_support_scan<0>> wt_type;
#else
  typedef sdsl::wt_pc<sdsl::balanced_shape,
      sdsl::bit_vector,
      sdsl::rank_support_scan<>,
      sdsl::select_support_scan<>,
      sdsl::select_support_scan<0>> wt_type;
#endif
  string input_file = "@input.sdsl"; 
  sdsl::store_to_file(input, input_file);
  input = input_type();
  // Timing is done directly in construct method.
  wt_type wt;
  sdsl::construct(wt, input_file);
  wt = wt_type();
  std::vector<float> times;
  for (int i=0; i < rounds; i++) {
    auto begin_time = std::chrono::high_resolution_clock::now();
    sdsl::construct(wt, input_file);
    wt = wt_type();
    auto end_time = std::chrono::high_resolution_clock::now();
    times.emplace_back(static_cast<float>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - begin_time).count()));
  }

  std::cout << "memory=" << getPeakRSS() << ' '
            << "runs=" << rounds << ' ';

  std::sort(times.begin(), times.end());
  std::cout << "median_time=";
  if (rounds % 2 == 0) {
    std::cout << (times[rounds >> 1] + times[(rounds >> 1) - 1]) / 2;
  } else {
    std::cout << times[rounds >> 1]; 
  }
  std::cout << " input=" << inFile << ' '
            << "characters=" << n << ' '
            << "sigma=" << sigma << ' '
            << "word_width=" << sizeof(symbol) <<' '
            << "threads=1" << std::endl;

  if(check) {
	  std::cout << "Checking is not implemented";
	  //TODO 
  }
  if (outFile != NULL) {
	sdsl::store_to_file(wt, outFile);
  }
}


int parallel_main(int argc, char* argv[]) {
  std::string iFile;
  std::string oFile;
  int rounds = 5;
  bool binary = false;
  int check = 0;
  uint64_t prefix_size = 0;

  tlx::CmdlineParser cp;
  cp.add_param_string("input", iFile, "Path to the input text");
  cp.add_string('o', "out_file", oFile,
                "Optional file where the output is stored");
  cp.add_int('r', "rounds", rounds, "Number of executions of the algorithm");
  cp.add_flag('b', "binary", binary, "Accept binary input");
  cp.add_int('c', "check", check,
             "Number of queries used to verify correctness of the computation");
  cp.add_bytes('l', "length", prefix_size,
                "Length of the prefix of the text that should be considered");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  if(binary) {
    ifstream in(iFile,ifstream::in |ios::binary);
    in.seekg(0,ios::end);
    long n = in.tellg();
    in.seekg(0);
    char* s = newA(char, n);
    in.read(s,n);
    in.close(); 
    timeWT((symbol*)s, n/sizeof(symbol), rounds, iFile.data(), oFile.data(),
           check);
    free(s);
  }
  else {
#ifdef INT
    _seq<uintT> S = readIntArrayFromFile<uintT>(iFile.data());
    uintT n = S.n;
    timeWT(S.A, n, rounds, fFile.data(), oFile.data(), check);
#else
    _seq<char> S = readStringFromFile(iFile.data());
    uintT n = S.n;
    timeWT((unsigned char*) S.A, n, rounds, iFile.data(), oFile.data(), check);
#endif
    S.del();
  }
}

