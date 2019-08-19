#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "src/bytearray.hpp"
#include "src/cstringutil.hpp"
#include "src/strstream.hpp"
#include "src/wavtree.hpp"

#include <tlx/cmdline_parser.hpp>

void online(std::string const& iFile, int const rounds,
            uint64_t const prefix_size) {
  char* filename;
  std::copy(iFile.begin(), iFile.end(), filename);
  clock_t start, finish;
  strstream *sst; 
  sst = strstream_open_file(filename);
  wavtree *wt;
  std::cout << "RESULT algo=wt_sort ";
  start = clock();
  wt = wavtree_new_online_from_stream(sst, WT_BALANCED);
  finish = clock();
  printf("wavelet_tree created\n");
  fprintf(stderr, "%.4f sec\n", (double)(finish - start) / (double)CLOCKS_PER_SEC);
  //strstream_close(sst);
  //wavtree_free(wt);
}

int main(int argc, char *argv[])
{
  std::string iFile;
  int rounds = 5;
  uint64_t prefix_size = 0;

  tlx::CmdlineParser cp;
  cp.add_param_string("input", iFile, "Path to the input text");
  cp.add_int('r', "rounds", rounds, "Number of executions of the algorithm");
  cp.add_bytes('l', "length", prefix_size,
                "Length of the prefix of the text that should be considered");

  if (!cp.process(argc, argv)) {
    return -1;
  }

  online(iFile, rounds, prefix_size);
  return 0;
}
