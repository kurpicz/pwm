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
#include <iostream>

#include "parallel.hpp"
#include "utils.hpp"
#include "WT.hpp"

using namespace std;

struct chunkInfo {
  uintT start, length, nodeID, range;
};

pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma) {
  symbol* s1 = newA(symbol,n);
  symbol* s2 = newA(symbol,n);
  WTnode* nodes = newA(WTnode,sigma*2);
#ifdef POINTERS
  nodes[0].parent = UINT_T_MAX;
#endif
  chunkInfo* offsets = newA(chunkInfo,sigma), *offsets2 = newA(chunkInfo,sigma);
  offsets[0].start=0; offsets[0].length=n; offsets[0].nodeID=0; offsets[0].range=1<<utils::log2Up(sigma);
  int levels = max(1,utils::log2Up(sigma));
  long* wt = newA(long,((long)n*levels+63)/64);
  for(long i=0;i<(n*levels+63)/64; i++) wt[i] = 0;
  uintT numChunks = 1;
  for(int l = 0; l < levels; l++) {
    long currPos = l*n;
    int mask = (long)1 << (levels - l - 1);
    uintT nextNumChunks = 0;
    symbol* source = (l == 0) ? s : s1;
    for(long j=0;j<numChunks;j++) {
      uintT left = 0, right = 0, start = offsets[j].start, 
	          length = offsets[j].length, nodeID = offsets[j].nodeID, range = offsets[j].range;
      nodes[nodeID].bitmapPtr = currPos; nodes[nodeID].length = length;
      //cout << "("<<start<<","<<length<<","<<nodeID<<","<<range<<")\n";
      if(range <= 2) {
	      for(uintT i=start;i<start+length;i++) {
	        if(source[i] & mask) wt[currPos/64] |= (long)1 << (currPos % 64);
	        currPos++;
	      }
        #ifdef POINTERS
	        nodes[nodeID].leftChild = nodes[nodeID].rightChild = UINT_T_MAX;
        #endif
      } else {
	      for(uintT i=start;i<start+length;i++) {
	        if(!(source[i] & mask)) right++;
              }
	        for(uintT i=start;i<start+length;i++) {
  	        if(source[i] & mask) {
  	          wt[currPos/64] |= (long)1 << (currPos % 64);
  	          s2[start+right++] = source[i];
  	        } else s2[start+left++] = source[i];
  	        currPos++;
	        } 
	    if(left) {
        #ifdef POINTERS
        	nodes[2*nodeID+1].parent = nodeID;
        	nodes[nodeID].leftChild = 2*nodeID+1;
        #endif
	      offsets2[nextNumChunks].start=start; offsets2[nextNumChunks].length=left; offsets2[nextNumChunks].nodeID=2*nodeID+1; offsets2[nextNumChunks++].range=(long)1 << (utils::log2Up(range)-1);
	    } 
      #ifdef POINTERS
      	else nodes[nodeID].leftChild = UINT_T_MAX;
      #endif
	    if(right-left) {
        #ifdef POINTERS
          nodes[2*nodeID+2].parent = nodeID;
          nodes[nodeID].rightChild = 2*nodeID+2;
        #endif
	      offsets2[nextNumChunks].start=start+left; offsets2[nextNumChunks].length =right-left; offsets2[nextNumChunks].nodeID = 2*nodeID+2; offsets2[nextNumChunks++].range = range - ((long)1 << (utils::log2Up(range)-1));
	    } 
      #ifdef POINTERS
      	else nodes[nodeID].rightChild = UINT_T_MAX;
      #endif
      }
    }
    numChunks = nextNumChunks;
    swap(offsets,offsets2);
    swap(s1,s2);
  }
  free(s1); free(s2); free(offsets); free(offsets2);
  return make_pair(nodes,wt);
}
