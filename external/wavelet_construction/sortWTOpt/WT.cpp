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
#include "../parallel.hpp"
#include "../utils.hpp"
#include <iostream>
#include "../sequence.hpp"
#include "WT.hpp"
#include "blockRadixSort.hpp"
using namespace std;

#define THRESHOLD 10000

template <class E>
struct topL {
  uint shift;
  topL(int _s) { shift = _s; }
  E operator() (E a) {return a >> shift;} ; 
};

inline void writeOr(long *a, long b) {
  volatile long newV, oldV; 
  do {oldV = *a; newV = oldV | b;}
  while ((oldV != newV) && !utils::CAS(a, oldV, newV));
}

struct notMax { bool operator() (uintT i) {return i != UINT_T_MAX;}};

pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma) {
  int levels = max(1,utils::log2Up(sigma));
  //cout << "levels = " << levels << endl;

  uintT numNodes = (long)1 << levels;
  WTnode* nodes = newA(WTnode,numNodes);
  nodes[0].length = n; nodes[0].bitmapPtr = 0;
#ifdef POINTERS
  nodes[0].parent = UINT_T_MAX;
  parallel_for(uint64_t i=1;i<numNodes;i++) nodes[i].parent = UINT_T_MAX-1; //indicates not in tree yet
#endif
  long* wt = newA(long,((long)n*levels+63)/64);
  parallel_for(long i=0;i<((long)n*levels+63)/64; i++) wt[i] = 0;
  //cout << "wt length = " << (n*levels+63)/64 << endl;;
  //cout << "sigma = " << sigma << endl;

  uintT* space = newA(uintT,(long)n*(levels-1)); //temp space
  uintT* space2 = newA(uintT,(long)n*(levels-1)); //temp space

  symbol* s2 = newA(symbol,(long)n*(levels-1)); //space for reordered character string

  parallel_for_1(long i=0;i<levels-1;i++) {
    symbol* start = s2+i*n;
    parallel_for(uint64_t j=0;j<n;j++) start[j] = s[j];
  }

  parallel_for_1(int l = 0; l < levels; l++) {
    //cout << "l = " << l << endl;
    int mask = (long)1 << (levels - l - 1);

    if(l == 0) {
      uintT e = n-1, word = e/64;
      while(e % 64 != 63) {
	if(s[e] & mask) writeOr(&wt[word],(long)1 << e % 64);	
	e--;
      }
      parallel_for(uintT k=0;k<e/64;k++) {
	uintT b = 64*k;
	for(uintT i=b;i<b+64;i++) 
	  if(s[i] & mask) wt[k] |= (long)1 << i % 64;
      }
    } else {
      intOffset levelOffset = (l-1)*n;
      uintT* mySpace = space+levelOffset, *mySpace2 = space2+levelOffset;
      symbol* mys2 = s2+levelOffset;
      intOffset wtOffset = l*n;

      uintT numNodesLevel = (long)1 << l;
      uint shift = levels-l;
      
      //stably sort on top l bits to put characters in correct position
      intSort::iSort(mys2, n, numNodesLevel, topL<symbol>(shift));

      //create bitmaps for this level
      intOffset s = wtOffset, e = wtOffset+n-1, word = s/64;
      while(s % 64) {
	if(mys2[s-wtOffset] & mask) writeOr(&wt[word],(long)1 << s % 64);	
	s++;
      } 
      if(word != e/64) {
	word = e/64;
	while(e % 64 != 63) {
	  if(mys2[e-wtOffset] & mask) writeOr(&wt[word],(long)1 << e % 64);	
	  e--;
	} 
      }
      intOffset startWord = s/64, endWord = e/64;
      parallel_for(intOffset k=startWord;k<=endWord;k++) {
	intOffset b = 64*k;
	for(intOffset i=b;i<b+64;i++) 
	  if(mys2[i-wtOffset] & mask) wt[k] |= (long)1 << i % 64;
      }

      //determine offsets and compute number of elements per node
      mySpace[0] = 0;
      parallel_for(uint64_t i=1;i<n;i++) {
	mySpace[i] = ((mys2[i] >> shift) != (mys2[i-1] >> shift)) ? i : UINT_T_MAX; 
      }

      uintT actualNumNodes = sequence::filter(mySpace,mySpace2,n,notMax());
      //cout << "numNodesLevel = " << numNodesLevel << " actualNumNodes = " << actualNumNodes << endl;
      if(actualNumNodes < THRESHOLD) {
	for(uint64_t i=0;i<actualNumNodes;i++) {
	  uintT nodeID = numNodesLevel-1+(mys2[mySpace2[i]] >> shift);
	  uintT length = (i == actualNumNodes-1) ? (n-mySpace2[i]) : (mySpace2[i+1]-mySpace2[i]);
	  nodes[nodeID].length = length;
	  nodes[nodeID].bitmapPtr = wtOffset+mySpace2[i];
#ifdef POINTERS
	  nodes[nodeID].parent = (nodeID-1)/2; //put node in tree
#endif
	}
      } 
      else {
	parallel_for(uint64_t i=0;i<actualNumNodes;i++) {
	  uintT nodeID = numNodesLevel-1+(mys2[mySpace2[i]] >> shift);
	  uintT length = (i == actualNumNodes-1) ? (n-mySpace2[i]) : (mySpace2[i+1]-mySpace2[i]);
	  nodes[nodeID].length = length;
	  nodes[nodeID].bitmapPtr = wtOffset+mySpace2[i];
#ifdef POINTERS
	  nodes[nodeID].parent = (nodeID-1)/2; //put node in tree
#endif
	}
      }
    }
    //cout << "finishing level " << l << endl;
  }

#ifdef POINTERS
  //fix child pointers of tree
  parallel_for_1(uint64_t i=0;i<numNodes/2-1;i++) {
    if(nodes[2*i+1].parent == i) nodes[i].leftChild = 2*i+1; 
    else nodes[i].leftChild = UINT_T_MAX;
    if(nodes[2*i+2].parent == i) nodes[i].rightChild = 2*i+2; 
    else nodes[i].rightChild = UINT_T_MAX;
  }
#endif

  free(space); free(space2); 
  free(s2); 
  return make_pair(nodes,wt);
}
