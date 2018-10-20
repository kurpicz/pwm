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
using namespace std;

#define THRESHOLD 10000

inline void writeOr(long *a, long b) {
  volatile long newV, oldV; 
  do {oldV = *a; newV = oldV | b;}
  while ((oldV != newV) && !utils::CAS(a, oldV, newV));
}

struct chunkInfo {
  uintT start, length, nodeID, range;
};

struct isChild { bool operator() (chunkInfo i) {return i.start != UINT_T_MAX;}};

pair<WTnode*,long*> WT(symbol* s, uintT n, uintT sigma) {
  symbol* s1 = newA(symbol,n);
  symbol* s2 = newA(symbol,n);
  WTnode* nodes = newA(WTnode,sigma*2);
#ifdef POINTERS
  nodes[0].parent = UINT_T_MAX;
#endif
  chunkInfo* offsets = newA(chunkInfo,sigma),
    *offsets2 = newA(chunkInfo,sigma);
  offsets[0].start=0; offsets[0].length=n; offsets[0].nodeID=0; offsets[0].range=(1 << utils::log2Up(sigma));
  int levels = max(1,utils::log2Up(sigma));
  long* wt = newA(long,((long)n*levels+63)/64);
  parallel_for(long i=0;i<((long)n*levels+63)/64; i++) wt[i] = 0;
  uintT numChunks = 1;

  uintT* space = newA(uintT,n); //temp space

  for(int l = 0; l < levels; l++) {
    intOffset levelOffset = l*n;
    int mask = (long)1 << (levels - l - 1);
    symbol* source = (l == 0) ? s : s1;
    parallel_for(uintT j=0;j<numChunks;j++) {
      uintT start = offsets[j].start, length = offsets[j].length, nodeID = offsets[j].nodeID, range = offsets[j].range;
      intOffset myOffset = (intOffset)start+levelOffset;
      nodes[nodeID].bitmapPtr = myOffset; nodes[nodeID].length = length;
      
      intOffset sourceOffset = start-myOffset;
      intOffset wtBegin = myOffset, wtEnd = myOffset+length-1;
      intOffset word = wtBegin/64;
      if(range <= 2) {
	if(length < 128) {
	  for(intOffset i=wtBegin;i<=wtEnd;i++) {
	    if(source[sourceOffset+i] & mask)
	      writeOr(&wt[i/64],(long)1 << (i % 64));
	  }
	}
	else if(length < THRESHOLD) {
	  while(wtBegin%64) { //shares first word with someone else
	    if(source[sourceOffset+wtBegin] & mask) writeOr(&wt[word],(long)1 << (wtBegin % 64));
	    wtBegin++; 
	  }
	  if(word != wtEnd/64) { //check if first word is same as last word
	    word = wtEnd/64;
	    while(wtEnd%64 != 63) { //shares last word with someone else
	      if(source[sourceOffset+wtEnd] & mask) writeOr(&wt[word],(long)1 << (wtEnd % 64));
	      wtEnd--;
	    }
	  }
	  for(intOffset i=wtBegin;i<=wtEnd;i++) {
	    if(source[sourceOffset+i] & mask)
	      wt[i/64] |= (long)1 << (i % 64);
	  }
	} 
	else {
	  while(wtBegin%64) { //shares first word with someone else
	    if(source[sourceOffset+wtBegin] & mask) writeOr(&wt[word],(long)1 << (wtBegin % 64));
	    wtBegin++; 
	  }
	  word = wtEnd/64;
	  while(wtEnd%64 != 63) { //shares last word with someone else
	    if(source[sourceOffset+wtEnd] & mask) writeOr(&wt[word],(long)1 << (wtEnd % 64));
	    wtEnd--;
	  }
	  intOffset startWord = wtBegin/64, endWord = wtEnd/64;
	  parallel_for(intOffset k=startWord;k<=endWord;k++) {
	    intOffset b = 64*k;
	    for(intOffset i=b;i<b+64;i++) 
	      if(source[sourceOffset+i] & mask) wt[k] |= (long)1 << (i % 64);
	  }
	}
#ifdef POINTERS
	nodes[nodeID].leftChild = nodes[nodeID].rightChild = UINT_T_MAX;
#endif
	offsets2[j*2].start = offsets2[j*2+1].start = UINT_T_MAX;
      } else {
	//count up number of elts on the left
	if(length < THRESHOLD) {
	  for(uintT i=start;i<start+length;i++) {
	    space[i] = (!(source[i] & mask)) ? 1 : 0;
	  }
	}
	else {
	  parallel_for(uintT i=start;i<start+length;i++) {
	    space[i] = (!(source[i] & mask)) ? 1 : 0;
	  }
	}
	// perform prefix sum on one bits
	intOffset rightStart = sequence::plusScan(space+start,space+start,length);

	if(length < 128) {
	  for(intOffset i=wtBegin;i<=wtEnd;i++) {
	    if(source[sourceOffset+i] & mask) {
	      writeOr(&wt[i/64],(long)1 << (i % 64));
	      s2[rightStart+sourceOffset+i-space[sourceOffset+i]] = source[sourceOffset+i];
	    }
	    else s2[start+space[sourceOffset+i]] = source[sourceOffset+i];
	  }
	}
	else if(length < THRESHOLD) {
	  while(wtBegin%64) { //shares first word with someone else
	    if(source[sourceOffset+wtBegin] & mask) { writeOr(&wt[word],(long)1 << (wtBegin % 64));
	      //position on right is equal to index minus position on left
	      s2[rightStart+sourceOffset+wtBegin-space[sourceOffset+wtBegin]] = source[sourceOffset+wtBegin]; }
	    else {s2[start+space[sourceOffset+wtBegin]] = source[sourceOffset+wtBegin]; }
	    wtBegin++;
	  }
	  word = wtEnd/64;
	  while(wtEnd%64 != 63) { //shares last word with someone else
	    if(source[sourceOffset+wtEnd] & mask) { writeOr(&wt[word],(long)1 << (wtEnd % 64));
	      //position on right is equal to index minus position on left
	      s2[rightStart+sourceOffset+wtEnd-space[sourceOffset+wtEnd]] = source[sourceOffset+wtEnd]; }
	    else {s2[start+space[sourceOffset+wtEnd]] = source[sourceOffset+wtEnd]; }
	    wtEnd--;
	  }
	  for(intOffset i=wtBegin;i<=wtEnd;i++) {
	    if(source[sourceOffset+i] & mask) {
	      wt[i/64] |= (long)1 << (i % 64);
	      s2[rightStart+sourceOffset+i-space[sourceOffset+i]] = source[sourceOffset+i];
	    }
	    else {s2[start+space[sourceOffset+i]] = source[sourceOffset+i]; }
	  }
	} else {
	  while(wtBegin%64) { //shares first word with someone else
	    if(source[sourceOffset+wtBegin] & mask) { writeOr(&wt[word],(long)1 << (wtBegin % 64));
	      //position on right is equal to index minus position on left
	      s2[rightStart+sourceOffset+wtBegin-space[sourceOffset+wtBegin]] = source[sourceOffset+wtBegin]; }
	    else {s2[start+space[sourceOffset+wtBegin]] = source[sourceOffset+wtBegin]; }
	    wtBegin++; 
	  }
	  word = wtEnd/64;
	  while(wtEnd%64 != 63) { //shares last word with someone else
	    if(source[sourceOffset+wtEnd] & mask) { writeOr(&wt[word],(long)1 << (wtEnd % 64));
	      //position on right is equal to index minus position on left
	      s2[rightStart+sourceOffset+wtEnd-space[sourceOffset+wtEnd]] = source[sourceOffset+wtEnd]; }
	    else {s2[start+space[sourceOffset+wtEnd]] = source[sourceOffset+wtEnd]; }
	    wtEnd--;
	  }
	  //loop over rest in chunks
	  intOffset startWord = wtBegin/64, endWord = wtEnd/64;
	  parallel_for(intOffset k=startWord;k<=endWord;k++) {
	    intOffset b = 64*k;
	    for(intOffset i=b;i<b+64;i++) {
	      if(source[sourceOffset+i] & mask) { wt[k] |= (long)1 << (i % 64);
	  	s2[rightStart+sourceOffset+i-space[sourceOffset+i]] = source[sourceOffset+i];
	      } else {s2[start+space[sourceOffset+i]] = source[sourceOffset+i];}
	    }
	  }
	}
	if(rightStart) {
#ifdef POINTERS
	  nodes[2*nodeID+1].parent = nodeID;
	  nodes[nodeID].leftChild = 2*nodeID+1;
#endif
	  uintT j2 = j*2;
	  offsets2[j2].start=start; offsets2[j2].length=rightStart; offsets2[j2].nodeID=2*nodeID+1; offsets2[j2].range=(long)1 << (utils::log2Up(range)-1);
	} else { offsets2[j*2].start = UINT_T_MAX; 
#ifdef POINTERS
	  nodes[nodeID].leftChild = UINT_T_MAX; 
#endif
	}
	if(length-rightStart) {
#ifdef POINTERS
	  nodes[2*nodeID+2].parent = nodeID;
	  nodes[nodeID].rightChild = 2*nodeID+2;
#endif
	  uintT j2p1 = j*2+1;
	  offsets2[j2p1].start=start+rightStart; offsets2[j2p1].length =length-rightStart; offsets2[j2p1].nodeID = 2*nodeID+2; offsets2[j2p1].range = range - ((long)1 << (utils::log2Up(range)-1));
	} else { offsets2[j*2+1].start = UINT_T_MAX; 
#ifdef POINTERS
	  nodes[nodeID].rightChild = UINT_T_MAX;
#endif
	}
      }
    }

    numChunks = sequence::filter(offsets2,offsets,2*numChunks,isChild());
    swap(s1,s2);
    
  }
  free(space); 
  free(offsets); free(offsets2);
  free(s2); 
  free(s1);
  return make_pair(nodes,wt);
}
