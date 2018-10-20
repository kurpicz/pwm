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
  unsigned long start, length, nodeID, range;
};

struct isChild { bool operator() (chunkInfo i) {return i.start != UINT_T_MAX;}};



void recursiveWT(symbol* s, long n, uintT sigma, chunkInfo& current_chunk, symbol* source, symbol* destination, WTnode* nodes, long* wt, int levels, int l) { 
	if (l == 0)
		std::swap(s,source);
		
	intOffset levelOffset = l*n;
	int mask = (long)1 << (levels - l - 1);
	long start = current_chunk.start, length = current_chunk.length, nodeID = current_chunk.nodeID, range = current_chunk.range;
	intOffset myOffset = (intOffset)start+levelOffset;
	nodes[nodeID].bitmapPtr = myOffset; nodes[nodeID].length = length;
	chunkInfo left_child,right_child;
      
	intOffset sourceOffset = start-myOffset;
      intOffset wtBegin = myOffset, wtEnd = myOffset+length-1;
      intOffset word = wtBegin/64;
      cout << "range: " << range << endl;
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
	cout << "set left_child.start: " << left_child.start << " vs. " << UINT_T_MAX << std::endl;
	cout << "set right_child.start: " << right_child.start << " vs. " << UINT_T_MAX << std::endl;
#endif
	left_child.start = right_child.start = UINT_T_MAX;
      } else {
	intOffset wtBegin_ = wtBegin;
	intOffset wtEnd_ = wtEnd;	
	if(length < 128) {
	  for(intOffset i=wtBegin;i<=wtEnd;i++) {
	    if(source[sourceOffset+i] & mask)
	      writeOr(&wt[i/64],(long)1 << (i % 64));
	  }
	}
	else if(length < THRESHOLD) {
	  while(wtBegin%64) { //shares first word with someone else
	    if(source[sourceOffset+wtBegin] & mask)  
		    writeOr(&wt[word],(long)1 << (wtBegin % 64));
	    wtBegin++;
	  }
	  word = wtEnd/64;
	  while(wtEnd%64 != 63) { //shares last word with someone else
	    if(source[sourceOffset+wtEnd] & mask) 
		    writeOr(&wt[word],(long)1 << (wtEnd % 64));
	      //position on right is equal to index minus position on left
	    wtEnd--;
	  }
	  for(intOffset i=wtBegin;i<=wtEnd;i++) {
	    if(source[sourceOffset+i] & mask) 
	      wt[i/64] |= (long)1 << (i % 64);
	  }
	} else {
	  while(wtBegin%64) { //shares first word with someone else
	    if(source[sourceOffset+wtBegin] & mask) 
		    writeOr(&wt[word],(long)1 << (wtBegin % 64));
	      //position on right is equal to index minus position on left
	    wtBegin++; 
	  }
	  word = wtEnd/64;
	  while(wtEnd%64 != 63) { //shares last word with someone else
	    if(source[sourceOffset+wtEnd] & mask) 
		    writeOr(&wt[word],(long)1 << (wtEnd % 64));
	      //position on right is equal to index minus position on left
	    wtEnd--;
	  }
	  //loop over rest in chunks
	  intOffset startWord = wtBegin/64, endWord = wtEnd/64;
	  parallel_for(intOffset k=startWord;k<=endWord;k++) {
	    intOffset b = 64*k;
	    for(intOffset i=b;i<b+64;i++) {
	      if(source[sourceOffset+i] & mask) 
		      wt[k] |= (long)1 << (i % 64);
	    }
	  }
	}
	// Reorder s1 s2 correctly
	// Check that input and flags are indexed the same, but output is indexed starting from 0
	intOffset rightStart = sequence::pack2Bit(source-levelOffset , destination+start , wt, wtBegin_, wtEnd_+1); 
	if(rightStart) {
#ifdef POINTERS
	  nodes[2*nodeID+1].parent = nodeID;
	  nodes[nodeID].leftChild = 2*nodeID+1;
#endif
	  left_child.start = start; left_child.length = rightStart; left_child.nodeID = 2*nodeID+1; left_child.range = (long)1 << (utils::log2Up(range)-1);	  
	} else { left_child.start = UINT_T_MAX; 
#ifdef POINTERS
	  nodes[nodeID].leftChild = UINT_T_MAX; 
#endif
	}
	if(length-rightStart) {
#ifdef POINTERS
	  nodes[2*nodeID+2].parent = nodeID;
	  nodes[nodeID].rightChild = 2*nodeID+2;
#endif
	  right_child.start = start+rightStart; right_child.length = length-rightStart; right_child.nodeID = 2*nodeID+2; right_child.range = range - ((long)1 << (utils::log2Up(range)-1));
	} else { right_child.start = UINT_T_MAX; 
#ifdef POINTERS
	  nodes[nodeID].rightChild = UINT_T_MAX;
#endif
	}
      }
      // Spawn children
  l++;
  if (l < levels) {
	  if (l == 1)
		  swap(source,s);

		cout << "left_child.start: " << left_child.start << " vs. " << UINT_T_MAX << std::endl;
		cout << "right_child.start: " << right_child.start << " vs. " << UINT_T_MAX << std::endl;
	  if (left_child.start != UINT_T_MAX)
		  cilk_spawn recursiveWT(s, n, sigma, left_child, destination, source, nodes, wt, levels, l);
	  if (right_child.start != UINT_T_MAX) 
		  recursiveWT(s, n, sigma, right_child, destination, source, nodes, wt, levels, l);
  }
}


pair<WTnode*,long*> WT(symbol* s, long n, uintT sigma) {
  symbol* s1 = newA(symbol,n);
  symbol* s2 = newA(symbol,n);
  WTnode* nodes = newA(WTnode,sigma*2);
#ifdef POINTERS
  nodes[0].parent = UINT_T_MAX;
#endif
  chunkInfo start_chunk;
  start_chunk.start=0; start_chunk.length=n; start_chunk.nodeID=0; start_chunk.range=((long)1 << utils::log2Up(sigma));
  int levels = max(1,utils::log2Up(sigma));
  long* wt = newA(long,((long)n*levels+63)/64);
  parallel_for(long i=0;i<((long)n*levels+63)/64; i++) wt[i] = 0;
  recursiveWT(s, n, sigma, start_chunk, s1, s2, nodes, wt, levels, 0);
  free(s2); 
  free(s1);
  return make_pair(nodes,wt);
}
