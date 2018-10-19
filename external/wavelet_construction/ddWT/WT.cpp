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

struct chunkInfo {
  uintT start, length, nodeID, range;
};
// nodes: Array of atleast [2*sigma] nodes
// wt: array of atleast [log(sigma)*n] bits
// Remarks: nodes, and wt have to be allocated by the caller and are filled by the function 
void serialWT(symbol* s, uintT n, uintT sigma, WTnode* nodes, long* wt) {
  symbol* s1 = newA(symbol,n);
  symbol* s2 = newA(symbol,n);
#ifdef POINTERS
  nodes[0].parent = UINT_T_MAX;
#endif
  chunkInfo* offsets = newA(chunkInfo,sigma),
    *offsets2 = newA(chunkInfo,sigma);
  offsets[0].start=0; offsets[0].length=n; offsets[0].nodeID=0; offsets[0].range=1<<utils::log2Up(sigma);
  int levels = max(1,utils::log2Up(sigma));
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
	for(uintT i=start;i<start+length;i++)
	  if(!(source[i] & mask)) right++;
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
}


void copy_wt(long* wt_in, uintT s_in, long* wt_out, uintT s_out, uintT len) {
	for (uintT i = 0; i < len; ++i) {
		if (wt_in[(s_in + i) / 64] & (long)1<<((s_in + i) % 64))
			wt_out[(s_out + i) / 64] |= (long)1<<((s_out + i) % 64);
	}
}


		

pair<WTnode*,long*> WT(symbol* input_str, uintT n, uintT sigma) {
  int num_threads = getWorkers();
  // Allocate needed space
  uintT num_total_nodes = 2*sigma*num_threads;
  WTnode* prenodes = newA(WTnode, num_total_nodes);
  int levels = max(1,utils::log2Up(sigma));
  long* prewt = newA(long,((long)n*levels+63)/64);

  // On default nodes are empty:
  WTnode* nodes = newA(WTnode, 2*sigma);
  parallel_for (uintT i = 0; i < num_total_nodes; ++i)
	  prenodes[i].length = 0;
  // Default wt bit is 0 
  long* wt = newA(long, ((long)n*levels+63)/64);
  // Bock_size is n/p rounded up and then rounded to next multiple of 64
  uintT block_size = ((n+num_threads-1)  / num_threads + 63) / 64 * 64;
  parallel_for(long i=0;i<(n*levels+63)/64; i++) wt[i] = 0;
  // Build wt per thread
  {
	  // Every thread builds one wt
	  parallel_for (int p = 0; p < num_threads;p++) {
		// Important that no word overlaps between the wt arrays
		uintT s = p * block_size;		
		uintT e = std::min((p+1)*block_size, n);
		// for the last argument it is importtant that s is always multiple of 64
		serialWT(input_str + s, (e-s), sigma, prenodes + p*2*sigma, prewt + (s*levels/64));	
	  }
  }
  // Arrays used for reordering
  uintT* node_lengths = newA(uintT, num_total_nodes);
  uintT* node_bitmapPtrs_old = newA(intOffset, num_total_nodes);
  uintT* node_bitmapPtrs_new = newA(intOffset, num_total_nodes);
  parallel_for (uintT i = 0; i < num_total_nodes; ++i) {
	  uintT p = i % num_threads;
	  uintT node_id = i / num_threads;
	  uintT i_old = p*2*sigma+node_id;
	  node_lengths[i] = prenodes[i_old].length;
	  node_bitmapPtrs_new[i] = prenodes[i_old].length;
	  node_bitmapPtrs_old[i] = prenodes[i_old].bitmapPtr + p*block_size*levels;
  }
  
  // Prefix sum on node_lengths to get new offset into wt array
  uintT result_size = sequence::prefixSum<uintT>(node_bitmapPtrs_new, 0, num_total_nodes); 
  // Set new nodes information
  parallel_for (uintT node_id = 0; node_id < 2*sigma; ++node_id) {
	uintT s = node_id*num_threads;
	uintT e = s + num_threads -1;
	nodes[node_id].length = node_bitmapPtrs_new[e] - node_bitmapPtrs_new[s] + node_lengths[e];
	nodes[node_id].bitmapPtr = node_bitmapPtrs_new[s];
 } 
 // Add pointers
#ifdef POINTERS
 nodes[0].parent = UINT_T_MAX;
 parallel_for (uintT node_id = 0; node_id < 2*sigma; ++node_id) {
	nodes[node_id].leftChild = nodes[node_id].rightChild = UINT_T_MAX;
	if (2*node_id+2 < 2*sigma) {
		// Left child
		if (nodes[2*node_id+1].length > 0) {
			nodes[2*node_id+1].parent = node_id;
			nodes[node_id].leftChild = 2*node_id+1;
		}
		// Right child
		if (nodes[2*node_id+2].length > 0) {
			nodes[2*node_id+2].parent = node_id;
			nodes[node_id].rightChild = 2*node_id+2;
		}
	}
 }
#endif

  // Actually copy the wt data in large blocks
  {
	  uintT block_size = result_size / (num_threads*4) / 64 *64;
	  uintT num_blocks = (result_size + block_size - 1) / block_size;
	  parallel_for (uintT b = 0; b < num_blocks; ++b) {
		uintT s = block_size * b;	
		uintT e = std::min(block_size * (b+1), result_size);
		// Binary search for start of the block
		uintT* array_end = node_bitmapPtrs_new + num_total_nodes;
		uintT* cur_pos_ptr = std::upper_bound(node_bitmapPtrs_new, array_end, s);
		uintT cur_pos = cur_pos_ptr - node_bitmapPtrs_new - 1; // Position of input node in node array
		while (cur_pos < num_total_nodes) {
			if (node_bitmapPtrs_new[cur_pos] >= e) // Done with the block
				break;		
			uintT s_new = std::max(s, node_bitmapPtrs_new[cur_pos]); 	
			uintT e_new = std::min(e, node_bitmapPtrs_new[cur_pos] + node_lengths[cur_pos]);
			uintT s_old = node_bitmapPtrs_old[cur_pos];
			// If node was splitted, move copy position same offset forward
			if (s_new > node_bitmapPtrs_new[cur_pos])
				s_old += s_new - node_bitmapPtrs_new[cur_pos];
			copy_wt(prewt, s_old, wt, s_new, e_new - s_new);
			cur_pos++;
		}
	  }
  }
  free(node_lengths); 
  free(node_bitmapPtrs_old);
  free(node_bitmapPtrs_new);
  free(prenodes);
  free(prewt);
  return pair<WTnode*,long*>(nodes,wt);
}



