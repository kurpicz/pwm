/******************************************************************************
 * pwt.c
 *
 * Parallel per-level algorithm to the construction of Wavelet trees
 *
 ******************************************************************************
 * Copyright (C) 2015 José Fuentes-Sepúlveda <jfuentess@udec.cl>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <sys/time.h>
#include <time.h>

#ifdef NO_CILK
#define __cilkrts_get_nworkers() 1
#define cilk_for for
#define cilk_spawn 
#define cilk_sync 
#else
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/common.h>
#endif

#include "basic_wt.h"
#include "util.h"

#ifdef MALLOC_COUNT
#include "benchmark/malloc_count.h"
#endif

#define FACTOR 20 /* bitmap usage*/

/* Create the wavelet tree*/
BIT_ARRAY** wt_create(symbol* T, unsigned long n, unsigned int alphabet) {

  unsigned int levels = log2(alphabet);

  BIT_ARRAY** wtree = (BIT_ARRAY**)malloc(levels * sizeof(BIT_ARRAY*));
  unsigned int i = 0;
  cilk_for (i = 0; i < levels; i++) {
    /* Number of nodes in the i-th level*/
    unsigned int nnum = 1 << i;

    /* Allocating memory to save the number of symbols on each node of the wavelet tree*/
    unsigned int* counters = (unsigned int*)calloc(nnum, sizeof(unsigned int));
    wtree[i] = bit_array_create(n);

    symbol k;
    unsigned int schunk = 0;
    unsigned int desp1 = levels - i, desp2 = desp1 - 1;

    /* Counting the number of symbols on each node*/
    unsigned long j = 0;
    for (j = 0; j < n; j++) {
      k = T[j];
      schunk = (unsigned int)k >> desp1;
      counters[schunk]++;
    }

    /* Offsets*/
    unsigned int cnt = 0, aux = 0;
    unsigned int z = 0;
    for(z = 0; z < nnum; z++) {
      aux = counters[z];
      counters[z] = cnt;
      cnt += aux;
    }

    /* Filling the bit arrays*/
    for (j = 0; j < n; j++) {
      k = T[j];
      schunk = (unsigned int)k >> desp1;
      k &= 1 << desp2;
      
      if (k)
        bit_array_set_bit(wtree[i], counters[schunk]);
      counters[schunk]++;
    }

#ifdef ARCH64
    fprintf(stderr, "createBitRankW32Int() needs to be adapted to work on 64-bits mode");
    exit(1);
#endif

    /* wtree[i] = createBitRankW32Int(bit_array->words, bit_array->num_of_bits, 1, FACTOR); */

    /* free(bit_array); */
    free(counters);
  }

  return wtree;
}

unsigned int _wt_rank_0(bitRankW32Int* level, unsigned int idx) {
  return ((1 + idx) - rank(level, idx));
}
unsigned int _wt_rank_1(bitRankW32Int* level, unsigned int idx) {
  return rank(level, idx);
}

unsigned int wt_access(bitRankW32Int** wtree, unsigned int idx, unsigned int alphabet) {

  unsigned int lls = 0; /* Lower limit of the sequence*/
  unsigned int uls = lenght_in_bits(wtree[0]) - 1; /* Upper limit of the sequence*/

  unsigned int levels = log2(alphabet);
  unsigned int lla = 0; /* Lower limit of the alphabet*/
  unsigned int ula = alphabet - 1; /* Upper limit of the alphabet*/

  unsigned int crank = 0;
  unsigned int partial_crank = 0;

  unsigned int i = 0;
  for (i = 0; i < levels; i++) {
    if(isBitSet(wtree[i], idx + lls) != 1) {
      crank = _wt_rank_0(wtree[i], lls + idx) - _wt_rank_0(wtree[i], lls - 1);
      idx = crank - 1;
      ula = (lla + ula) / 2;
      uls -= _wt_rank_1(wtree[i], uls) - _wt_rank_1(wtree[i], lls - 1);

    } else {
      crank = _wt_rank_1(wtree[i], lls + idx) - _wt_rank_1(wtree[i], lls - 1); /* rank 1*/
      idx = crank - 1;
      lla = ((lla + ula) / 2) + 1;
      lls += _wt_rank_0(wtree[i], uls) - _wt_rank_0(wtree[i], lls - 1);
    }
  }

  return lla;
}

/* Print the wavelet tree*/
void wtree_to_string(bitRankW32Int** wtree, unsigned int alphabet) {
  /* Number of levels of the wavelet tree*/
  unsigned int levels = log2(alphabet);

  unsigned int i = 0;
  for (i = 0; i < levels; i++) {
    BIT_ARRAY p;
    p.words = wtree[i]->data;
    p.num_of_bits = wtree[i]->n;

    char* st = bit_array_to_string(&p);
    printf("%s\n", st);
    free(st);
  }
}

int main(int argc, char* argv[]) {

  if(argc != 3 && argc != 4){
    fprintf(stderr, "Execute: %s <input file> <alphabet size> [<validation_file>]\n", argv[0]);
    exit(-1);
  }

  unsigned long n; /* Size of the input sequence*/
  symbol* text = read_text_from_file(argv[1], &n); /* Input sequence*/
  unsigned int alphabet = (unsigned int)atoi(argv[2]); /* size of the alphabet*/

  /* Memory usage*/
#ifdef MALLOC_COUNT
  /* size_t s_total_memory = malloc_count_total(); */
  /* size_t s_current_memory = malloc_count_current(); */
  /* malloc_reset_peak(); */

  /* Running time. CLOCK_THREAD_CPUTIME_ID: Running time of the thread that call it (main thread in this case)*/
#else
  struct timespec stime, etime;
  double t;
  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID , &stime)) {
    fprintf(stderr, "clock_gettime failed");
    exit(-1);
  }
#endif

  /* Wavelet tree construction*/
  BIT_ARRAY** wtree = wt_create(text, n, alphabet);

#ifdef MALLOC_COUNT
  /* size_t e_total_memory = malloc_count_total(); */
  /* size_t e_current_memory = malloc_count_current(); */
  /* printf("%s, %u, %zu, %zu, %zu, %zu, %zu\n", argv[1], alphabet, s_total_memory, e_total_memory, malloc_count_peak(), s_current_memory, e_current_memory); */

#else
  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID , &etime)) {
    fprintf(stderr, "clock_gettime failed");
    exit(-1);
  }

  t = (etime.tv_sec - stime.tv_sec) + (etime.tv_nsec - stime.tv_nsec) / 1000000000.0;
  printf("%d,%s,%lu,%lf\n", __cilkrts_get_nworkers(), argv[1], n, t);
#endif

  free(text);

  return EXIT_SUCCESS;
}
