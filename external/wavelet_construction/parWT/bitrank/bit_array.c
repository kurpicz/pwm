/*
 * commit 640451233b7958e5ae85572ea570a86a15934915
 * Author: Diego Caro <diegocaro@gmail.com>
 * Date:   Mon May 21 20:28:13 2012 -0400
 * https://github.com/diegocaro/BitArray/
 */

/*
 bit_array.c
 project: bit array C library
 url: https://github.com/noporpoise/BitArray/
 Adapted from: http://stackoverflow.com/a/2633584/431087
 author: Isaac Turner <turner.isaac@gmail.com>

 Copyright (c) 2011, Isaac Turner
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <limits.h> 
#include <string.h>
#include <math.h>

#include "bit_array.h"

#define MIN(a, b)  (((a) <= (b)) ? (a) : (b))
#define MAX(a, b)  (((a) >= (b)) ? (a) : (b))


/**/
/* For internal use*/
/**/

/* sizeof gives size in bytes (8 bits per byte)*/
int WORD_SIZE = sizeof(word_t) * 8;

/* Index of word*/
__inline__ word_addr_t bindex(bit_index_t b) {
  return b / WORD_SIZE;
}

/* Offset within a word (values up to 64 most likely)*/
__inline__ unsigned int boffset(bit_index_t b) {
  return b % WORD_SIZE;
}

/* Number of words required to store so many bits*/
__inline__ word_addr_t nwords(bit_index_t b) {
  return (b + WORD_SIZE - 1) / WORD_SIZE;
}

/*#define BIT_MASK(length)  (((word_t)1 << (length))-1) // overflows*/
#define BIT_MASK(length)  ((word_t)ULONG_MAX >> (WORD_SIZE-(length)))

/**/
/* Constructor*/
/**/
BIT_ARRAY* bit_array_create(bit_index_t nbits) {
  BIT_ARRAY* bitarr = (BIT_ARRAY*) malloc(sizeof(BIT_ARRAY));

  if(bitarr == NULL) {
    /* error - could not allocate enough memory*/
    errno = ENOMEM;
    return NULL;
  }

  word_addr_t num_of_words = nwords(nbits);

#ifdef DEBUG
  printf("Creating BIT_ARRAY (bits: %lu; words: %lu; WORD_SIZE: %i)\n",
         (unsigned long)nbits, (unsigned long)num_of_words, WORD_SIZE);
#endif

  /*bitarr->num_of_words = num_of_words;*/
  bitarr->num_of_bits = nbits;
  bitarr->words = (word_t*) calloc(sizeof(word_t), num_of_words);

  if(bitarr->words == NULL) {
    /* error - could not allocate enough memory*/
    free(bitarr);
    errno = ENOMEM;
    return NULL;
  }

  /* Initialise to zero*/
  /* bit_array_fill_zeros(bitarr);*/

  return bitarr;
}

/**/
/* Destructor*/
/**/
void bit_array_free(BIT_ARRAY* bitarr) {
  free(bitarr->words);
  free(bitarr);
}


/**/
/* Methods*/
/**/
void parallel_bit_array_set_bit(BIT_ARRAY* bitarr, bit_index_t b) {

  word_t old, new;

  if ( b >= 0 && b < bitarr->num_of_bits ) {
#ifdef ARCH64
    do {
      old = bitarr->words[b >> 6];
      new = old | ((word_t)1 << (b & 63));

    } while(!(__atomic_compare_exchange_n(&bitarr->words[b >> 6], &old, new, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST )));
#else
    do {
      old = bitarr->words[b >> 5];
      new = old | ((word_t)1 << (b & 31));

    } while(!(__atomic_compare_exchange_n(&bitarr->words[b >> 5], &old, new, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST )));
#endif
  } else {
    /* out of bounds error*/
#ifdef ARCH64
    fprintf(stderr, "bit_array.c: bit_array_set_bit() - "
            "out of bounds error (index: %lu; length: %lu)\n",
            b, bitarr->num_of_bits);
#else
    fprintf(stderr, "bit_array.c: bit_array_set_bit() - "
            "out of bounds error (index: %lu; length: %lu)\n",
            b, bitarr->num_of_bits);
#endif
    errno = EDOM;

    exit(EXIT_FAILURE);
  }
}


void bit_array_set_bit(BIT_ARRAY* bitarr, bit_index_t b) {
  if ( b >= 0 && b < bitarr->num_of_bits ) {

#ifdef ARCH64
     bitarr->words[b >> 6] |= ((word_t)1 << (b & 63));
#else
     bitarr->words[b >> 5] |= ((word_t)1 << (b & 31));
#endif
  } else {
    /* out of bounds error*/
#ifdef ARCH64
    fprintf(stderr, "bit_array.c: bit_array_set_bit() - "
            "out of bounds error (index: %lu; length: %lu)\n",
            b, bitarr->num_of_bits);
#else
    fprintf(stderr, "bit_array.c: bit_array_set_bit() - "
            "out of bounds error (index: %lu; length: %lu)\n",
            b, bitarr->num_of_bits);
#endif
    errno = EDOM;

    exit(EXIT_FAILURE);
  }
}


void bit_array_clear_bit(BIT_ARRAY* bitarr, bit_index_t b) {

  if ( b >= 0 && b < bitarr->num_of_bits ) {
    bitarr->words[bindex(b)] &= ~((word_t)1 << (boffset(b)));
  } else {
    /* out of bounds error*/
    fprintf(stderr, "bit_array.c: bit_array_clear_bit() - "
            "out of bounds error (index: %lu; length: %lu)\n",
            b, bitarr->num_of_bits);

    errno = EDOM;

    exit(EXIT_FAILURE);
  }


}

char bit_array_get_bit(BIT_ARRAY* bitarr, bit_index_t b) {

  if ( b >= 0 && b < bitarr->num_of_bits ) {
    return (bitarr->words[bindex(b)] >> (boffset(b))) & 0x1;
  } else {
    /* out of bounds error*/
    fprintf(stderr, "bit_array.c: bit_array_get_bit() - "
            "out of bounds error (index: %lu; length: %lu)\n",
            b, bitarr->num_of_bits);

    errno = EDOM;

    exit(EXIT_FAILURE);
  }
}

/* set all elements of data to zero */
void bit_array_fill_zeros(BIT_ARRAY* bitarr) {
  /*size_t num_of_bytes = (bitarr->num_of_bits / 8) + 1;*/
  size_t num_of_bytes = nwords(bitarr->num_of_bits) * sizeof(word_t);
  memset(bitarr->words, 0, num_of_bytes);
}

/* set all elements of data to one */
void bit_array_fill_ones(BIT_ARRAY* bitarr) {
  size_t num_of_bytes = nwords(bitarr->num_of_bits) * sizeof(word_t);
  memset(bitarr->words, 0xFF, num_of_bytes);
}

/* To string method (remember to free the result!)*/
char* bit_array_to_string(BIT_ARRAY* bitarr) {
  char* str = (char*) malloc(sizeof(char) * (bitarr->num_of_bits + 1));

  bit_index_t i;
  int count = 0; /* BORRAR*/

  for(i = 0; i < bitarr->num_of_bits; i++)
/*  for(i = 0; i < bitarr->num_of_bits; i+=1000)*/
  {
    str[i] = bit_array_get_bit(bitarr, i) ? '1' : '0';
/*    str[count++] = bit_array_get_bit(bitarr, i) ? '1' : '0';*/
  }

  str[bitarr->num_of_bits] = '\0';
/*  str[count] = '\0';*/

  return str;
}

BIT_ARRAY* bit_array_clone(BIT_ARRAY* bitarr) {
  BIT_ARRAY* cpy = (BIT_ARRAY*) malloc(sizeof(BIT_ARRAY));

  word_addr_t num_of_words = nwords(bitarr->num_of_bits);

  cpy->num_of_bits = bitarr->num_of_bits;
  cpy->words = (word_t*) malloc(sizeof(word_t) * num_of_words);

  /* Copy across bits*/
  memcpy(cpy->words, bitarr->words, num_of_words * sizeof(word_t));

  return cpy;
}

/*
void bit_array_copy(BIT_ARRAY* dest, bit_index_t dstindx,
                    BIT_ARRAY* src, bit_index_t srcindx, bit_index_t length)
{

}
*/


/* Enlarge or shrink the size of a bit array*/
/* Shrinking will free some memory if it is large*/
/* Enlarging an array will add zeros to the end of it*/
/* returns 1 on success, 0 on failure*/
char bit_array_resize(BIT_ARRAY* bitarr, bit_index_t new_num_of_bits) {
  bit_index_t old_num_of_bits = bitarr->num_of_bits;
  bitarr->num_of_bits = new_num_of_bits;

  word_addr_t old_num_of_words = nwords(old_num_of_bits);
  word_addr_t new_num_of_words = nwords(new_num_of_bits);

  if(new_num_of_words != old_num_of_words) {
    /* Need to change the amount of memory used*/
    bitarr->words = realloc(bitarr->words, new_num_of_words * sizeof(word_t));

    if(bitarr->words == NULL) {
      /* error - could not allocate enough memory*/
      errno = ENOMEM;
      return 0;
    }
  }

  /* if we are growing - need to zero new bits*/
  if(new_num_of_bits > old_num_of_bits) {
    /* zero entire words first*/
    if(new_num_of_words > old_num_of_words) {
      memset(bitarr->words + old_num_of_words, 0,
             (new_num_of_words - old_num_of_words) * sizeof(word_t));
    }

    /* zero bits on the end of what used to be the last word*/
    unsigned int bits_on_last_word = boffset(old_num_of_bits);

    if(bits_on_last_word > 0) {
      bitarr->words[old_num_of_words - 1] &= BIT_MASK(bits_on_last_word);
    }
  }

  return 1;
}

/**/
/* Logic operators*/
/**/

void bit_array_and(BIT_ARRAY* dest, BIT_ARRAY* src1, BIT_ARRAY* src2) {
  if(dest->num_of_bits == src1->num_of_bits &&
      src1->num_of_bits == src2->num_of_bits) {
    word_addr_t num_of_words = nwords(src1->num_of_bits);

    word_addr_t i;
    for(i = 0; i < num_of_words; i++) {
      dest->words[i] = src1->words[i] & src2->words[i];
    }
  } else {
    /* error*/
    fprintf(stderr, "bit_array.c: bit_array_and() : "
            "dest, src1 and src2 must be of the same length\n");
    exit(EXIT_FAILURE);
  }


}

void bit_array_or(BIT_ARRAY* dest, BIT_ARRAY* src1, BIT_ARRAY* src2) {
  if(dest->num_of_bits == src1->num_of_bits &&
      src1->num_of_bits == src2->num_of_bits) {
    word_addr_t num_of_words = nwords(src1->num_of_bits);

    word_addr_t i;
    for(i = 0; i < num_of_words; i++) {
      dest->words[i] = src1->words[i] | src2->words[i];
    }

  } else {
    /* error*/
    fprintf(stderr, "bit_array.c: bit_array_and() : "
            "dest, src1 and src2 must be of the same length\n");
    exit(EXIT_FAILURE);
  }


}

void bit_array_xor(BIT_ARRAY* dest, BIT_ARRAY* src1, BIT_ARRAY* src2) {
  if(dest->num_of_bits == src1->num_of_bits &&
      src1->num_of_bits == src2->num_of_bits) {
    word_addr_t num_of_words = nwords(src1->num_of_bits);

    word_addr_t i;
    for(i = 0; i < num_of_words; i++) {
      dest->words[i] = src1->words[i] ^ src2->words[i];
    }
  } else {
    /* error*/
    fprintf(stderr, "bit_array.c: bit_array_and() : "
            "dest, src1 and src2 must be of the same length\n");
    exit(EXIT_FAILURE);
  }


}

void bit_array_not(BIT_ARRAY* dest, BIT_ARRAY* src) {
  if(dest->num_of_bits != src->num_of_bits) {
    /* error*/
    fprintf(stderr, "bit_array.c: bit_array_and() : "
            "dest and src1 must be of the same length\n");
    exit(EXIT_FAILURE);
  }

  word_addr_t num_of_words = nwords(dest->num_of_bits);

  word_addr_t i;
  for(i = 0; i < num_of_words; i++) {
    dest->words[i] = ~(src->words[i]);
  }
}


__inline__ word_t get_word(BIT_ARRAY* bitarr, word_addr_t word_index,
                       word_addr_t nwords) {
  if(word_index >= nwords) {
    return 0;
  }

  unsigned int offset;

  if(word_index == nwords - 1 &&
      (offset = boffset(bitarr->num_of_bits)) != 0) {
    /* get masked*/
    return bitarr->words[word_index] & BIT_MASK(offset);
  } else {
    return bitarr->words[word_index];
  }
}

/* Compare two bit arrays by value stored*/
/* arrays do not have to be the same length (e.g. 101 (5) > 00000011 (3))*/
int bit_array_compare(BIT_ARRAY* bitarr1, BIT_ARRAY* bitarr2) {
  word_addr_t nwords1 = nwords(bitarr1->num_of_bits);
  word_addr_t nwords2 = nwords(bitarr2->num_of_bits);

  word_addr_t max_words = MAX(nwords1, nwords2);

  word_addr_t i;
  word_t word1, word2;

  for(i = max_words - 1; i >= 0; i--) {
    word1 = get_word(bitarr1, i, nwords1);
    word2 = get_word(bitarr2, i, nwords2);

    if(word1 > word2) {
      return 1;
    } else if(word1 < word2) {
      return -1;
    }
  }

  return 0;
}

/* Return 0 if there was an overflow error, 1 otherwise*/
char bit_array_add(BIT_ARRAY* dest, BIT_ARRAY* src1, BIT_ARRAY* src2) {
  word_addr_t nwords1 = nwords(src1->num_of_bits);
  word_addr_t nwords2 = nwords(src2->num_of_bits);

  word_addr_t max_words = MAX(nwords1, nwords2);

  word_addr_t dest_words = nwords(dest->num_of_bits);

  char carry = 0;

  word_addr_t i;
  word_t word1, word2;

  for(i = 0; i < max_words; i++) {
    word1 = get_word(src1, i, nwords1);
    word2 = get_word(src2, i, nwords2);

    word_t result = word1 + word2 + carry;
    carry = (result < src1->words[i] || result < src2->words[i]) ? 1 : 0;

    /* Check we can store this result*/
    if(i < dest_words - 1) {
      dest->words[i] = result;
    } else if(i >= dest_words || carry) {
      /* overflow error*/
      return 0;
    } else {
      /* Check last word (i == dest_words-1)*/
      unsigned int bits_on_last_word = boffset(dest->num_of_bits);

      if(bits_on_last_word > 0 && BIT_MASK(bits_on_last_word) < result) {
        /* overflow error*/
        return 0;
      }
    }
  }

  if(carry) {
    dest->words[max_words] = 1;
  }

  /* Zero the rest of dest*/
  for(i = max_words + 1; i < dest_words; i++) {
    dest->words[i] = 0;
  }

  return 1;
}

bit_index_t bit_array_length(BIT_ARRAY* bit_arr) {
  return bit_arr->num_of_bits;
}

/* If there is an overflow, bit array will be set to all 1s and 0 is returned*/
/* Returns 0 if there was an overflow, 1 otherwise*/
char bit_array_increment(BIT_ARRAY* bitarr) {
  word_addr_t num_of_words = nwords(bitarr->num_of_bits);

  char carry = 1;

  word_addr_t i;
  for(i = 0; i < num_of_words - 1; i++) {
    if(bitarr->words[i] + 1 < bitarr->words[i]) {
      /* Carry continues*/
      bitarr->words[i] = 0;
    } else {
      /* Carry is absorbed*/
      bitarr->words[i]++;
      carry = 0;
      break;
    }
  }

  /* Deal with last word*/
  if(carry) {
    unsigned int bits_in_last_word = boffset(bitarr->num_of_bits - 1) + 1;
    word_t mask = BIT_MASK(bits_in_last_word);
    word_t prev_last_word = bitarr->words[num_of_words - 1] & mask;

    /* Increment*/
    bitarr->words[num_of_words - 1] = (prev_last_word + 1) & mask;

    if(prev_last_word == mask) {
      /* overflow*/
      return 0;
    }
  }

  return 1;
}

/* If there is an underflow, bit array will be set to all 0s and 0 is returned*/
/* Returns 0 if there was an underflow, 1 otherwise*/
char bit_array_decrement(BIT_ARRAY* bitarr) {
  word_addr_t num_of_words = nwords(bitarr->num_of_bits);

  word_addr_t i;
  for(i = 0; i < num_of_words - 1; i++) {
    if(bitarr->words[i] > 0) {
      bitarr->words[i]--;

      i--;
      while(i >= 0) {
        bitarr->words[i--] = (unsigned int) ULONG_MAX;
      }

      return 1;
    }
  }

  /* Must subtract from last word*/
  unsigned int bits_in_last_word = boffset(bitarr->num_of_bits - 1) + 1;
  word_t mask = BIT_MASK(bits_in_last_word);
  word_t prev_last_word = bitarr->words[num_of_words - 1] & mask;

  if(prev_last_word == 0) {
    /* underflow*/
    /* number unchanged*/
    return 0;
  } else {
    bitarr->words[num_of_words - 1] = prev_last_word - 1;
    return 1;
  }
}

long bit_array_get_long(BIT_ARRAY* bitarr, bit_index_t start) {
  /* Bounds checking*/
  if(start >= bitarr->num_of_bits) {
    fprintf(stderr, "bit_array.c: bit_array_get_long() - out of bounds error "
            "(index: %lu, length: %lu)\n", start, bitarr->num_of_bits);
    exit(EXIT_FAILURE);
  }

  word_addr_t num_of_words = nwords(bitarr->num_of_bits);

  word_addr_t word_index = bindex(start);
  unsigned int start_offset = boffset(start);

  long result = get_word(bitarr, word_index++, num_of_words) >> start_offset;

  unsigned int offset = WORD_SIZE - start_offset;

  /* 64 bits in a long*/
  while(offset < 64) {
    result |= get_word(bitarr, word_index++, num_of_words) << offset;
    offset += WORD_SIZE;
  }

  return result;
}

int bit_array_get_int(BIT_ARRAY* bitarr, bit_index_t start) {
  /* Bounds checking*/
  if(start >= bitarr->num_of_bits) {
    fprintf(stderr, "bit_array.c: bit_array_get_long() - out of bounds error "
            "(index: %lu, length: %lu)\n", start, bitarr->num_of_bits);
    exit(EXIT_FAILURE);
  }

  word_addr_t num_of_words = nwords(bitarr->num_of_bits);

  word_addr_t word_index = bindex(start);
  unsigned int start_offset = boffset(start);

  int result
  = ((get_word(bitarr, word_index++, num_of_words) >> start_offset) & UINT_MAX);

  unsigned int offset = WORD_SIZE - start_offset;

  /* 32 bits in an int*/
  while(offset < 32) {
    result |= (get_word(bitarr, word_index++, num_of_words) << offset) & UINT_MAX;
    offset += WORD_SIZE;
  }

  return result;
}

char bit_array_get_char(BIT_ARRAY* bitarr, bit_index_t start) {
  /* Bounds checking*/
  if(start >= bitarr->num_of_bits) {
    fprintf(stderr, "bit_array.c: bit_array_get_long() - out of bounds error "
            "(index: %lu, length: %lu)\n", start, bitarr->num_of_bits);
    exit(EXIT_FAILURE);
  }

  word_addr_t num_of_words = nwords(bitarr->num_of_bits);

  word_addr_t word_index = bindex(start);
  unsigned int start_offset = boffset(start);

  char result = get_word(bitarr, word_index, num_of_words) >> start_offset;
  result |= get_word(bitarr, word_index + 1, num_of_words) << (WORD_SIZE - start_offset);

  return result;
}


void bit_array_save(BIT_ARRAY* bitarr, FILE* f) {
  size_t num_of_bytes = nwords(bitarr->num_of_bits);

  fwrite(&num_of_bytes, sizeof(size_t), 1, f);

  fwrite(&bitarr->num_of_bits, sizeof(unsigned int), 1, f);

  fwrite(bitarr->words, sizeof(word_t), num_of_bytes, f);

}


BIT_ARRAY* bit_array_load(FILE* f) {
  BIT_ARRAY* bitarr;
  size_t num_of_bytes;
  size_t read;

  bitarr = malloc(sizeof(BIT_ARRAY));

  read = fread(&num_of_bytes, sizeof(size_t), 1, f);
  if(read != 1){
    fprintf(stderr, "Error reading the file.\n");
    exit(-1);
  }

  bitarr->words = malloc(sizeof(word_t) * num_of_bytes);

  read = fread(&bitarr->num_of_bits, sizeof(unsigned int), 1, f);
  if(read != 1){
    fprintf(stderr, "Error reading the file.\n");
    exit(-1);
  }

  read = fread(bitarr->words, sizeof(word_t), num_of_bytes, f);
  if(read != 1){
    fprintf(stderr, "Error reading the file.\n");
    exit(-1);
  }

  return bitarr;
}


/* Thread-safe */
void parallel_bit_array_concat_from_to(BIT_ARRAY* dest, BIT_ARRAY* source, unsigned long shift, unsigned long from, unsigned long len) {

  if((source->num_of_bits <= 0) || (len <= 0)) return;

  /* unsigned int desp, corr, sdesp, scorr, num_shifts, i; */
  /* unsigned int len_aux = len; */
  unsigned long desp, corr, sdesp, scorr, num_shifts, i;
  unsigned long len_aux = len;
  word_t *d, *f;
  word_t old_d, new_d; /* Thread safe shifting*/
  unsigned long aux = 0;
  /*  unsigned int aux = 0;*/

  d = &dest->words[nwords(shift + 1) - 1];
  f = &source->words[nwords(from + 1) - 1];

  desp = shift % WORD_SIZE; /* Position to start to write in dest (lower value: 0, upper value: 31)*/
  corr = WORD_SIZE - desp; /* Number of free positions in dest including the start position (lower value: 1, upper value: 32)*/

  sdesp = from % WORD_SIZE; /* Position to start to read from the source (lower value: 0, upper value: 31)*/
  scorr = WORD_SIZE - sdesp; /* Number of positions since the start position of the source, including the start position (lower value:1, upper value: 32)*/


  /* The first word needs to be copied atomically*/

  if(sdesp > 0) {
    if(corr > scorr) {
      if(len_aux < scorr) {
        do{
          old_d = *d;
          new_d = old_d | (((*f << (scorr - len_aux)) >> (WORD_SIZE - len_aux)) << desp);
        }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
        return;
      }

      do{
        old_d = *d;
        new_d = old_d | ((*f >> sdesp) << desp);
      }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
      desp = desp + scorr;
      corr = WORD_SIZE - desp;
    }
    else if(corr < scorr) {
      if(len_aux < corr) {
        do{
          old_d = *d;
          new_d = old_d | (((*f << (scorr - len_aux)) >> (WORD_SIZE - len_aux)) << desp);
        }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
        return;
      }

      do{
        old_d = *d;
        new_d = old_d | ((*f >> sdesp) << desp);
      }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));

      *d++;

      if(len_aux < scorr) {
        if((WORD_SIZE  - len_aux + corr) < WORD_SIZE ) {
          do{
            old_d = *d;
            new_d = old_d | ((*f << (scorr - len_aux)) >> (WORD_SIZE  - len_aux + corr));
          }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
        }
        return;
      }

      if((sdesp + corr) < WORD_SIZE) {
        do{
          old_d = *d;
          new_d = old_d | (*f >> (sdesp + corr));
        }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
      }
      desp = scorr - corr;
      corr = WORD_SIZE - desp;
    }
    else if(corr == scorr) { /* Implies desp == sdesp*/
      if(len_aux < corr) {
        if((scorr - len_aux + sdesp) < WORD_SIZE) {
          do{
            old_d = *d;
            new_d = old_d | (((*f << (scorr - len_aux)) >> (scorr - len_aux + sdesp)) << desp);
          }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
        }
        return;
      }

      do{
        old_d = *d;
        new_d = old_d | ((*f >> sdesp) << desp);
      }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));

      *d++;
      desp = (desp + corr)%WORD_SIZE; /* Always desp = 0 in this case;*/
      corr = WORD_SIZE - desp;
    }

    len_aux -= scorr;
  }


  /* Until here, we force the case sdesp=0. Therefore, it's necessary redefine*/
  /* desp and corr, depending the case*/

  /* All words, except the first and the last, can be copied without using CAS operations*/

  num_shifts = ((len_aux + desp - 1)/WORD_SIZE) + 1;

  for(i = 0; i < num_shifts-1; i++) {
    /* General copy of bits using "desp" and "corr"*/


    if(len_aux <= 0) {
      return;
    }

    if((i>0) || (sdesp > 0)) /* It implies that it satisfied the previous condition (sdesp > 0)*/
      *f++;

    if(len_aux <= corr) {
      *d |= (((*f << ( WORD_SIZE - len_aux)) >> (WORD_SIZE - len_aux)) << desp);
      return;
    }

    /* Strange behavior with desp = WORD_SIZE (Is it necessary?)*/
    *d |= (*f << desp);

    len_aux -= corr;
    *d++;

    if(len_aux <= desp) {
      *d |= ((*f << (desp - len_aux)) >> (WORD_SIZE - len_aux));
      return;
    }

    /* Strange behavior with corr = WORD_SIZE*/
    if(corr < WORD_SIZE)
	    *d |= (*f >> corr);

    len_aux -= desp;

  }


  /* The last word needs to be copied atomically*/
  
  if(len_aux <= 0) {
    return;
  }
  
  if((i>0) || (sdesp > 0)) /* It implies that it satisfied the previous condition (sdesp > 0)*/
    *f++;
  
  if(len_aux <= corr) {
    do{
      old_d = *d;
      new_d = old_d | (((*f << ( WORD_SIZE - len_aux)) >> (WORD_SIZE - len_aux)) << desp);
    }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
    return;
  }
  
  do{
    old_d = *d;
    new_d = old_d | (*f << desp);
  }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
  
  len_aux -= corr;
  *d++;
  
  if(len_aux <= desp) {
    do{
      old_d = *d;
      new_d = old_d | ((*f << (desp - len_aux)) >> (WORD_SIZE - len_aux));
    }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
    return;
  }
  
  /* Strange behavior with corr = WORD_SIZE*/
  if(corr < WORD_SIZE) {
    do{
      old_d = *d;
      new_d = old_d | (*f >> corr);
    }while(!(__atomic_compare_exchange_n(d, &old_d, new_d, 0, __ATOMIC_SEQ_CST, __ATOMIC_SEQ_CST)));
  }
}
