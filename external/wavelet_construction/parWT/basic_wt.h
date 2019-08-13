#include "bitrank/bitrankw32int.h"
#include "bitrank/bit_array.h"
#include "util.h"


BIT_ARRAY** wt_create(symbol*, unsigned long, unsigned int);
unsigned int wt_access(bitRankW32Int**, unsigned int, unsigned int);
unsigned int wt_rank(bitRankW32Int**, unsigned int, unsigned int, unsigned int);
unsigned int wt_select(bitRankW32Int**, unsigned int, unsigned int, unsigned int);
void wtree_to_string(bitRankW32Int**, unsigned int);
