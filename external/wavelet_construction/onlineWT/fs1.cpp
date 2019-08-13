#include <stdlib.h>
#include <stdio.h>

#include "src/bytearray.hpp"
#include "src/cstringutil.hpp"
#include "src/strstream.hpp"
#include "src/wavtree.hpp"

void offline(char *filename)
{
	strstream *sst = strstream_open_file(filename);
    byte_t *ab_table = bytearr_new(256);
	size_t n=0, m=0;
    for (int c=0; (c=(strstream_getc(sst)))!=EOF; n++) {
    	if (!ab_table[c]) {
			ab_table[c]=1; 
			m++;
		}
	}
    char *ab_str = cstr_new(m);
    size_t k = 0;
    for (int i=0; i<256; i++)
        if (ab_table[i]==1) ab_str[k++] = (char)i;
    //printf("ALPHABET IS %s", ab_str);
    alphabet *ab = alphabet_new(m, ab_str);
	printf("creating wavelet_tree (offline)\n");
    wavtree *wt = wavtree_new_from_stream(ab, sst, WT_BALANCED);
	printf("wavelet_tree created\n");
}


int main(int argc, char *argv[])
{
	offline(argv[1]);
	return 0;
}
