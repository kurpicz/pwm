#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "src/bytearray.hpp"
#include "src/cstringutil.hpp"
#include "src/strstream.hpp"
#include "src/wavtree.hpp"

void online(char *filename) {
	clock_t start, finish;
	strstream *sst; 
	sst = strstream_open_file(filename);
	wavtree *wt;
	printf("creating wavelet_tree (online)\n");
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
	online(argv[1]);
	return 0;
}
