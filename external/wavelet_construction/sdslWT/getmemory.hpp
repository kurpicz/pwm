
#include <unistd.h>
#include <sys/resource.h>
#include <stdio.h>



size_t getPeakRSS() {
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage);
	return (size_t)(rusage.ru_maxrss * 1024L);
}

