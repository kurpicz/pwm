#ifndef __TIMING_H__
#define __TIMING_H__

#include <unistd.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/times.h>

double ticks;
double secs;
struct tms t1,t2;
struct timeval tv1, tv2;

void startTimer() {
	ticks= (double)sysconf(_SC_CLK_TCK);
	gettimeofday(&tv1,NULL);
	times (&t1);
}


double timeFromBegin() {
	times (&t2);
	return ((t2.tms_utime+t2.tms_stime)-(t1.tms_utime+t1.tms_stime))/ticks;
}

long realTimeFromBegin() {
	gettimeofday(&tv2,NULL);
	return (tv2.tv_sec-tv1.tv_sec)*1000000+(tv2.tv_usec-tv1.tv_usec);
}

#endif
