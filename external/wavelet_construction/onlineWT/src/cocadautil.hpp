#ifndef COCADAUTIL_H
#define COCADAUTIL_H

#include <stdbool.h>
#include <stdlib.h>

/**
 * @file cocadautil.h
 * @author Paulo Fonseca
 *
 * @brief Generic definitions, macros and utility functions.
 */


/**
 * NULL pointer constant.
 */
#if !defined(NULL)
#define NULL ((void *)0)
#endif

/**
 * Boolean data type.
 */

#define MIN( A, B )   ( (A)<(B) ? (A) : (B) )
#define MAX( A, B )   ( (A)>(B) ? (A) : (B) )

#define MIN3( A, B, C ) ((A) < (B) ? MIN(A, C) : MIN(B, C))
#define MAX3( A, B, C ) ((A) > (B) ? MAX(A, C) : MAX(B, C))

#define NEW( TYPE )   ((TYPE*)(malloc(sizeof(TYPE))))
#define FREE( PTR )   ({if ((PTR)!=NULL) free((PTR));})

/**
 * The maximum value of a pointer
 */
#define PTR_MAX SIZE_MAX



#endif
