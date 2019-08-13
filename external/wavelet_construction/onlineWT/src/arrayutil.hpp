#ifndef ARRAYUTIL_H
#define ARRAYUTIL_H

#include <stdlib.h>

#include "bitsandbytes.hpp"
#include "cocadautil.hpp"

/**
 * @file arrayutil.h
 *
 * @brief Assorted array utility macros and functions.
 */

/**
 * @brief Allocates a new array of N elements of a given TYPE.
 */
#define NEW_ARRAY( TYPE, N ) ((TYPE*)(malloc((N)*sizeof(TYPE))))

/**
 * @brief Fills an array ARR from position =FROM to position <TO with value VAL.
 */
#define FILL_ARRAY( ARR, FROM, TO, VAL ) \
    for(size_t _i=(FROM), _to=(TO); _i<_to; (ARR)[_i++]=(VAL))

/**
 * @brief Copies N elements from an array SRC from position =FROMSRC
 *        into an array DEST from position =FROMDEST.
 */
#define COPY_ARRAY( DEST, FROMDEST, SRC, FROMSRC, N ) \
    for(size_t _i=0, _n=(N), _fs=(FROMSRC), _fd=(FROMDEST); _i<_n; _i++) \
        (DEST)[_fd+_i]=(SRC)[_fs+_i]

/*
char*               char_arr_new       (size_t len);
short*              short_arr_new      (size_t len);
int*                int_arr_new        (size_t len);
long*               long_arr_new       (size_t len);
long long*          longlong_arr_new   (size_t len);
unsigned char*      uchar_arr_new      (size_t len);
unsigned short*     ushort_arr_new     (size_t len);
unsigned int*       uint_arr_new       (size_t len);
unsigned long*      ulong_arr_new      (size_t len);
unsigned long long* ulonglong_arr_new  (size_t len);
byte_t*             byte_arr_new       (size_t len);
size_t*             size_t_arr_new     (size_t len);
float*              float_arr_new      (size_t len);
double*             double_arr_new     (size_t len);
long double*        longdouble_arr_new (size_t len);

void char_arr_fill       (int *arr, size_t to, size_t from, int val);
void short_arr_fill      (short *arr, size_t to, size_t from, short val);
void int_arr_fill        (int *arr, size_t to, size_t from, int val);
void long_arr_fill       (long *arr, size_t to, size_t from, long val);
void longlong_arr_fill   (long long *arr, size_t to, size_t from, long long val);
void uchar_arr_fill      (unsigned char *arr, size_t to, size_t from, unsigned char val);
void ushort_arr_fill     (unsigned short *arr, size_t to, size_t from, unsigned short val);
void uint_arr_fill       (unsigned int *arr, size_t to, size_t from, unsigned int val);
void ulong_arr_fill      (unsigned long *arr, size_t to, size_t from, unsigned long val);
void ulonglong_arr_fill  (unsigned long long *arr, size_t to, size_t from, unsigned long long val);
void byte_arr_fill       (byte_t *arr, size_t to, size_t from, byte_t val);
void size_t_arr_fill     (size_t *arr, size_t to, size_t from, size_t val);
void float_arr_fill      (float *arr, size_t to, size_t from, float val);
void double_arr_fill     (double *arr, size_t to, size_t from, double val);
void longdouble_arr_fill (long double *arr, size_t to, size_t from, long double val);
*/
#endif
