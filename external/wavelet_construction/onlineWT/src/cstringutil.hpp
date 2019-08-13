#ifndef CSTRINGUTIL_H
#define CSTRINGUTIL_H

#include <stdio.h>
#include <stdlib.h>

#include "cocadautil.hpp"

/**
 * @file cstringutil.h
 * @author Paulo Fonseca
 *
 * @brief Plain C strings (char *) manipulation functions.
 */


/**
 * @brief Creates a new C string of a given length. 
 *        The physical length of the string will be @p len+1 and the 
 *        string will be initially filled with '\0'.
 *        The returned array will thus be capable of storing strings whose
 *        strlen varies from 0 to @p len.
 */
char *cstr_new(size_t len);


/**
 * @brief Sets @p str[j]=@p c for @p from <= j < @p to.
 */
void cstr_fill(char *str, size_t from, size_t to, char c);


/**
 * @brief Sets the first @p len positons of @p str to '\0'.
 */
void cstr_clear(char *str, size_t len);


/**
 * @brief 'Trims' the string @p str to @p str[@p from..@p to-1]. 
 *        The memory used by the parts of the string out of this interval
 *        will be freed. As part of the operation, the remaining 'trimmed'
 *        portion may be relocated.
 * @return The address of the trimmed string.
 */
char *cstr_trim(char *str, size_t from,  size_t to);


/**
 * @brief Same as cstr_trim(@p str, 0, @plen)
 */
char *cstr_trim_to_len(char *str, size_t len);


/**
 * @brief Reverts a string in place.
 */
void cstr_revert(char *str, size_t len);


/**
 * @brief Converts a size_t to a string.
 * @param dest The target string with enough space
 *        ( ceil(log(sizeof(size_t), b) + 1), where b is the base of the
 *        conversion
 * @param val the value to be converted
 * @param base 'b' (binary), 'o' (octal), 'd' (decimal=default), 
 *        'h' (hexadecimal)
 */
void sizet_to_cstr(char *dest, size_t val, char base);


/**
 * @brief Converts a pointer (memory address) to its hex value representation.
 * @param ptr  The pointer to be converted
 * @param dest The destination string or NULL, in which case, a new string 
 *             is allocated.
 * @return The value of @p dest, if it is not NULL, or the address of the 
 *         newly created string.
 */
void *ptr_to_cstr(void *ptr, char *dest);


#endif

