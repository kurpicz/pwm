#ifndef DYNSTR_H
#define DYNSTR_H

typedef struct _dynstr dynstr;

/**
 * @file dynstr.h
 * @author Paulo Fonseca
 * 
 * @brief Dynamic/variable-length string (a.k.a. String buffer) ADT.
 */


/**
 * @brief Creates an empty dynamic string.
 * @param init_capacity The initial capacity (# of chars)
 */
dynstr *dstr_new(size_t init_capacity);


/**
 * @brief Destructor.
 */
void dstr_free(dynstr *dstr);


/**
 * @brief Creates a new dynamic string from a source static string @src.
 * @warning The source string contents are simply copied onto the dynamic
 *          string and the former is left untouched.
 */
dynstr *dstr_new_from_str(char *src);


/**
 * @brief Returns the "logical" length of a given dynamic string.
 */
size_t dstr_len(dynstr *dstr);


/**
 * @brief Returns the capacity (i.e. the "physical" length)
 *        of a given dynamic string.
 */
size_t dstr_capacity(dynstr *dstr);


/**
 * @brief Returns the character at a given position.
 */
char dstr_get(dynstr *dstr, size_t pos);


/**
 * @brief Sets (overwrites) the character of a given position @p pos to @p c.
 */
void dstr_set(dynstr *dtsr, size_t pos, char c);


/**
 * @brief Appends a copy of the contents of a static string @p suff. 
 * @warning The source string @p suff is left untouched.
 */
void dstr_append(dynstr *dstr, char *suff);


/**
 * @brief Appends a character @p c.
 */
void dstr_append_char(dynstr *dstr, char c);


/**
 * @brief Returns a read-only reference to the current internal static string.
 *
 * @warning The internal string can change between calls and the returned 
 *          reference can become NULL or invalid.
 */
const char *dstr_as_str(dynstr *dstr);


/**
 * @brief Detaches and returns the current internal static string after 
 *        trimming (removal of trailing unused positions).
 * 
 * @see cstr_trim
 * @warning After this operation, the dynamic string @p dstr is destroyed.
 */
char *dstr_detach(dynstr *dstr);

#endif