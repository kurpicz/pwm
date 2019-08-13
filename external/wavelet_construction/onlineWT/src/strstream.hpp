#ifndef STRSTREAM_H
#define STRSTREAM_H

#include <stddef.h>

#include "cocadautil.hpp"

/**
 * @file strstream.h
 * @author Paulo Fonseca
 * 
 * @brief String stream.
 */


/**
 * String stream type
 */
typedef struct _strstream strstream;


/**
 * @brief Opens a stream for a in-memory source string.
 * @param str The source string.
 * @param slen The source string length.
 */
strstream *strstream_open_str(char *str, size_t slen);


/**
 * @brief Opens a stream for a source text file.
 */
strstream *strstream_open_file(char *filename);


/**
 * @brief Resets the stream, i.e. moves cursor to initial position.
 */
void strstream_reset(strstream *sst);


/**
 * @brief Tests whether a stream has reached its end.
 */
bool strstream_end(strstream *sst);


/**
 * @brief Reads the next char from a stream.
 * @returns The next character as an int, or EOF if the stream has 
 *          reached its end.
 * 
 * Example of usage: 
 * @code
 * strstream *fsst = strstream_open_file(filename);
 * for (int c; (c=strstream_getc(fsst)) != EOF;)
 *     printf ("Read c=%c\n", (char)c);
 * strstream_close(fsst);
 * @endcode
 */
int strstream_getc(strstream *sst);


/**
 * @brief Attempts to read the next @p n chars into the string *dest. 
 *        Less than @p n characters can be read if the stream reaches its end. 
 * @returns The number of chars actually read. 
 */
size_t strstream_reads(strstream *sst, char *dest, size_t n);


/**
 * @brief Closes the stream and disposes the stream object.
 */
void strstream_close(strstream *sst);


#endif
