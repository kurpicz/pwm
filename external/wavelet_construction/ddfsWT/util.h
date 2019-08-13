#ifndef UTIL_H_S3NZO2TM
#define UTIL_H_S3NZO2TM

#ifdef CHARS
typedef unsigned char symbol;
#elif SHORT
typedef unsigned short symbol;
#else
typedef unsigned int symbol;
#endif

symbol* read_text_from_file(const char* fn, unsigned long* n);

#endif /* end of include guard: UTIL_H_S3NZO2TM */
