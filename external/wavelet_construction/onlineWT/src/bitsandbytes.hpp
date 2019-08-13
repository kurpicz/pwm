#ifndef BITSANDBYTES_H
#define BITSANDBYTES_H

#include <endian.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>


/**
 * @file bitsandbytes.h
 * @author Paulo Fonseca
 *
 * @brief General definitions and functions necessary for low level
 *        bit/byte operatios.
 *
 * @note ATTENTION: The code herein is not portable. Check the docs for
 *       'porting' instructions and issues.
 */


/**
 * @brief A required unsigned byte type.
 *
 * The C11 standard  IEC 9899:2011 defines a "byte" as an addressable unit
 * of data storage large enough to hold any member of the basic character
 * set of the execution environment.
 *
 * It also defines a char as "single-byte" character and so sizeof(char)
 * should always return 1.
 *
 * Moreover, the standard library <limits.h> defines CHAR_BIT to be the
 * number of bits for smallest object that is not a bit-field ("byte")
 * and specifies a minimum size of 8 (eight).
 *
 * Although a byte may be composed of more that eight bits,
 * ISO Norm IEC 80000-13:2008 (item 13.9 c) suggests that the word "byte"
 * be defined as a synonymm of octet, i.e. a sequence of eight bits.
 *
 * A byte_t type is therefore defined as an alias for unsigned char.
 * A BYTESIZE constant is defined as a CHAR_BIT synonym, and
 * accordingly a maximum value constant BYTE_MAX is defined as
 * UCHAR_MAX synonym.
 */
typedef unsigned char byte_t;

#define BYTESIZE CHAR_BIT

#define BYTE_MAX UCHAR_MAX

/*
 * Although I have attempted to make this code robust to any value of
 * BYTESIZE, it has not been tested at all on any architecture with
 * a value other than 8.
 *
 * The following condition has thus been included as a safeguard.
 *
 * If this is removed to support larger bytes, the byte masks below
 * need also to be modified.
 */
#if BYTESIZE!=8
#error "This code has only been tested on architectures with bytes of 8 bits"
#endif

/*
 * 8-bit Byte masks:
 *
 * Least significant bits masks _LSBMASK[j] = 0^(8-j)1^j
 */
static byte_t _LSBMASK[9] = { 0x00, 0x01, 0x03, 0x07, 0x0f,
                              0x1f, 0x3f, 0x7f, 0xff
                            };
#define LSBMASK(N)  (_LSBMASK[N])
/*
 * Most significant bits masks _MSBMASK[j]=1^j0^(8-j)
 */
static byte_t _MSBMASK[9] = { 0x00, 0x80, 0xc0, 0xe0, 0xf0,
                              0xf8,0xfc, 0xfe, 0xff
                            };
#define MSBMASK(N)  (_MSBMASK[N])


/*
 * Some bitwise operations in this code can benefit from bit parallelism
 * at the processor level. It is useful to have the WORD_BIT/BYTEWORDSIZE
 * constants, which correspond to the size of a machine word in bits/bytes
 * known at compile time.
 *
 * We try to determins the values of these constants based on the processor
 * type, as shown in
 * http://nadeausoftware.com/articles/2012/02/c_c_tip_how_detect_processor_type_using_compiler_predefined_macros
 *
 * Since I have only tested the code on x86 architectures, I have included
 * a guard to prevent compilation in other machines and force the developer
 * to be aware of this and take the appropriate measures to make the code
 * compatible with other platforms.
 */

enum _cl_processor { x86_32bits, x86_64bits };

#if defined(__x86_64__) || defined(_M_X64)
#define PROCESSOR x86_64bits
#define WORD_BIT 64
#define BYTEWORDSIZE 8
typedef uint64_t word_t;
#elif defined(__i386) || defined(_M_IX86)
#define PROCESSOR x86_32bits
#define WORD_BIT 32
#define BYTEWORDSIZE 4
typedef uint32_t word_t;
/*
 * You might want to add support for your platform here
 */
#else
#error "This code has only been tested on 32/64 bits x86 platforms."
#endif

//#define WORD_MAX (~((word_t)0))


/*
 * A constant ENDIANNESS is required to indicate byte endianess.
 *
 * As of 2016, there is no portable, standard way to check for platform
 * endianess at compile time. Since much of the code on this library depends on
 * packing numbers in bit (byte) arrays, there is a necessity of knowing
 * the details of the binary representation for marshalling/unmarshalling.
 * This could be done at runtime but it would require some tests which would
 * degrade performance. Instead, we choose to sacrifice some portability and
 * rely on some specific compilers and platforms. In particular, this code
 * only works for if byte endianness is BIG or LITTLE, as defined in the
 * endian.h.
 *
 * Some of the code herein assumes ENDIANESS to be either LITTLE (0) of
 * BIG (1). The code will not compile as is on platforms that use any other
 * byte ordering alternative.
 */
enum _cl_endianness { BIG, LITTLE };
#define LITTLE 0
#define BIG 1
#ifdef __BYTE_ORDER
#if __BYTE_ORDER == __BIG_ENDIAN
#define ENDIANNESS BIG
#elif __BYTE_ORDER == __LITTLE_ENDIAN
#define ENDIANNESS LITTLE
#else
#error "This code requires BIG or LITTLE byte endianess.\n\
See release notes for porting issues."
#endif
#else
#error "Unable to identify byte endianness.\n\
This code requires BIG or LITTLE byte endianess.\n\
See release notes for porting issues."
#endif


/**
 * @brief Converts a byte to a binary string.
 */
void byte_to_str(byte_t b, char *dest);


/**
 * @brief Reverts the bits of a byte in-place.
 */
void byte_reverse(byte_t *b);


/**
 * @brief Same as byte_bitcount(n, 0)
 * @see byte_bitcount
 */
size_t byte_bitcount0(byte_t n);


/**
 * @brief Same as byte_bitcount(n, 1)
 * @see byte_bitcount
 */
size_t byte_bitcount1(byte_t n);


/**
 * @brief Returns the number of bits with value==@p bit of a given byte.
 */
size_t byte_bitcount(byte_t n, byte_t bit);


/**
 * @brief Same as uint16_bitcount(n, 0)
 * @see uint16_bitcount
 */
size_t uint16_bitcount0(uint16_t n);


/**
 * @brief Same as uint16_bitcount(n, 1)
 * @see uint16_bitcount
 */
size_t uint16_bitcount1(uint16_t n);


/**
 * @brief Returns the number of bits with value==@p bit of a given 16-bit uint.
 */
size_t uint16_bitcount(uint16_t n, byte_t bit);


/**
 * @brief Same as uint32_bitcount(n, 0)
 * @see uint32_bitcount
 */
size_t uint32_bitcount0(uint32_t n);


/**
 * @brief Same as uint32_bitcount(n, 1)
 * @see uint32_bitcount
 */
size_t uint32_bitcount1(uint32_t n);


/**
 * @brief Returns the number of bits with value==@p bit of a given 32-bit uint.
 */
size_t uint32_bitcount(uint32_t n, byte_t bit);


/**
 * @brief Same as uint64_bitcount(n, 0)
 * @see uint64_bitcount
 */
size_t uint64_bitcount0(uint64_t n);


/**
 * @brief Same as uint64_bitcount(n, 1)
 * @see uint64_bitcount
 */
size_t uint64_bitcount1(uint64_t n);


/**
 * @brief Returns the number of bits with value==@p bit of a given 64-bit uint.
 */
size_t uint64_bitcount(uint64_t n, byte_t bit);


/**
 * @brief Same as byte_rank(@p b, @p pos, 0)
 * @see byte_rank
 */
size_t byte_rank0(byte_t b, size_t pos);


/**
 * @brief Same as byte_rank(@p b, @p pos, 1)
 * @see byte_rank
 */
size_t byte_rank1(byte_t n, size_t i);


/**
 * @brief Computes rank_@p bit(@p b, @p pos) = # positions j<=@p pos
 * s.t. @p b[j]==@p bit, 
 * where @p b[j] denotes the jth bit of byte @p b from the left.
 * If i>=BYTESIZE, returns the total number of positions with value == @p bit.
 */
size_t byte_rank(byte_t b, size_t pos, byte_t bit);


/**
 * @brief Same as byte_select(@p b, @p rank, 0)
 * @see byte_select
 */
size_t byte_select0(byte_t b, size_t rank);


/**
 * @brief Same as byte_select(@p b, @p rank, 1)
 * @see byte_select
 */
size_t byte_select1(byte_t b, size_t rank);


/**
 * @brief Computes select_@p bit(@p b, @p rank) = j s.t.
 * @p b[j]==@p bit and rank_@p bit(@p b, j)=@p rank,
 * where @p b[j] denotes the jth bit of byte @p b from the left.
 * If no such position exists, returns BYTESIZE.
 */
size_t byte_select(byte_t b, size_t rank, byte_t bit);

#endif