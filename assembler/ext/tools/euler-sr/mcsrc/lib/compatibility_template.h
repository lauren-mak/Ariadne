/***************************************************************************
 * Title:          compatibility_template.h and compatibility.h
 * Author:         Glenn Tesler
 * Created:        2009
 * Last modified:  06/02/2010
 *
 * Copyright (c) 2009-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

// Define data types, constants, and functions that may vary between platforms,
// word sizes, etc. 

#ifndef _COMPATIBILITY_H_
#define _COMPATIBILITY_H_

// Some compilers (incl. g++) need _STDC_*_MACROS set to bring in
// certain definitions when inttypes.h and stdint.h are loaded.
// You must make sure to put '#include "compatibility.h"' *before*
// #include <inttypes.h>, <stdint.h>, or anything that #includes them, since
// if they're already included, "#include <inttypes.h>" below will be ignored.

#ifndef __STDC_CONSTANT_MACROS
   // needed for macros like UINT64_C(x)
#  define __STDC_CONSTANT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
   // needed for macros like PRIu64
#  define __STDC_FORMAT_MACROS
#endif
#ifndef __STDC_LIMIT_MACROS
   // needed for macros like UINT64_MAX
#  define __STDC_LIMIT_MACROS
#endif

#include <limits.h>
#include <stdlib.h>
#include <inttypes.h>


/***************************************************************************
 * Code that is customized by CreateCompatibility.cpp
 * It will be injected after the line "// --AUTOGEN--"
 ***************************************************************************/

// The "--AUTOGEN--" line below will generate code similar to this:
// #define _SIZET_IS_LONG_
// #define PRI_PID "%d"
// #define PRI_SZT "%zu"
// #define PRI_SSZT "%zd"
// #define LongTupleWord_less_size_t
//
// The _SIZET_IS_LONG_ line may be replaced by _SIZET_IS_INT_, depending on
// platform.  If there are other possibilities, they will have to be added
// in CreateCompatibility.cpp and in this file.
//
// Usage of PRI_*: It would be better to use C++ style "std::cout << ...",
// but if you must use printf:
//
//      pid_t pid = getpid();
//      printf("pid = " PRI_PID "\n", pid);
//
//      std::vector<E> *edges;
//      ...
//      printf("number of edges = " PRI_SZT "\n", edges->size());
//
//      size_t x = ...;
//      ssize_t y = ...;
//      printf("x = " PRI_SZT "  y = " PRI_SSZT "\n", x, y);

// --AUTOGEN--



/***************************************************************************
 * Data types and widths
 ***************************************************************************/

#if !defined(_SIZET_IS_LONG_) && !defined(_SIZET_IS_INT_)
#  define _SIZET_IS_LONG_
#endif

#define _NBITS_T_(_t) (sizeof(_t) * CHAR_BIT)

#define _IS_SIGNED_T_(_t) (! ((_t) 0 < (_t) -1))
#define _MIN_T_(_t) ((_t) \
										     (_IS_SIGNED_T_(_t) \
													? ((~ (_t) 0) << (_NBITS_T_(_t) - 1))	\
													: (_t) 0))
#define _MAX_T_(_t) ((_t) ((~ (_t) 0) - _MIN_T_(_t)))

//#ifndef PTRDIFF_MAX
//# define PTRDIFF_MAX _MAX_T_(ptrdiff_t)
//#endif
//#ifndef OFF_MAX
//# define OFF_MAX _MAX_T_(off_t)
//#endif
#ifndef SIZE_MAX
#   define SIZE_MAX _MAX_T_(size_t)
#endif
#ifndef SSIZE_MAX
#   define SSIZE_MAX _MAX_T_(ssize_t)
#endif
#ifndef SIZE_BITS
#   define SIZE_BITS _NBITS_T_(size_t)
#endif

#define _INT_BITS_ _NBITS_T_(int)


/***************************************************************************
 * LongTuple specification
 *
 * A "LongTuple" is an array used to store a k-mer in a bitfields, 2 bits
 * per nucleotide, defined as
 *
 *     typedef LongTupleWord LongTuple[LongTupleNumWords];
 *
 * To be flexible with different platforms, these items must be declared here:
 *
 * LongTupleWord: any unsigned integer type, provided:
 *      has an integer number of bytes
 *      has an even number of bits
 * Those are true for all known types besides bitfields.
 * A type of the same size as size_t (or a multiple of it) is probably best,
 * but one can also use smaller types and then form an array of them.
 *
 * LongTupleNumWords: A small integer, at least 1
 *
 * LongTupleCountBits: Set to 0 for now.
 * Reserved for possible future use to put the count in a bitfield
 * (in CountedIntegralTuple), as was previously done.
 ***************************************************************************/

/* Note: To increase the range of k:
 * The bound on k is given by
 *    2*(k+1) <= LongTupleNumWords * sizeof(LongTupleWord) * CHAR_BIT
 *    2*(k+1) <= 2 * 8 * 8 = 128
 *    k <= 63
 * To increase the range of k, edit the defintions of LongTupleWord
 * and LongTupleNumWords below.
 *
 * The default settings
 *    typedef uint64_t LongTupleWord;
 *    #define LongTupleNumWords 2
 * give a structure that is 64*2=128 bits long.
 * Each nucleotide is stored in 2 bits, allowing k+1 <= 64, so k <= 63.
 */

typedef uint64_t LongTupleWord;
#define LongTupleNumWords 2
typedef LongTupleWord LongTuple[LongTupleNumWords];

// One of the following will be generated by CreateCompatibility.cpp:
// #define LongTupleWord_less_size_t // for sizeof(LongTupleWord) < sizeof(size_t)
// #undef LongTupleWord_less_size_t // for sizeof(LongTupleWord) >= sizeof(size_t)

// TESTING:
// typedef uint8_t LongTupleWord;
// #define LongTupleNumWords 8
// typedef LongTupleWord LongTuple[LongTupleNumWords];

// LongTupleMaxK: maximum K user can enter on command line
// But internally, we use both K (vertices) and K+1 (edges), so need
// to be able to store (K+1)-mers as well.
//    divide by 2: 2 bits per nucleotide
//    subtract 1: both K and K+1 should fit in a LongTuple
// LongTupleCountBits:
//    RESERVED: may use some bits for bit field for count, like we used to do
#define LongTupleCountBits 0
#define LongTupleMaxK (((_NBITS_T_(LongTuple) - LongTupleCountBits) / 2) - 1)

// usage:
//     LongTuple longtuple = ...;
//     printf("longtuple = " PRId_longTuple "\n", longtuple);  // decimal
//     printf("longtuple = " PRIx_longTuple "\n", longtuple);  // hex
#define PRId_longTuple PRId64
#define PRIx_longTuple PRIx64


/***************************************************************************
 * Protect some types from the automated scripts that would otherwise change
 *   int -> ssize_t
 *   unsigned int -> size_t
 *
 * and the reverse
 *   ssize_t -> int
 *   size_t -> unsigned int
 ***************************************************************************/

#define _INT_ int
#define _UINT_ unsigned int

#define _SZT_ size_t
#define _SSZT_ ssize_t

#define _LONG_ long
#define _ULONG_ unsigned long


/***************************************************************************
 * C functions whose names change for different types of parameters
 * int, long, (s)size_t, etc.
 ***************************************************************************/

#if defined(_SIZET_IS_LONG_)
#   define szabs labs
#   define atosz atol
#elif defined(_SIZET_IS_INT_)
#   define szabs abs
#   define atosz atoi
#else
#   error "Need to define type of size_t in compatibility.h"
#endif


/***************************************************************************
 * Additional constants and functions with platform differences
 ***************************************************************************/

// Define RANDOM_MAX to be max of random() function:
//   On 64-bit SPARC, RAND_MAX = 2^15 - 1 is max of rand() function,
//   but random() goes up to 2^31 - 1, and doesn't have a named constant.
//   On other platforms tested so far, random() goes up to RAND_MAX.
#if defined(__sun__)
#  define RANDOM_MAX 0x7fffffffL
#else
#  define RANDOM_MAX RAND_MAX
#endif


#if defined(__sun__) && !defined(NEED_LIBSUNMATH)

#  ifndef isnan
#     define isnan(x) \
          (sizeof (x) == sizeof (long double) ? isnan_ld (x) \
           : sizeof (x) == sizeof (double) ? isnan_d (x) \
           : isnan_f (x))
      static inline int isnan_f  (float       x) { return x != x; }
      static inline int isnan_d  (double      x) { return x != x; }
      static inline int isnan_ld (long double x) { return x != x; }
#  endif

#  define isNan isnan

#else
#  define isNan std::isnan
#endif


/***************************************************************************
 * Bitfields
 *
 * Bitfields may have platform or wordsize dependencies, but are still
 * better to define in a class/structure definition in the relevant
 * file rather than here.  Please follow these conventions:
 *
 * 1. Mark *all* definitions of classes/structures containing bitfields
 *    with this comment:
 *       "// HAS BITFIELD"
 *
 * 2. For one bit bitfields, use "bool" if possible, rather than other
 *    integer types; this will eliminate compilation warnings
 *    about size incompatibilities when compiling with all warnings on:
 *
 *       bool fooA : 1;
 *       bool fooB : 1;
 *
 * 3. Instead of using explicit numbers for bit field widths such that they
 *    add up to a word size such as 32 or 64, try to express it in terms
 *    of SIZE_BITS or _INT_BITS_, etc.  E.g., in ReadPaths.h:
 *        _SSZT_ index:(SIZE_BITS-1);
 *        _SSZT_ set:1;
 *
 * 4. If there are constants that can be used to define bitfield
 *    sizes across multiple classes, place them here.
 ***************************************************************************/



#endif // _COMPATIBILITY_H_
