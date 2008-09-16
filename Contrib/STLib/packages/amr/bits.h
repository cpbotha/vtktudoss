// -*- C++ -*-

/*! 
  \file amr/integer/bits.h
  \brief Integer information and bit operations.
*/

#if !defined(__amr_bits_h__)
#define __amr_bits_h__

#include "defs.h"

#include <boost/static_assert.hpp>

#include <tr1/array>

#include <limits>

#include <cassert>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_bits)
#define DEBUG_amr_bits
#endif

BEGIN_NAMESPACE_AMR

namespace {

  template<int enumeration>
  struct
  IntegerTypes;

  template<>
  struct
  IntegerTypes<1> {
    typedef long Type;
  };

  template<>
  struct
  IntegerTypes<2> {
    typedef int Type;
  };

  template<>
  struct
  IntegerTypes<3> {
    typedef short Type;
  };

  template<>
  struct
  IntegerTypes<4> {
    typedef char Type;
  };

}

//! Select an integer type with at least the specified number of bits.
/*!
  The Boost library has a more sophisticated solution, but this is sufficient
  for my purposes.
*/
template<int Bits>
struct
Integer {
  typedef typename IntegerTypes<
    (Bits - 1 <= std::numeric_limits<long>::digits) +
      (Bits - 1 <= std::numeric_limits<int>::digits) +
      (Bits - 1 <= std::numeric_limits<short>::digits) +
      (Bits - 1 <= std::numeric_limits<char>::digits)>::Type Type;
};



namespace {

  template<int enumeration>
  struct
  UnsignedIntegerTypes;

  template<>
  struct
  UnsignedIntegerTypes<1> {
    typedef unsigned long Type;
  };

  template<>
  struct
  UnsignedIntegerTypes<2> {
    typedef unsigned int Type;
  };

  template<>
  struct
  UnsignedIntegerTypes<3> {
    typedef unsigned short Type;
  };

  template<>
  struct
  UnsignedIntegerTypes<4> {
    typedef unsigned char Type;
  };

}

//! Select an unsigned integer type with at least the specified number of bits.
/*!
  The Boost library has a more sophisticated solution, but this is sufficient
  for my purposes.
*/
template<int Bits>
struct
UnsignedInteger {
  typedef typename UnsignedIntegerTypes<
    (Bits - 1 <= std::numeric_limits<long>::digits) +
      (Bits - 1 <= std::numeric_limits<int>::digits) +
      (Bits - 1 <= std::numeric_limits<short>::digits) +
      (Bits - 1 <= std::numeric_limits<char>::digits)>::Type Type;
};


//! Extract bits, starting with the least significant.
template<typename _Integer, std::size_t _Size>
inline
void
getBits(_Integer n, std::tr1::array<bool, _Size>* bits) {
  BOOST_STATIC_ASSERT(_Size > 0);
  (*bits)[0] = n % 2;
  for (std::size_t i = 1; i != _Size; ++i) {
    n /= 2;
    (*bits)[i] = n % 2;
  }
}

//! Reverse the bits of an unsigned integer type.
template<typename _Integer>
inline
_Integer
reverseBits(_Integer source) {
  // Get the least significant bit.
  _Integer reversed = source;
  // The number of shifts we will make.
  int shift = std::numeric_limits<_Integer>::digits - 1;
  // Loop while there are non-zero bits left in the source.
  for (source >>= 1; source; source >>= 1) {
    reversed <<= 1;
    reversed |= source & 1;
    --shift;
  }
  // Do a shift when some of the source's most significant bits are zero.
  reversed <<= shift;

  return reversed;
}

//! Reverse the n least significant bits of an unsigned integer type.
/*!
  The more significant bits will be zero.
*/
template<typename _Integer>
inline
_Integer
reverseBits(_Integer source, int n) {
#ifdef DEBUG_amr_bits
  assert(n >= 0);
#endif
  _Integer reversed = 0;
  for (; n; --n) {
    reversed <<= 1;
    reversed |= source & 1;
    source >>= 1;
  }
  return reversed;
}

//! Interlace the n least significant bits.
/*!
  The return type must be specified explicitly.  It cannot be deduced from
  the arguments.
*/
template<typename _ResultInteger, typename _ArgumentInteger, std::size_t _N>
inline
_ResultInteger
interlaceBits(std::tr1::array<_ArgumentInteger, _N> sources, const int n) {
#ifdef DEBUG_amr_bits
  assert(n >= 0);
#endif
  _ResultInteger reversed = 0;
  // Interlace in reverse order.
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != _N; ++j) {
      reversed <<= 1;
      reversed |= sources[j] & 1;
      sources[j] >>= 1;
    }
  }
  // Reverse the bits.
  return reverseBits(reversed, n * _N);
}

//! Unlace the n least significant bits in each coordinate.
/*!
  The return type must be specified explicitly.  It cannot be deduced from
  the arguments.
*/
template<typename _SourceInteger, typename _TargetInteger, std::size_t _N>
inline
void
unlaceBits(_SourceInteger source, const std::size_t n, 
	   std::tr1::array<_TargetInteger, _N>* targets) {
#ifdef DEBUG_amr_bits
  assert(n >= 0);
  assert(std::size_t(std::numeric_limits<_SourceInteger>::digits) >= n * _N);
#endif

  // Clear the target coordinates.
  // CONTINUE: assign is broken in GCC 4.0.
  //targets->assign(_TargetInteger(0));
  std::fill(targets->begin(), targets->end(), _TargetInteger(0));
  // Unlace into reverse order targets.
  for (std::size_t i = 0; i != n; ++i) {
    for (std::size_t j = 0; j != _N; ++j) {
      (*targets)[j] <<= 1;
      (*targets)[j] |= source & 1;
      source >>= 1;
    }
  }
  // Reverse the bits.
  for (int j = 0; j != _N; ++j) {
    (*targets)[j] = reverseBits((*targets)[j], n);
  }
}

END_NAMESPACE_AMR

#endif
