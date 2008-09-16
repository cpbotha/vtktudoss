// -*- C++ -*-

/*!
  \file array/IndexTypes.h
  \brief Types for array indexing.
*/

#if !defined(__array_IndexTypes_h__)
#define __array_IndexTypes_h__

#include "MultiArrayStorage.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_IndexTypes)
#define DEBUG_array_IndexTypes
#endif

BEGIN_NAMESPACE_ARRAY

//! Types for array indexing.
template<std::size_t _Dimension>
class
IndexTypes {
  //
  // Enumerations.
  //
public:

  //! The number of dimensions.
  enum {Dimension = _Dimension};

  //
  // Public types.
  //
public:

  // Types for STL compliance.

  //! The size type.
  typedef std::size_t size_type;
  //! Pointer difference type.
  typedef std::ptrdiff_t difference_type;

  // Other types.

  //! An array index is a signed integer.
  typedef difference_type Index;
  //! A list of sizes.
  typedef Array<size_type, Dimension> SizeList;
  //! A list of indices.
  typedef Array<Index, Dimension> IndexList;
  //! The storage order.
  typedef MultiArrayStorage<Dimension> Storage;
};

END_NAMESPACE_ARRAY

#endif
