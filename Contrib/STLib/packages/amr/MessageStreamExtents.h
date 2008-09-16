// -*- C++ -*-

/*!
  \file amr/MessageStreamExtents.h
  \brief Message stream extents.
*/

#if !defined(__amr_MessageStreamExtents_h__)
#define __amr_MessageStreamExtents_h__

#include "defs.h"

#include <cassert>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_MessageStreamExtents)
#define DEBUG_amr_MessageStreamExtents
#endif

BEGIN_NAMESPACE_AMR

//! Return the padding required for the specified number of bytes. If the data is followed by the specified type.
template<typename _Type, typename _NextType = double>
class
MessageStreamPadding {
public:

  static
  int
  get(const int n) {
    return (sizeof(_NextType) - (n * sizeof(_Type) % sizeof(_NextType))) %
      sizeof(_NextType);
  }
};

//! Message stream extents.
/*!
*/
template<typename _T>
class
MessageStreamExtents {
public:

  //! The extent for a single element.
  enum {Extent = ((sizeof(_T) - 1) / sizeof(double) + 1) * sizeof(double)};

  //! The extent for an array of elements.
  static
  int
  getExtent(const int n) {
    return ((n * sizeof(_T) - 1) / sizeof(double) + 1) * sizeof(double);
  }
};

END_NAMESPACE_AMR

#endif
