// -*- C++ -*-

/*! 
  \file array/defs.h
  \brief Definitions for the array package.
*/

#if !defined(__array_defs_h__)
//! Include guard.
#define __array_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_array)
#define DEBUG_array
#endif

//! Begin the array namespace.
#define BEGIN_NAMESPACE_ARRAY namespace array {
//! End the array namespace.
#define END_NAMESPACE_ARRAY }

//! All classes and functions in the array package are defined in the array namespace.
namespace array {}

#endif
