// -*- C++ -*-

/*! 
  \file packages/hj/defs.h
  \brief Definitions for the Hamilton-Jacobi package.
*/

#if !defined(__hj_defs_h__)
//! Include guard.
#define __hj_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_hj)
#define DEBUG_hj
#endif

//! Begin the hj namespace.
#define BEGIN_NAMESPACE_HJ namespace hj {
//! End the hj namespace.
#define END_NAMESPACE_HJ }

//! All classes and functions in the Hamilton-Jacobi package are in the hj namespace.
BEGIN_NAMESPACE_HJ
END_NAMESPACE_HJ

#endif
