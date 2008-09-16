// -*- C++ -*-

/*! 
  \file concurrent/defs.h
  \brief Definitions for the concurrent algorithms package.
*/

#if !defined(__concurrent_defs_h__)
//! Include guard.
#define __concurrent_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_concurrent)
#define DEBUG_concurrent
#endif

//! Start the concurrent namespace.
#define BEGIN_NAMESPACE_CONCURRENT namespace concurrent {
//! End the concurrent namespace.
#define END_NAMESPACE_CONCURRENT }

//! All classes and functions in the Concurrent Algorithms package are defined in the concurrent namespace.
BEGIN_NAMESPACE_CONCURRENT
END_NAMESPACE_CONCURRENT

#endif
