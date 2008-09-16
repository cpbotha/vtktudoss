// -*- C++ -*-

/*! 
  \file mst/defs.h
  \brief Definitions for the MST package.
*/

#if !defined(__mst_defs_h__)
//! Include guard.
#define __mst_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_mst)
#define DEBUG_mst
#endif

//! Begin the mst namespace.
#define BEGIN_NAMESPACE_MST namespace mst {
//! End the mst namespace.
#define END_NAMESPACE_MST }

//! All classes and functions in the MST package are defined in the mst namespace.
BEGIN_NAMESPACE_MST
END_NAMESPACE_MST

#endif
