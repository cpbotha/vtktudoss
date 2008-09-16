// -*- C++ -*-

/*! 
  \file elc/defs.h
  \brief Definitions for the elc package.
*/

#if !defined(__elc_defs_h__)
//! Include guard.
#define __elc_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_elc)
#define DEBUG_elc
#endif

//! Begin the elc namespace.
#define BEGIN_NAMESPACE_ELC namespace elc {
//! End the elc namespace.
#define END_NAMESPACE_ELC }

//! All classes in the Eulerian-Lagrangian Coupling package are defined in the elc namespace.
BEGIN_NAMESPACE_ELC
END_NAMESPACE_ELC

#endif
