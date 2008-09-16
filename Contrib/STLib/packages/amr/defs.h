// -*- C++ -*-

/*! 
  \file amr/defs.h
  \brief Definitions for the adaptive mesh refinement package.
*/

#if !defined(__amr_defs_h__)
//! Include guard.
#define __amr_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_amr)
#define DEBUG_amr
#endif

//! Begin the amr namespace.
#define BEGIN_NAMESPACE_AMR namespace amr {
//! End the amr namespace.
#define END_NAMESPACE_AMR }

//! All classes and functions in the amr package are defined in the amr namespace.
namespace amr {}

#endif
