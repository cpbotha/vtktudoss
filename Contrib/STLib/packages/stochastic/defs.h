// -*- C++ -*-

/*! 
  \file stochastic/defs.h
  \brief Definitions for the computational stochasticetry package.
*/

#if !defined(__stochastic_defs_h__)
//! Include guard.
#define __stochastic_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_stochastic)
#define DEBUG_stochastic
#endif

//! Begin the stochastic namespace.
#define BEGIN_NAMESPACE_STOCHASTIC namespace stochastic {
//! End the stochastic namespace.
#define END_NAMESPACE_STOCHASTIC }

//! All classes and functions in the stochastic package are defined in the stochastic namespace.
namespace stochastic {}

#endif
