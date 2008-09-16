// -*- C++ -*-

/*! 
  \file stochastic/tauLeapingConcurrent.h
  \brief The tau-leaping method for SSA.
*/

#if !defined(__stochastic_tauLeapingConcurrent_h__)
#define __stochastic_tauLeapingConcurrent_h__

#include "tauLeaping.h"

#include <mpi.h>

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_tauLeapingConcurrent)
#define DEBUG_stochastic_tauLeapingConcurrent
#endif

BEGIN_NAMESPACE_STOCHASTIC

/*!
  \page page_stochastic_tauLeapingConcurrent The concurrent tau-leaping method for SSA.

  CONTINUE.
*/

//-----------------------------------------------------------------------------
/*! \defgroup stochastic_tauLeapingConcurrent The concurrent tau-leaping method for the SSA. */
// @{

//! Advance the state to the specified time.
/*!
  \return The number of steps.
*/
template<typename UnaryFunctor, typename T>
int
computeTauLeapingConcurrentSsa(const MPI::Intracomm& intracomm,
			       State<T>* state, T epsilon, T maximumTime, 
			       int seed = 1);

//! Take the specified number of steps.
template<typename UnaryFunctor, typename T>
void
computeTauLeapingConcurrentSsa(const MPI::Intracomm& intracomm,
			       State<T>* state, T epsilon, int numberOfSteps, 
			       int seed = 1);

// @}
  
END_NAMESPACE_STOCHASTIC

#define __stochastic_tauLeapingConcurrent_ipp__
#include "tauLeapingConcurrent.ipp"
#undef __stochastic_tauLeapingConcurrent_ipp__

#endif
