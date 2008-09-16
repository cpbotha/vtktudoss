// -*- C++ -*-

/*! 
  \file stochastic/reactionPropensityInfluence.h
  \brief Function to compute the reaction-propensity influence.
*/

#if !defined(__stochastic_reactionPropensityInfluence_h__)
#define __stochastic_reactionPropensityInfluence_h__

#include "Reaction.h"

#include "../ads/array/StaticArrayOfArrays.h"

#include <set>
#include <map>

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_reactionPropensityInfluence)
#define DEBUG_stochastic_reactionPropensityInfluence
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Compute the reaction-propensity influence.
/*!
  \param numberOfSpecies The number of species.
  \param reactionsBeginning The beginning of the reactions.
  \param reactionsEnd The end of the reactions.
  \param influence For each reaction, this will list the reaction 
  propensities that may be affected if the reaction fires.
  \param includeSelf Whether a reaction should include itself in the set
  of influenced reactions.
*/
template<typename ForwardIterator>
void
computeReactionPropensityInfluence(int numberOfSpecies,
				   ForwardIterator reactionsBeginning,
				   ForwardIterator reactionsEnd,
				   ads::StaticArrayOfArrays<int>* influence,
				   bool includeSelf);

END_NAMESPACE_STOCHASTIC

#define __stochastic_reactionPropensityInfluence_ipp__
#include "reactionPropensityInfluence.ipp"
#undef __stochastic_reactionPropensityInfluence_ipp__

#endif
