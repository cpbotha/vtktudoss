// -*- C++ -*-

#if !defined(__stochastic_Direct_ipp__)
#error This file is an implementation detail of Direct.
#endif

BEGIN_NAMESPACE_STOCHASTIC

// Try to take a step with the direct SSA method.
// Return true if a step is taken.
template<class _DiscreteFiniteGenerator,
	 class _ExponentialGenerator,
	 class _PropensitiesFunctor,
	 class _State>
template<typename _TerminationCondition>
inline
bool
Direct<_DiscreteFiniteGenerator,
       _ExponentialGenerator,
       _PropensitiesFunctor,
       _State>::
step(_TerminationCondition terminationCondition) {
  // If we have reached the termination condition.
  if (_tau == std::numeric_limits<Number>::max() ||
      terminationCondition(_state, _tau)) {
    return false;
  }

  // Determine the reaction to fire.
  const int reactionIndex = _discreteFiniteGenerator();
#ifdef DEBUG_stochastic_Direct
  assert(_discreteFiniteGenerator.getPmf(reactionIndex) > 0);
#endif

  // Fire the reaction.
  _state.advanceTime(_tau);
  _state.fireReaction(reactionIndex);

  // Recompute the propensities and update the discrete, finite generator.
  updatePropensities(reactionIndex);

  // Compute the next time step.
  _tau = computeTau();

  return true;
}

END_NAMESPACE_STOCHASTIC
