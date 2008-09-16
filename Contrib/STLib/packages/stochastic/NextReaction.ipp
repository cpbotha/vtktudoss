// -*- C++ -*-

#if !defined(__stochastic_NextReaction_ipp__)
#error This file is an implementation detail of NextReaction.
#endif

BEGIN_NAMESPACE_STOCHASTIC

//----------------------------------------------------------------------------
// Implementation details.
//----------------------------------------------------------------------------

// Try to take a step with the next reaction SSA method.
// Return true if a step is taken.
template<class _ReactionPriorityQueue,
	 class _PropensitiesFunctor,
	 class _State>
template<typename _TerminationCondition>
inline
bool
NextReaction<_ReactionPriorityQueue, _PropensitiesFunctor, _State>::
step(_TerminationCondition terminationCondition) {
  // The index of the reaction that will fire next.
  if (_reactionIndex == -1) {
    _reactionIndex = _reactionPriorityQueue->top();
  }
  // The time at which the next reaction occurs.
  const Number reactionTime = _reactionPriorityQueue->get(_reactionIndex);

  // If no more reactions can fire.
  if (reactionTime == _reactionPriorityQueue->infinity()) {
    return false;
  }

  // The time step.
  const Number tau = reactionTime - _state.getTime();

  // If we have reached the termination condition.
  if (terminationCondition(_state, tau)) {
    return false;
  }

  // Fire the reaction.
  _state.setTime(reactionTime);
  _state.fireReaction(_reactionIndex);

  updatePropensitiesAndReactionTimes(_reactionIndex);
  _reactionIndex = -1;

  return true;
}


template<class _ReactionPriorityQueue,
	 class _PropensitiesFunctor,
	 class _State>
inline
void
NextReaction<_ReactionPriorityQueue, _PropensitiesFunctor, _State>::
_updatePropensitiesAndReactionTimes(Loki::Int2Type<false> /*dummy*/,
				    const int reactionIndex) {
  // Recompute the propensities.
  _propensities.swap(_oldPropensities);
  computePropensities();

  // Update the reaction times for the reactions whose propensities
  // were influenced.
  for (typename ads::StaticArrayOfArrays<int>::const_iterator 
	 i = _reactionInfluence.begin(reactionIndex); 
       i != _reactionInfluence.end(reactionIndex); ++i) {
    // Update the reaction time.
    _reactionPriorityQueue->update(*i, _state.getTime(), _oldPropensities[*i],
				   _propensities[*i]);
  }

  // Compute a new reaction time for the fired reaction.
  _reactionPriorityQueue->pushTop(_state.getTime(),
				  _propensities[reactionIndex]);
}


template<class _ReactionPriorityQueue,
	 class _PropensitiesFunctor,
	 class _State>
inline
void
NextReaction<_ReactionPriorityQueue, _PropensitiesFunctor, _State>::
_updatePropensitiesAndReactionTimes(Loki::Int2Type<true> /*dummy*/,
				    const int reactionIndex) {
  // Update the reaction times for the reactions whose propensities
  // were influenced.
  for (typename ads::StaticArrayOfArrays<int>::const_iterator 
	 i = _reactionInfluence.begin(reactionIndex); 
       i != _reactionInfluence.end(reactionIndex); ++i) {
    // Compute the new propensity for the influenced reaction.
    const Number newPropensity = 
      _propensitiesFunctor(*i, _state.getPopulations());
    // Update the reaction time.
    _reactionPriorityQueue->update(*i, _state.getTime(), _propensities[*i],
				   newPropensity);
    // Update the propensity.
    _propensities[*i] = newPropensity;
  }

#if 0
  // This method is typically slower.  It is faster for partition, cost 
  // adaptive, though.  I have no idea why.
  typename ads::StaticArrayOfArrays<int>::const_iterator i,
    begin = _reactionInfluence.begin(reactionIndex),
    end = _reactionInfluence.end(reactionIndex);
  {
    typename std::vector<Number>::iterator old = _oldPropensities.begin();
    for (i = begin; i != end; ++i) {
      // Record the old propesity.
      *old++ = _propensities[*i];
      // Compute the new propensity for the influenced reaction.
      _propensities[*i] = _propensitiesFunctor(*i, _state.getPopulations());
    }
  }
#endif

  // Compute the new propensity for the fired reaction.
  _propensities[reactionIndex] = 
    _propensitiesFunctor(reactionIndex, _state.getPopulations());

#if 0
  {
    typename std::vector<Number>::const_iterator old =
      _oldPropensities.begin();
    for (i = begin; i != end; ++i) {
      // Compute the new propensity for the influenced reaction.
      // Update the reaction time.
      _reactionPriorityQueue->update(*i, _state.getTime(), *old++,
				     _propensities[*i]);
    }
  }
#endif

  // Compute a new reaction time for the fired reaction.
  _reactionPriorityQueue->pushTop(_state.getTime(),
				  _propensities[reactionIndex]);
}

END_NAMESPACE_STOCHASTIC
