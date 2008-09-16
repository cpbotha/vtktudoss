// -*- C++ -*-

#if !defined(__stochastic_HybridDirectTauLeaping_ipp__)
#error This file is an implementation detail of HybridDirectTauLeaping.
#endif

BEGIN_NAMESPACE_STOCHASTIC

// Initialize the state with the initial populations and time.
template<class _State>
inline
void
HybridDirectTauLeaping<_State>::
initialize(const typename State::PopulationsContainer& populations,
	   const Number time) {
  // Initialize the state.
  _state.setPopulations(populations);
  _state.setTime(time);
  _state.resetReactionCounts();
  _directStepCount = 0;
  _tauLeapingStepCount = 0;

  // Initially put all of the reactions in the slow group.
  _discreteFiniteGenerator.clear();
  _tauLeaping.clear();
  for (std::size_t i = 0; i != _propensities.size(); ++i) {
    _discreteFiniteGenerator.insert(i);
  }

  // Compute the propensities.
  computePropensities();
  // Compute the initial exponential deviate.
  _exponentialDeviate = _exponentialGenerator();
}

// CONTINUE: Can I avoid computing certain propensities?
// Try to take a step with the direct/tau-leaping method.
// Return true if a step is taken.
template<class _State>
template<typename _TerminationCondition>
inline
bool
HybridDirectTauLeaping<_State>::
step(MemberFunctionPointer method, _TerminationCondition terminationCondition) {
  const Number minimumReactions = 0.1;

  // If we have reached the termination condition.
  if (terminationCondition(_state)) {
    return false;
  }

  // Compute the propensities. We will use these in computing tau and in 
  // firing reactions.
  computePropensities();

  // Compute the time leap.
  const Number originalTau = 
    _tauLeaping.computeStep(_state.getStateChangeVectors(), _propensities,
			    _state.getPopulations());
  Number tau = originalTau;

  // Move the volatile reactions to the direct group.
  moveVolatileAndSlowToDirect(minimumReactions / tau);

  const Number directPmfSum = computeDirectPmfSum();

  bool shouldTakeDirectStep = false;
  // If the tau-leaping step as large as the direct step. (The direct step is
  // _exponentialDeviate / directPmfSum.)
  if (directPmfSum != 0 && (tau == std::numeric_limits<Number>::max() ||
			    directPmfSum * tau >= _exponentialDeviate)) {
    shouldTakeDirectStep = true;
    // The tau-leaping step is the direct step.
    tau = _exponentialDeviate / directPmfSum;
  }

  // Determine the new time.
  Number time = _state.getTime() + tau;
  // If the time leap will take us past the end time.
  if (time > terminationCondition.getEndTime()) {
    tau = terminationCondition.getEndTime() - _state.getTime();
    // Advance the time to the ending time.
    time = terminationCondition.getEndTime();
    shouldTakeDirectStep = false;
  }

  // Advance the state by firing the tau-leaping reactions.
  if (! (this->*method)(tau)) {
    return false;
  }
  
  // Add the propensities contribution for this time step to the discrete
  // generator. 
  _discreteFiniteGenerator.addPmf(_propensities, tau);

  // If a slow reaction fires at the end of the step.
  if (shouldTakeDirectStep) {
    // Determine the reaction to fire.
    const int reactionIndex = _discreteFiniteGenerator();
    // Fire the reaction.
    _state.fireReaction(reactionIndex);
    // If the fired reaction is not volatile and is fast, move it to the 
    // tau-leaping group.
    // CONTINUE
    if (! isVolatile(reactionIndex)) {
      // If there are no reactions in the tau-leaping group and this is the
      // fastest reaction.
      // If the reaction would fire a sufficient number of times per step
      // in the tau-leaping group.
      if ((_tauLeaping.empty() && 
	   _discreteFiniteGenerator.isMaximum(reactionIndex)) || 
	  (_propensities[reactionIndex] / tau) * originalTau > 
	  minimumReactions) {
	_discreteFiniteGenerator.erase(reactionIndex);
	_tauLeaping.insert(reactionIndex);
      }
    }
    // Compute a new exponential deviate.
    _exponentialDeviate = _exponentialGenerator();
    ++_directStepCount;
  }
  else {
    _exponentialDeviate -= tau * directPmfSum;
  }

  // If the state is not valid, do not advance the time and return false.
  if (! _state.isValid()) {
    return false;
  }

  // Advance the time.
  _state.setTime(time);
  return true;
}


template<class _State>
inline
bool
HybridDirectTauLeaping<_State>::
stepForward(const Number tau) {
  const std::vector<int>& tauLeapingReactions = 
    _tauLeaping.getActiveReactions();
  _propensities *= tau;
  if (! tauLeapingReactions.empty()) {
    for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
      const int m = tauLeapingReactions[i];
      _state.fireReaction(m, _tauLeaping.generatePoisson(_propensities[m]));
    }
    ++_tauLeapingStepCount;
  }
  return _state.isValid();
}


template<class _State>
inline
bool
HybridDirectTauLeaping<_State>::
stepMidpoint(const Number tau) {
  const std::vector<int>& tauLeapingReactions = 
    _tauLeaping.getActiveReactions();
  if (tauLeapingReactions.empty()) {
    _propensities *= tau;
    return true;
  }

  ++_tauLeapingStepCount;
  // Now the propensities have been calculated at the beginning of the 
  // time interval.
#ifdef DEBUG_stlib
  assert(_state.getNumberOfSpecies() == _p.size());
#endif

  // Determine the midpoint populations.
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  const Number half = 0.5 * tau;
  for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
    const int m = tauLeapingReactions[i];
    // Note: The reaction counts are not integer.
    _state.fireReaction(&_p, m, _propensities[m] * half);
  }
  // Check the populations.
  for (std::size_t i = 0; i != _p.size(); ++i) {
    if (_p[i] < 0) {
      return false;
    }
  }
  // Determine the midpoint propensities.
  computePropensities(_p);

  // Take a step with the midpoint propensities.  
  _propensities *= tau;
  for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
    const int m = tauLeapingReactions[i];
    _state.fireReaction(m, _tauLeaping.generatePoisson(_propensities[m]));
  }
  // Return true if the state is valid.
  return _state.isValid();
}


template<class _State>
inline
bool
HybridDirectTauLeaping<_State>::
stepRungeKutta4(const Number tau) {
  const std::vector<int>& tauLeapingReactions = 
    _tauLeaping.getActiveReactions();
  if (tauLeapingReactions.empty()) {
    return true;
  }

  ++_tauLeapingStepCount;
  // Now the propensities have been calculated at the beginning of the 
  // time interval.
#ifdef DEBUG_stlib
  assert(_state.getNumberOfSpecies() == _p.size());
  assert(_propensities.size() == _k1.size());
#endif

  // k1
  for (std::size_t i = 0; i != _k1.size(); ++i) {
    _k1[i] = tau * _propensities[i];
  }

  // k2
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
    const int m = tauLeapingReactions[i];
    _state.fireReaction(&_p, m, 0.5 * _k1[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k2.size(); ++i) {
    _k2[i] = tau * _propensities[i];
  }

  // k3
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
    const int m = tauLeapingReactions[i];
    _state.fireReaction(&_p, m, 0.5 * _k2[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k3.size(); ++i) {
    _k3[i] = tau * _propensities[i];
  }

  // k4
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
    const int m = tauLeapingReactions[i];
    _state.fireReaction(&_p, m, _k3[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k4.size(); ++i) {
    _k4[i] = tau * _propensities[i];
  }

  // Average propensities times tau.
  for (std::size_t i = 0; i != _propensities.size(); ++i) {
    _propensities[i] = (1. / 6.) * (_k1[i] + 2 * (_k2[i] + _k3[i]) + _k4[i]);
  }

  // Take a step with the average propensities.  
  for (std::size_t i = 0; i != tauLeapingReactions.size(); ++i) {
    const int m = tauLeapingReactions[i];
    _state.fireReaction(m, _tauLeaping.generatePoisson(_propensities[m]));
  }

  // Return true if the state is valid.
  return _state.isValid();
}


// Move the volatile reactions to the direct group.
template<class _State>
inline
void
HybridDirectTauLeaping<_State>::
moveVolatileAndSlowToDirect(const Number minimumPropensity) {
  int reaction;
  const std::vector<int>& tauLeapingReactions = 
    _tauLeaping.getActiveReactions();
  std::size_t i = 0;
  while (i != tauLeapingReactions.size()) {
    reaction = tauLeapingReactions[i];
    if (isVolatile(reaction) || _propensities[reaction] < minimumPropensity) {
      _tauLeaping.erase(reaction);
      _discreteFiniteGenerator.insert(reaction);
    }
    else {
      ++i;
    }
  }
}

template<class _State>
inline
bool
HybridDirectTauLeaping<_State>::
isVolatile(const int index) {
  const ads::SparseArray<1, int>& reactants = 
    _reactionSet.getReaction(index).getReactants();
  // For each reactant.
  for (int i = 0; i != reactants.size(); ++i) {
    // If the stoichiometry times the population of the species is less small
    // then the reaction is volatile.
    if (reactants[i] * _state.getPopulations()[reactants.getIndex(i)] < 
	_volatileLimit) {
      return true;
    }
  }
  // CONTINUE: Check the modifiers when I implement that.
  return false;
}

END_NAMESPACE_STOCHASTIC
