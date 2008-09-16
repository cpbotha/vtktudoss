// -*- C++ -*-

#if !defined(__stochastic_TauLeaping_ipp__)
#error This file is an implementation detail of TauLeaping.
#endif

BEGIN_NAMESPACE_STOCHASTIC

// Constructor.
template<class _State, bool _MinimumAllowedStep>
inline
TauLeaping<_State, _MinimumAllowedStep>::
TauLeaping(const State& state, const ReactionSet& reactionSet) :
  // Copy.
  _state(state),
  _reactionSet(reactionSet),
  _propensities(state.getNumberOfReactions()),
  // Construct.
  _discreteUniformGenerator(),
  _normalGenerator(&_discreteUniformGenerator),
  // CONTINUE: normal threshhold.
  _poissonGenerator(&_normalGenerator, 1000),
  _mu(state.getNumberOfSpecies()),
  _sigmaSquared(state.getNumberOfSpecies()),
  _highestOrder(state.getNumberOfSpecies()),
  _highestIndividualOrder(state.getNumberOfSpecies()),
  _stepCount(0) {
}

template<class _State, bool _MinimumAllowedStep>
inline
void
TauLeaping<_State, _MinimumAllowedStep>::
initialize(const typename State::PopulationsContainer& populations,
	   const Number time) {
  // Initialize the state.
  _state.setPopulations(populations);
  _state.setTime(time);
  _state.resetReactionCounts();

  //
  // Compute the orders for the species.
  //
  typedef typename ReactionSet::Reaction Reaction;
  // Initialize the arrays.
  std::fill(_highestOrder.begin(), _highestOrder.end(), 0);
  std::fill(_highestIndividualOrder.begin(), _highestIndividualOrder.end(), 0);

  int sum, index, order;
  // Loop over the reactions.
  for (std::size_t n = 0; n != _state.getNumberOfReactions(); ++n) {
    const Reaction& reaction = _reactionSet.getReaction(n);
    // The sum of the reactant coefficients.
    sum = ads::computeSum(reaction.getReactants());
    // Loop over the reactant species.
    for (int i = 0; i != reaction.getReactants().size(); ++i) {
      index = reaction.getReactants().getIndex(i);
      order = reaction.getReactants()[i];
      if (sum > _highestOrder[index]) {
	_highestOrder[index] = sum;
      }
      if (order > _highestIndividualOrder[index]) {
	_highestIndividualOrder[index] = order;
      }
    }
  }

  _stepCount = 0;
}

// Try to take a step.  Return true if a step is taken.
template<class _State, bool _MinimumAllowedStep>
template<typename _TerminationCondition>
inline
bool
TauLeaping<_State, _MinimumAllowedStep>::
step(MemberFunctionPointer method, const Number epsilon,
     _TerminationCondition terminationCondition) {
  // If we have reached the termination condition.
  if (terminationCondition(_state)) {
    return false;
  }

  // Compute the propensities. We will use these in computing tau and in 
  // firing reactions.
  computePropensities();

  // Compute the time leap.
  Number tau = computeTau(epsilon);
  // Advance the time by tau.
  Number time = _state.getTime() + tau;
  // If the time leap will take us past the end time.
  if (time > terminationCondition.getEndTime()) {
    tau = terminationCondition.getEndTime() - _state.getTime();
    // Advance the time to the ending time.
    time = terminationCondition.getEndTime();
  }

  // If the state is not valid, do not advance the time and return false.
  if (! (this->*method)(tau)) {
    return false;
  }

  // Advance the time.
  _state.setTime(time);
  return true;
}

// Try to take a step.  Return true if a step is taken.
template<class _State, bool _MinimumAllowedStep>
template<typename _TerminationCondition>
inline
bool
TauLeaping<_State, _MinimumAllowedStep>::
stepFixed(MemberFunctionPointer method, Number tau,
	  _TerminationCondition terminationCondition) {
  // If we have reached the termination condition.
  if (terminationCondition(_state)) {
    return false;
  }

  // Compute the propensities.
  computePropensities();

  // Advance the time by tau.
  Number time = _state.getTime() + tau;
  // If the time leap will take us past the end time.
  if (time > terminationCondition.getEndTime()) {
    tau = terminationCondition.getEndTime() - _state.getTime();
    // Advance the time to the ending time.
    time = terminationCondition.getEndTime();
  }

  // Advance the state.
  // If the state is not valid, do not advance the time and return false.
  if (! (this->*method)(tau)) {
    return false;
  }

  // Advance the time.
  _state.setTime(time);
  return true;
}


template<class _State, bool _MinimumAllowedStep>
inline
bool
TauLeaping<_State, _MinimumAllowedStep>::
stepForward(const Number tau) {
  ++_stepCount;
  // Advance the state.
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, _poissonGenerator(_propensities[m] * tau));
  }
  // Return true if the state is valid.
  return _state.isValid();
}

template<class _State, bool _MinimumAllowedStep>
inline
bool
TauLeaping<_State, _MinimumAllowedStep>::
stepMidpoint(const Number tau) {
  ++_stepCount;
  // Now the propensities have been calculated at the beginning of the 
  // time interval.
#ifdef DEBUG_stlib
  assert(_state.getNumberOfSpecies() == _p.size());
#endif

  // Determine the midpoint populations.
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  const Number half = 0.5 * tau;
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
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
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, _poissonGenerator(_propensities[m] * tau));
  }
  // Return true if the state is valid.
  return _state.isValid();
}

template<class _State, bool _MinimumAllowedStep>
inline
bool
TauLeaping<_State, _MinimumAllowedStep>::
stepRungeKutta4(const Number tau) {
  ++_stepCount;
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
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, 0.5 * _k1[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k2.size(); ++i) {
    _k2[i] = tau * _propensities[i];
  }

  // k3
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, 0.5 * _k2[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k3.size(); ++i) {
    _k3[i] = tau * _propensities[i];
  }

  // k4
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
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
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, _poissonGenerator(_propensities[m]));
  }
  // Return true if the state is valid.
  return _state.isValid();
}

template<class _State, bool _MinimumAllowedStep>
inline
void
TauLeaping<_State, _MinimumAllowedStep>::
computePropensities(const PopulationsContainer& populations) {
  for (std::size_t m = 0; m < _propensities.size(); ++m) {
    _propensities[m] = _reactionSet.computePropensity(m, populations);
  }
}

template<class _State, bool _MinimumAllowedStep>
inline
typename TauLeaping<_State, _MinimumAllowedStep>::Number
TauLeaping<_State, _MinimumAllowedStep>::
computeTau(const Number epsilon) {
#ifdef DEBUG_stochastic_TauLeaping
  assert(epsilon > 0);
#endif

  // Compute mu and sigmaSquared.
  computeMuAndSigmaSquared();

  // Initialize tau to infinity.
  Number tau = std::numeric_limits<Number>::max();

  // Take the minimum over the species.
  Number numerator, temp, a, b;
  for (std::size_t n = 0; n != _mu.size(); ++n) {
    // If the n_th species is not a reactant in any reaction.
    if (_highestOrder[n] == 0) {
      // This species does not affect tau.
      continue;
    }
    if (_MinimumAllowedStep) {
      numerator = std::max(epsilon * _state.getPopulation(n) / computeG(n), 
			   1.0);
    }
    else {
      numerator = epsilon * std::max(_state.getPopulation(n), 
				     PopulationType(1)) / computeG(n);
    }

    if (_mu[n] != 0) {
      a = numerator / std::abs(_mu[n]);
    }
    else {
      a = std::numeric_limits<Number>::max();
    }

    if (_sigmaSquared[n] != 0) {
      b = numerator * numerator / _sigmaSquared[n];
    }
    else {
      b = std::numeric_limits<Number>::max();
    }
    temp = std::min(a, b);
    if (temp < tau) {
      tau = temp;
    }
  }

  return tau;
}


// Compute the g described in "Efficient step size selection for the 
// tau-leaping simulation method".
template<class _State, bool _MinimumAllowedStep>
inline
typename TauLeaping<_State, _MinimumAllowedStep>::Number
TauLeaping<_State, _MinimumAllowedStep>::
computeG(const int speciesIndex) const {
  const int order = _highestOrder[speciesIndex];

  if (order == 1) {
    return 1.0;
  }
  else if (order == 2) {
    if (_highestIndividualOrder[speciesIndex] == 1) {
      return 2.0;
    }
    else if (_highestIndividualOrder[speciesIndex] == 2) {
      return 2.0 + 1.0 / (_state.getPopulation(speciesIndex) - 1);
    }
  }
  else if (order == 3) {
    if (_highestIndividualOrder[speciesIndex] == 1) {
      return 3.0;
    }
    else if (_highestIndividualOrder[speciesIndex] == 2) {
      return 1.5 * (2.0 + 1.0 / (_state.getPopulation(speciesIndex) - 1.0));
    }
    else if (_highestIndividualOrder[speciesIndex] == 3) {
      return 3.0 + 1.0 / (_state.getPopulation(speciesIndex) - 1)
	+ 2.0 / (_state.getPopulation(speciesIndex) - 2);
    }
  }

  // Catch any other cases with an assertion failure.
  assert(false);
  return 0;
}

template<class _State, bool _MinimumAllowedStep>
inline
void
TauLeaping<_State, _MinimumAllowedStep>::
computeMuAndSigmaSquared() {
  typedef typename _State::ScvContainer::value_type StateChangeVector;
  // Initialize.
  std::fill(_mu.begin(), _mu.end(), 0.0);
  std::fill(_sigmaSquared.begin(), _sigmaSquared.end(), 0.0);

  Number propensity, value;
  int index;
  // Loop over the reactions.
  for (std::size_t m = 0; m < _propensities.size(); ++m) {
    propensity = _propensities[m];
#ifdef DEBUG_stlib
    assert(propensity >= 0);
#endif
    const StateChangeVector& stateChange = 
      _state.getStateChangeVector(m);
    // For each species that is modified by this reaction.
    for (int i = 0; i != stateChange.size(); ++i) {
      index = stateChange.getIndex(i);
      value = stateChange[i];
      _mu[index] += value * propensity;
      _sigmaSquared[index] += value * value * propensity;
    }
  }
}

END_NAMESPACE_STOCHASTIC
