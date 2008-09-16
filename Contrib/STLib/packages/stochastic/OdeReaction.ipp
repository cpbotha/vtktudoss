// -*- C++ -*-

#if !defined(__stochastic_OdeReaction_ipp__)
#error This file is an implementation detail of OdeReaction.
#endif

BEGIN_NAMESPACE_STOCHASTIC

// Constructor.
template<class _State>
inline
OdeReaction<_State>::
OdeReaction(const State& state, const ReactionSet& reactionSet) :
  // Copy.
  _state(state),
  _reactionSet(reactionSet),
  _propensities(state.getNumberOfReactions()),
  _stepCount(0) {
}

template<class _State>
inline
void
OdeReaction<_State>::
initialize(const typename State::PopulationsContainer& populations,
	   const Number time) {
  // Initialize the state.
  _state.setPopulations(populations);
  _state.setTime(time);
  _state.resetReactionCounts();
  _stepCount = 0;
}

template<class _State>
inline
void
OdeReaction<_State>::
setupRungeKuttaCashKarp() {
  setupFixedRungeKuttaCashKarp();
  _error.resize(_state.getNumberOfSpecies());
}

// Simulate until the termination condition is reached.
template<class _State>
template<typename _TerminationCondition>
inline
bool
OdeReaction<_State>::
simulateRungeKuttaCashKarp(const Number epsilon, 
			   _TerminationCondition terminationCondition) {
  const Number endTime = terminationCondition.getEndTime();
  // Rough guess for an initial dt.
  Number dt = 0.01 * (endTime - _state.getTime());
  Number nextDt;

  // Step until we reach the termination condition.
  while (! terminationCondition(_state)) {
    const Number initialDt = dt;
    bool finish = false;
    // Reduce dt if it will take us past the end time.
    if (_state.getTime() + dt > endTime) {
      dt = endTime - _state.getTime();
      finish = true;
    }
    // Try a step.
    if (! stepRungeKuttaCashKarp(&dt, &nextDt, epsilon)) {
      // If we encounter an error return false.
      return false;
    }
    // If this is the last step.
    if (finish && dt == initialDt) {
      // Do this to avoid problems with round-off errors.
      _state.setTime(endTime);
    }
  }
  return true;
}

template<class _State>
inline
void
OdeReaction<_State>::
setupFixedForward() {
}

// Simulate with fixed size steps until the termination condition is reached.
template<class _State>
template<typename _TerminationCondition>
inline
bool
OdeReaction<_State>::
simulateFixedForward(const Number dt,
		     _TerminationCondition terminationCondition) {
  // Step until we reach the termination condition.
  while (! terminationCondition(_state)) {
    // Try a step.
    if (! stepFixed(&OdeReaction::stepForward, dt, 
		    terminationCondition.getEndTime())) {
      // If we encounter an error return false.
      return false;
    }
  }
  return true;
}

template<class _State>
inline
void
OdeReaction<_State>::
setupFixedMidpoint() {
  _p.resize(_state.getNumberOfSpecies());
}

// Simulate with fixed size steps until the termination condition is reached.
template<class _State>
template<typename _TerminationCondition>
inline
bool
OdeReaction<_State>::
simulateFixedMidpoint(const Number dt, 
		      _TerminationCondition terminationCondition) {
  // Step until we reach the termination condition.
  while (! terminationCondition(_state)) {
    // Try a step.
    if (! stepFixed(&OdeReaction::stepMidpoint, dt, 
		    terminationCondition.getEndTime())) {
      // If we encounter an error return false.
      return false;
    }
  }
  return true;
}

template<class _State>
inline
void
OdeReaction<_State>::
setupFixedRungeKutta4() {
  _p.resize(_state.getNumberOfSpecies());
  _k1.resize(_state.getNumberOfReactions());
  _k2.resize(_state.getNumberOfReactions());
  _k3.resize(_state.getNumberOfReactions());
  _k4.resize(_state.getNumberOfReactions());
}

//! Simulate with fixed size steps until the termination condition is reached.
template<class _State>
template<typename _TerminationCondition>
inline
bool
OdeReaction<_State>::
simulateFixedRungeKutta4(const Number dt, 
			 _TerminationCondition terminationCondition) {
  // Step until we reach the termination condition.
  while (! terminationCondition(_state)) {
    // Try a step.
    if (! stepFixed(&OdeReaction::stepRungeKutta4, dt, 
		    terminationCondition.getEndTime())) {
      // If we encounter an error return false.
      return false;
    }
  }
  return true;
}

template<class _State>
inline
void
OdeReaction<_State>::
setupFixedRungeKuttaCashKarp() {
  _p.resize(_state.getNumberOfSpecies());
  _k1.resize(_state.getNumberOfReactions());
  _k2.resize(_state.getNumberOfReactions());
  _k3.resize(_state.getNumberOfReactions());
  _k4.resize(_state.getNumberOfReactions());
  _k5.resize(_state.getNumberOfReactions());
  _k6.resize(_state.getNumberOfReactions());
}

//! Simulate with fixed size steps until the termination condition is reached.
template<class _State>
template<typename _TerminationCondition>
inline
bool
OdeReaction<_State>::
simulateFixedRungeKuttaCashKarp(const Number dt, 
				_TerminationCondition terminationCondition) {
  // Step until we reach the termination condition.
  while (! terminationCondition(_state)) {
    // Try a step.
    if (! stepFixed(&OdeReaction::stepRungeKuttaCashKarp, dt, 
		    terminationCondition.getEndTime())) {
      // If we encounter an error return false.
      return false;
    }
  }
  return true;
}

#if 0
// Take a step. Return true if the state is valid.
template<class _State>
inline
bool
OdeReaction<_State>::
step(MemberFunctionPointer method, Number dt, const Number endTime) {
  CONTNUE;
  // Decrease dt if it will take us past the end time.
  Number time = _state.getTime() + dt;
  // If the time leap will take us past the end time.
  if (time > endTime) {
    dt = endTime - _state.getTime();
    // Advance the time to the ending time.
    time = endTime;
  }

  // Advance the state.
  (this->*method)(dt);
  // Advance the time.
  _state.setTime(time);
  // Return true if there are no negative populations.
  return _state.isValid();
}
#endif

// Take a step. Return true if the state is valid.
template<class _State>
inline
bool
OdeReaction<_State>::
stepFixed(MemberFunctionPointer method, Number dt, const Number endTime) {
  // Advance the time by dt.
  Number time = _state.getTime() + dt;
  // If the time leap will take us past the end time.
  if (time > endTime) {
    dt = endTime - _state.getTime();
    // Advance the time to the ending time.
    time = endTime;
  }

  // Advance the state.
  (this->*method)(dt);
  // Advance the time.
  _state.setTime(time);
  // Return true if there are no negative populations.
  return _state.isValid();
}


template<class _State>
inline
void
OdeReaction<_State>::
stepForward(const Number dt) {
  ++_stepCount;
  // Compute the propensities.
  computePropensities();
  // Advance the state.
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, _propensities[m] * dt);
  }
}

template<class _State>
inline
void
OdeReaction<_State>::
stepMidpoint(const Number dt) {
#ifdef DEBUG_stlib
  assert(_state.getNumberOfSpecies() == _p.size());
#endif
  ++_stepCount;

  // Compute the propensities.
  computePropensities();
  // Determine the midpoint populations.
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  const Number half = 0.5 * dt;
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, _propensities[m] * half);
  }

  // Determine the midpoint propensities.
  computePropensities(_p);
  // Take a step with the midpoint propensities.  
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, _propensities[m] * dt);
  }
}

template<class _State>
inline
void
OdeReaction<_State>::
stepRungeKutta4(const Number dt) {
#ifdef DEBUG_stlib
  assert(_state.getNumberOfSpecies() == _p.size());
  assert(_propensities.size() == _k1.size());
#endif
  ++_stepCount;
  // Compute the propensities.
  computePropensities();

  // k1
  for (std::size_t i = 0; i != _k1.size(); ++i) {
    _k1[i] = dt * _propensities[i];
  }

  // k2
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, 0.5 * _k1[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k2.size(); ++i) {
    _k2[i] = dt * _propensities[i];
  }

  // k3
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, 0.5 * _k2[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k3.size(); ++i) {
    _k3[i] = dt * _propensities[i];
  }

  // k4
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, _k3[m]);
  }
  computePropensities(_p);
  for (std::size_t i = 0; i != _k4.size(); ++i) {
    _k4[i] = dt * _propensities[i];
  }

  // Take a step with the average propensities.  
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, (1. / 6.) * (_k1[m] + 
					2 * (_k2[m] + _k3[m]) + _k4[m]));
  }
}

template<class _State>
inline
void
OdeReaction<_State>::
computeRungeKuttaCashKarp(const Number dt) {
#ifdef DEBUG_stlib
  assert(_state.getNumberOfSpecies() == _p.size() &&
	 _propensities.size() == _k1.size() &&
	 _propensities.size() == _k2.size() &&
	 _propensities.size() == _k3.size() &&
	 _propensities.size() == _k4.size() &&
	 _propensities.size() == _k5.size() &&
	 _propensities.size() == _k6.size())
#endif
  const Number 
    // Not used because there is no explicit time dependence.
    //a2 = 0.2,
    //a3 = 0.3,
    //a4 = 0.6,
    //a5 = 1.0,
    //a6 = 0.875,
    b21 = 0.2,
    b31 = 3./40.,
    b32 = 9./40.,
    b41 = 0.3,
    b42 = -0.9,
    b43 = 1.2,
    b51 = -11./54.,
    b52 = 2.5,
    b53 = -70./27.,
    b54 = 35./27.,
    b61 = 1631./55296.,
    b62 = 175./512.,
    b63 = 575./13824.,
    b64 = 44275./110592.,
    b65 = 253./4096.;
  
  ++_stepCount;

  // Initial propensities.
  computePropensities();
  std::copy(_propensities.begin(), _propensities.end(), _k1.begin());

  // First step.
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, dt * b21 * _k1[m]);
  }

  // Second step.
  computePropensities(_p);
  std::copy(_propensities.begin(), _propensities.end(), _k2.begin());
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, dt * (b31 * _k1[m] + b32 * _k2[m]));
  }

  // Third step.
  computePropensities(_p);
  std::copy(_propensities.begin(), _propensities.end(), _k3.begin());
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, dt * (b41 * _k1[m] + b42 * _k2[m] +
				      b43 * _k3[m]));
  }

  // Fourth step.
  computePropensities(_p);
  std::copy(_propensities.begin(), _propensities.end(), _k4.begin());
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, dt * (b51 * _k1[m] + b52 * _k2[m] +
				      b53 * _k3[m] + b54 * _k4[m]));
  }

  // Fifth step.
  computePropensities(_p);
  std::copy(_propensities.begin(), _propensities.end(), _k5.begin());
  std::copy(_state.getPopulations().begin(), _state.getPopulations().end(), 
	    _p.begin());
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_p, m, dt * (b61 * _k1[m] + b62 * _k2[m] +
				      b63 * _k3[m] + b64 * _k4[m] +
				      b65 * _k5[m]));
  }

  // Sixth step.
  computePropensities(_p);
  std::copy(_propensities.begin(), _propensities.end(), _k6.begin());
}

template<class _State>
inline
void
OdeReaction<_State>::
solutionRungeKuttaCashKarp(const Number dt) {
  const Number
    c1 = 37./378.,
    c3 = 250./621.,
    c4 = 125./594.,
    c6 = 512./1771.;

  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(m, dt * (c1 * _k1[m] + c3 * _k3[m] + c4 * _k4[m] +
				 c6 * _k6[m]));
  }
}

template<class _State>
inline
typename OdeReaction<_State>::Number
OdeReaction<_State>::
errorRungeKuttaCashKarp(const Number dt) {
  const Number
    c1 = 37./378.,
    c3 = 250./621.,
    c4 = 125./594.,
    c6 = 512./1771.,
    dc1 = c1 - 2825./27648.,
    dc3 = c3 - 18575./48384.,
    dc4 = c4 - 13525./55296.,
    dc5 = -277./14336.,
    dc6 = c6 - 0.25;

  std::fill(_error.begin(), _error.end(), 0.);
  for (std::size_t m = 0; m != _state.getNumberOfReactions(); ++m) {
    _state.fireReaction(&_error, m, dt * (dc1 * _k1[m] + dc3 * _k3[m] +
					  dc4 * _k4[m] + dc5 * _k5[m] +
					  dc6 * _k6[m]));
  }
#if 0
  // Here is a straight-forward method for determining the maximum relative
  // error.
  Number error = 0;
  for (std::size_t n = 0; n != _error.size(); ++n) {
    error = std::max(error, std::abs(_error[n]) / 
		     std::max(1., _state.getPopulation(n)));
  }
  return error
#endif
  // Here is a more efficient method that avoids costly divisions.
  // e[i] / s[i] < e[j] / s[j];
  // e[i] * s[j] < e[j] * s[i];
  Number error = std::abs(_error[0]);
  Number scale = std::max(1., _state.getPopulation(0));
  for (std::size_t n = 1; n != _error.size(); ++n) {
    const Number e = std::abs(_error[n]);
    const Number s = std::max(1., _state.getPopulation(n));
    if (error * s < e * scale) {
      error = e;
      scale = s;
    }
  }
  return error / scale;
}

template<class _State>
inline
bool
OdeReaction<_State>::
stepRungeKuttaCashKarp(Number* dt, Number* nextDt, const Number epsilon) {
  const Number 
    Safety = 0.9,
    PGrow = -0.2,
    PShrink = -0.25,
    // (5 / Safety)^(1 / PGrow)
    ErrorCondition = 1.89e-4;

  Number scaledError;
  while (true) {
    // Determine the maximum relative error with the current step size.
    computeRungeKuttaCashKarp(*dt);
    scaledError = errorRungeKuttaCashKarp(*dt) / epsilon;
    // If the error is acceptable.
    if (scaledError <= 1.) {
      // Store the populations and reaction counts.
      _backup = _state;
      // Take the step.
      solutionRungeKuttaCashKarp(*dt);
      // If there are no negative populations.
      if (_state.isValid()) {
	// Accept the step.
	_state.advanceTime(*dt);
	break;
      }
      else {
	// Roll back.
	_state = _backup;
	// Reduce the step size.
	*dt *= 0.5;
      }
    }
    else {
      // Reduce the step size.
      const Number candidateDt = 
	Safety * (*dt) * std::pow(scaledError, PShrink);
      // Note: dt could be negative.
      *dt = (*dt >= 0. ? std::max(candidateDt, 0.1 * (*dt)) : 
	     std::min(candidateDt, 0.1 * (*dt)));
    }
    // Check for step size underflow.
    if (_state.getTime() + *dt == _state.getTime()) {
      return false;
    }
  }
  // Compute the next step size. Allow no more than a factor of 5 increase.
  if (scaledError > ErrorCondition) {
    *nextDt = Safety * (*dt) * std::pow(scaledError, PGrow);
  }
  else {
    *nextDt = 5. * (*dt);
  }
  return true;
}

template<class _State>
inline
void
OdeReaction<_State>::
computePropensities(const PopulationsContainer& populations) {
  for (std::size_t m = 0; m < _propensities.size(); ++m) {
    _propensities[m] = _reactionSet.computePropensity(m, populations);
  }
}

END_NAMESPACE_STOCHASTIC
