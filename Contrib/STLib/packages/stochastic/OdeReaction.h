// -*- C++ -*-

/*!
  \file stochastic/OdeReaction.h
  \brief ODE integration of each reaction.
*/

#if !defined(__stochastic_OdeReaction_h__)
#define __stochastic_OdeReaction_h__

#include "State.h"
#include "ReactionSet.h"
#include "EssTerminationCondition.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_OdeReaction)
#define DEBUG_stochastic_OdeReaction
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Perform ODE integration of each reaction.
/*!
  \param _State The state of the simulation: reactions, populations, and time.
*/
template<class _State = State<> >
class OdeReaction {
  //
  // Public types.
  //
public:

  //! The state.
  typedef _State State;
  //! The number type.
  typedef typename State::Number Number;
  //! The number type for populations.
  typedef typename State::PopulationType PopulationType;
  //! The populations container.
  typedef typename State::PopulationsContainer PopulationsContainer;
  //! The set of reactions.
  typedef ReactionSet<Number> ReactionSet;

  //
  // Private types.
  //
private:

  typedef void (OdeReaction::*MemberFunctionPointer)(Number);

  //
  // Member data.  
  //
private:

  State _state;
  ReactionSet _reactionSet;
  std::vector<Number> _propensities;

  //! The number of steps.
  std::size_t _stepCount;
  //! Used for Runge-Kutta and midpoint.
  PopulationsContainer _p;
  //! Used for Runge-Kutta Cash-Karp
  PopulationsContainer _error;
  //! Used for Runge-Kutta.
  std::vector<Number> _k1, _k2, _k3, _k4, _k5, _k6;
  //! Used for Runge-Kutta Cash-Karp
  typename State::Base _backup;


  //
  // Not implemented.
  //
private:

  //! Default constructor not implemented.
  OdeReaction();
  //! Copy constructor not implemented.
  OdeReaction(const OdeReaction&);
  //! Assignment operator not implemented.
  OdeReaction&
  operator=(const OdeReaction&);
  
  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct.
  OdeReaction(const State& state, const ReactionSet& reactionSet);

  //! Destruct.
  ~OdeReaction()
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Simulation.
  //@{
public:

  //! Initialize the state with the initial populations and time.
  void
  initialize(const typename State::PopulationsContainer& populations,
	     const Number time = 0);

  //! Setup for a simulation with adaptive step size, Cash-Karp Runge-Kutta.
  void
  setupRungeKuttaCashKarp();

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  bool
  simulateRungeKuttaCashKarp(Number epsilon, 
			     _TerminationCondition terminationCondition);

  //! Setup for a simulation with fixed step, forward.
  void
  setupFixedForward();

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  bool
  simulateFixedForward(Number dt, _TerminationCondition terminationCondition);

  //! Setup for a simulation with fixed step, midpoint.
  void
  setupFixedMidpoint();

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  bool
  simulateFixedMidpoint(Number dt, _TerminationCondition terminationCondition);

  //! Setup for a simulation with fixed step, fourth order Runge-Kutta.
  void
  setupFixedRungeKutta4();

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  bool
  simulateFixedRungeKutta4(Number dt,
			   _TerminationCondition terminationCondition);

  //! Setup for a simulation with fixed step, Cash-Karp Runge-Kutta.
  void
  setupFixedRungeKuttaCashKarp();

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  bool
  simulateFixedRungeKuttaCashKarp(Number dt,
				  _TerminationCondition terminationCondition);

  //! Try to take a step.  Return false if the termination condition has been reached.
  template<typename _TerminationCondition>
  bool
  step(MemberFunctionPointer method, const Number epsilon, Number endTime);

  //! Take a step. Return true if the state is valid.
  bool
  stepFixed(MemberFunctionPointer method, Number dt, Number endTime);

private:

  void
  stepForward(Number dt);

  void
  stepMidpoint(Number dt);

  void
  stepRungeKutta4(Number dt);

  //! Perform computations that enable computing a step or the error.
  void
  computeRungeKuttaCashKarp(Number dt);

  //! Update the solution in the step.
  void
  solutionRungeKuttaCashKarp(Number dt);

  //! Compute the relative error if we take the step.
  /*! The relative error is the maximum over the specise of 
    error[n] / max(1, x[n]) where x[n] is the n_th species population and 
    error[n] is the truncation error.
  */
  Number
  errorRungeKuttaCashKarp(Number dt);

  //! Take a step with the specified time step.
  void
  stepRungeKuttaCashKarp(const Number dt) {
    computeRungeKuttaCashKarp(dt);
    solutionRungeKuttaCashKarp(dt);
  }

  //! Try to take a step with the specified time step and error tolerance.
  /*! If the step is too large, reduce until an acceptable step is found.
    Return true if a step can be taken. 
  */
  bool
  stepRungeKuttaCashKarp(Number* dt, Number* nextDt, Number epsilon);

  void
  computePropensities() {
    computePropensities(_state.getPopulations());
  }

  void
  computePropensities(const PopulationsContainer& populations);

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Return a const reference to the state.
  const State&
  getState() const {
    return _state;
  }

  //! Return the number of steps taken.
  std::size_t
  getStepCount() const {
    return _stepCount;
  }

  //@}
};

END_NAMESPACE_STOCHASTIC

#define __stochastic_OdeReaction_ipp__
#include "OdeReaction.ipp"
#undef __stochastic_OdeReaction_ipp__

#endif
