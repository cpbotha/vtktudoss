// -*- C++ -*-

/*! 
  \file stochastic/HybridDirectTauLeaping.h
  \brief The direct method for SSA.
*/

#if !defined(__stochastic_HybridDirectTauLeaping_h__)
#define __stochastic_HybridDirectTauLeaping_h__

#include "EssTerminationCondition.h"
#include "State.h"
#include "ReactionSet.h"
#include "TauLeapingDynamic.h"

#include "../numerical/random/exponential/ExponentialGeneratorZiggurat.h"
#include "../numerical/random/discreteFinite/DiscreteFiniteGeneratorDynamic.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_HybridDirectTauLeaping)
#define DEBUG_stochastic_HybridDirectTauLeaping
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Perform a stochastic simulation using Gillespie's direct method.
/*!
  \param _State The state of the simulation: reactions, populations, and time.

  
*/
template<class _State = State<> >
class HybridDirectTauLeaping {
  //
  // Public types.
  //
public:

  //! The state.
  typedef _State State;
  //! The number type for populations.
  typedef typename State::PopulationType PopulationType;
  //! The number type.
  typedef typename State::Number Number;
  //! The populations container.
  typedef typename State::PopulationsContainer PopulationsContainer;
  //! The set of reactions.
  typedef ReactionSet<Number> ReactionSet;
  //! The tau-leaping algorithm.
  typedef TauLeapingDynamic<State> TauLeaping;

  //! The exponential generator.
  typedef numerical::ExponentialGeneratorZiggurat<> ExponentialGenerator;
  //! The discrete, uniform generator.
  typedef typename ExponentialGenerator::DiscreteUniformGenerator
  DiscreteUniformGenerator;
  //! The discrete, finite generator.
  typedef numerical::DiscreteFiniteGeneratorDynamic<DiscreteUniformGenerator>
  DiscreteFiniteGenerator;

  //
  // Private types.
  //
private:

  typedef bool (HybridDirectTauLeaping::*MemberFunctionPointer)(Number);

  //
  // Member data.
  //
private:

  // State.
  State _state;
  ReactionSet _reactionSet;
  std::vector<Number> _propensities;

  // Random number generators.
  DiscreteUniformGenerator _discreteUniformGenerator;
  ExponentialGenerator _exponentialGenerator;
  DiscreteFiniteGenerator _discreteFiniteGenerator;

  // Direct method.
  Number _exponentialDeviate;
  std::size_t _directStepCount;
  // Reactions that will exhaust their reactants in few steps must be in the
  // direct group.
  PopulationType _volatileLimit;

  // Tau-leaping.
  TauLeaping _tauLeaping;
  std::size_t _tauLeapingStepCount;
  //! Used for Runge-Kutta 4 and midpoint.
  PopulationsContainer _p;
  //! Used for Runge-Kutta 4.
  std::vector<Number> _k1, _k2, _k3, _k4;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  HybridDirectTauLeaping();
  //! Copy constructor not implemented.
  HybridDirectTauLeaping(const HybridDirectTauLeaping&);
  //! Assignment operator not implemented.
  HybridDirectTauLeaping&
  operator=(const HybridDirectTauLeaping&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct.
  HybridDirectTauLeaping(const State& state, const ReactionSet& reactionSet,
			 const Number epsilon) :
    // State.
    _state(state),
    _reactionSet(reactionSet),
    _propensities(state.getNumberOfReactions()),
    // Random number generators.
    _discreteUniformGenerator(),
    _exponentialGenerator(&_discreteUniformGenerator),
    _discreteFiniteGenerator(&_discreteUniformGenerator),
    // Direct method.
    _exponentialDeviate(-1),
    _directStepCount(0),
    // CONTINUE: This is reasonable, but I should experiment a bit.
    _volatileLimit(1.0 / epsilon),
    // Tau-leaping.
    _tauLeaping(state, reactionSet, &_discreteUniformGenerator, epsilon),
    _tauLeapingStepCount(0),
    _p() {
  }

  //! Destruct.
  ~HybridDirectTauLeaping() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Simulation.
  //@{
public:

  //! Seed the random number generator.
  void
  seed(typename DiscreteUniformGenerator::result_type s) {
    _discreteUniformGenerator.seed(s);
  }

  //! Initialize the state with the initial populations and time.
  void
  initialize(const typename State::PopulationsContainer& populations,
	     const Number time = 0);

#if 0
  //! Simulate until the specified time is reached.
  void
  simulate(const Number time) {
    EssTerminationConditionEndTime<Number> terminationCondition(time);
    simulate(terminationCondition);
  }

  //! Fire the specified number of reactions.
  void
  simulate(const std::size_t reactionLimit) {
    EssTerminationConditionReactionCount<Number>
      terminationCondition(reactionLimit);
    simulate(terminationCondition);
  }

  //! Simulate until no more reactions can fire.
  void
  simulate() {
    EssTerminationConditionExhaust<Number> terminationCondition;
    simulate(terminationCondition);
  }
#endif

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateForward(_TerminationCondition terminationCondition) {
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(&HybridDirectTauLeaping::stepForward, terminationCondition)) {
    }
  }

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateMidpoint(_TerminationCondition terminationCondition) {
    _p.resize(_state.getNumberOfSpecies());
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(&HybridDirectTauLeaping::stepMidpoint, terminationCondition)) {
    }
  }

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateRungeKutta4(_TerminationCondition terminationCondition) {
    _p.resize(_state.getNumberOfSpecies());
    _k1.resize(_state.getNumberOfReactions());
    _k2.resize(_state.getNumberOfReactions());
    _k3.resize(_state.getNumberOfReactions());
    _k4.resize(_state.getNumberOfReactions());
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(&HybridDirectTauLeaping::stepRungeKutta4,
		terminationCondition)) {
    }
  }

private:

  //! Try to take a step.  Return true if a step is taken.
  /*! Multiply _propensities by the step, tau. */
  template<typename _TerminationCondition>
  bool
  step(MemberFunctionPointer method, 
       _TerminationCondition terminationCondition);

  bool
  stepForward(Number tau);

  bool
  stepMidpoint(Number tau);

  bool
  stepRungeKutta4(Number tau);

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

  //! Return a const reference to the discrete, uniform generator.
  const DiscreteUniformGenerator&
  getDiscreteUniformGenerator() const {
    return _discreteUniformGenerator;
  }

  //! Get the number of direct steps.
  std::size_t
  getDirectCount() const {
    return _directStepCount;
  }

  //! Get the number of tau-leaping steps.
  std::size_t
  getTauLeapingCount() const {
    return _tauLeapingStepCount;
  }

private:
  
  // Return true if the reaction is volatile.
  bool
  isVolatile(int index);

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Return a reference to the discrete, uniform generator.
  DiscreteUniformGenerator&
  getDiscreteUniformGenerator() {
    return _discreteUniformGenerator;
  }
  
private:

  // Move the volatile and slow reactions to the direct group.
  void
  moveVolatileAndSlowToDirect(Number minimumPropensity);

  // Compute the sum of the PMF for the direct group reactions.
  Number
  computeDirectPmfSum() {
    const std::vector<int>& events = _discreteFiniteGenerator.getEvents();
    Number sum = 0;
    for (std::size_t i = 0; i != events.size(); ++i) {
      sum += _propensities[events[i]];
    }
    return sum;
  }

  // CONTINUE
#if 0
  // Compute the time to the next direct reaction.
  Number
  computeDirectStep() {
    Number pmfSum = 0;
    for (std::size_t i = 0; i != _discreteFiniteGenerator.getEvents().size();
	 ++i) {
      pmfSum += _propensities[_discreteFiniteGenerator.getEvents()[i]];
    }
    // If any reactions can fire.
    if (pmfSum != 0) {
      // Return the time to the next reaction.
      return _exponentialGenerator() / pmfSum;
    }
    // Otherwise return infinity.
    return std::numeric_limits<Number>::max();
  }
#endif

  void
  computePropensities() {
    computePropensities(_state.getPopulations());
  }

  void
  computePropensities(const PopulationsContainer& populations) {
    for (std::size_t m = 0; m < _propensities.size(); ++m) {
      _propensities[m] = _reactionSet.computePropensity(m, populations);
    }
  }

  //@}
};

//@}
  
END_NAMESPACE_STOCHASTIC

#define __stochastic_HybridDirectTauLeaping_ipp__
#include "HybridDirectTauLeaping.ipp"
#undef __stochastic_HybridDirectTauLeaping_ipp__

#endif
