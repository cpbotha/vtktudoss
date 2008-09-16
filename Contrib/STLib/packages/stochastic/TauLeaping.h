// -*- C++ -*-

/*! 
  \file stochastic/TauLeaping.h
  \brief The tau-leaping method for SSA.
*/

#if !defined(__stochastic_TauLeaping_h__)
#define __stochastic_TauLeaping_h__

#include "State.h"
#include "ReactionSet.h"
#include "EssTerminationCondition.h"

#include "../ads/array/StaticArrayOfArrays.h"
#include "../numerical/random/poisson/PoissonGeneratorInvAcNormSure.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_TauLeaping)
#define DEBUG_stochastic_TauLeaping
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Perform a stochastic simulation using the tau-leaping method.
/*!
  \param _State The state of the simulation: reactions, populations, and time.
  \param _MinimumAllowedStep If true, the tau will be chosen so that on 
  average at least one reaction fires.
*/
template<class _State, bool _MinimumAllowedStep = false>
class TauLeaping {
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
  //! The Poisson generator.
  typedef numerical::PoissonGeneratorInvAcNormSure<> PoissonGenerator;
  //! The discrete uniform generator.
  typedef typename PoissonGenerator::DiscreteUniformGenerator
  DiscreteUniformGenerator;
  //! The normal generator.
  typedef typename PoissonGenerator::NormalGenerator NormalGenerator;

  //
  // Private types.
  //
private:

  typedef bool (TauLeaping::*MemberFunctionPointer)(Number);

  //
  // Member data.  
  //
private:

  State _state;
  ReactionSet _reactionSet;
  std::vector<Number> _propensities;
  DiscreteUniformGenerator _discreteUniformGenerator;
  NormalGenerator _normalGenerator;
  PoissonGenerator _poissonGenerator;

  // The mean population change.
  std::vector<Number> _mu;
  // The variance in the population change.
  std::vector<Number> _sigmaSquared;
  std::vector<int> _highestOrder;
  std::vector<int> _highestIndividualOrder;
  //! The number of steps.
  std::size_t _stepCount;
  //! Used for Runge-Kutta 4 and midpoint.
  PopulationsContainer _p;
  //! Used for Runge-Kutta 4.
  std::vector<Number> _k1, _k2, _k3, _k4;


  //
  // Not implemented.
  //
private:

  //! Default constructor not implemented.
  TauLeaping();
  //! Copy constructor not implemented.
  TauLeaping(const TauLeaping&);
  //! Assignment operator not implemented.
  TauLeaping&
  operator=(const TauLeaping&);
  
  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct.
  TauLeaping(const State& state, const ReactionSet& reactionSet);

  //! Destruct.
  ~TauLeaping()
  {}

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

  //! Initialize the state with the initial populations and time. Compute the orders for the species.
  /*!
    Let \c highestOrder be the highest order reaction in which the species 
    appears.  Let \c highestIndividualOrder be the highest order of the
    species in a reaction.
    Suppose that the reactions are the following.  (We only use the reactants
    in computing the orders.)
    \f[
    x0 \rightarrow \cdots, \quad
    x1 + x2 \rightarrow \cdots, \quad
    x2 + 2 x3 \rightarrow \cdots, \quad
    3 x4 \rightarrow \cdots 
    \f]
    Then the orders are the following.
    \verbatim
    highestOrder == {1, 2, 3, 3, 3}
    highestIndividualOrder == {1, 1, 1, 2, 3}
    \endverbatim
  */
  void
  initialize(const typename State::PopulationsContainer& populations,
	     const Number time = 0);

#if 0
  //! Simulate until the specified time is reached.
  void
  simulate(const Number epsilon, const Number time) {
    EssTerminationConditionEndTime<Number> terminationCondition(time);
    simulate(epsilon, terminationCondition);
  }

  //! Fire the specified number of reactions.
  void
  simulate(const Number epsilon, const std::size_t reactionLimit) {
    EssTerminationConditionReactionCount<Number>
      terminationCondition(reactionLimit);
    simulate(epsilon, terminationCondition);
  }

  //! Simulate until no more reactions can fire.
  void
  simulate(const Number epsilon) {
    EssTerminationConditionExhaust<Number> terminationCondition;
    simulate(epsilon, terminationCondition);
  }
#endif

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateForward(const Number epsilon, 
		  _TerminationCondition terminationCondition) {
    setError(epsilon);
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(&TauLeaping::stepForward, epsilon, terminationCondition)) {
    }
  }

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateMidpoint(const Number epsilon, 
		   _TerminationCondition terminationCondition) {
    setError(epsilon);
    _p.resize(_state.getNumberOfSpecies());
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(&TauLeaping::stepMidpoint, epsilon, terminationCondition)) {
    }
  }

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateRungeKutta4(const Number epsilon, 
		      _TerminationCondition terminationCondition) {
    setError(epsilon);
    _p.resize(_state.getNumberOfSpecies());
    _k1.resize(_state.getNumberOfReactions());
    _k2.resize(_state.getNumberOfReactions());
    _k3.resize(_state.getNumberOfReactions());
    _k4.resize(_state.getNumberOfReactions());
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(&TauLeaping::stepRungeKutta4, epsilon, terminationCondition)) {
    }
  }

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateFixedForward(const Number tau, 
		       _TerminationCondition terminationCondition) {
    setError(0.01);
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (stepFixed(&TauLeaping::stepForward, tau, terminationCondition)) {
    }
  }

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateFixedMidpoint(const Number tau, 
			_TerminationCondition terminationCondition) {
    setError(0.01);
    _p.resize(_state.getNumberOfSpecies());
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (stepFixed(&TauLeaping::stepMidpoint, tau, terminationCondition)) {
    }
  }

  //! Simulate with fixed size steps until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulateFixedRungeKutta4(const Number tau, 
			   _TerminationCondition terminationCondition) {
    setError(0.01);
    _p.resize(_state.getNumberOfSpecies());
    _k1.resize(_state.getNumberOfReactions());
    _k2.resize(_state.getNumberOfReactions());
    _k3.resize(_state.getNumberOfReactions());
    _k4.resize(_state.getNumberOfReactions());
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (stepFixed(&TauLeaping::stepRungeKutta4, tau, terminationCondition)) {
    }
  }

  //! Try to take a step.  Return true if a step is taken.
  template<typename _TerminationCondition>
  bool
  step(MemberFunctionPointer method, const Number epsilon,
       _TerminationCondition terminationCondition);

  //! Try to take a step.  Return true if a step is taken.
  template<typename _TerminationCondition>
  bool
  stepFixed(MemberFunctionPointer method, Number tau,
	    _TerminationCondition terminationCondition);

private:

  void
  setError(const Number error) {
    // The relative error in the mean is less than 0.1 * error.
    // continuityError / mean < 0.1 * error
    // 1 / mean < 0.1 * error
    // mean > 10 / error
    const Number t = 10. / error;
    _poissonGenerator.setNormalThreshhold(t);
    // The relative error in neglecting the standard deviation is less 
    // than 0.1 * error.
    // sqrt(mean) / mean < 0.1 * error
    // mean > 100 / error^2
    _poissonGenerator.setSureThreshhold(t * t);
  }
    

  bool
  stepForward(Number tau);

  bool
  stepMidpoint(Number tau);

  bool
  stepRungeKutta4(Number tau);

  void
  computePropensities() {
    computePropensities(_state.getPopulations());
  }

  void
  computePropensities(const PopulationsContainer& populations);

  Number
  computeTau(Number epsilon);

  //! Compute mu and sigma squared.
  void
  computeMuAndSigmaSquared();

  //! Compute the g described in "Efficient step size selection for the tau-leaping simulation method".
  Number
  computeG(int speciesIndex) const;

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
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Return a reference to the discrete, uniform generator.
  DiscreteUniformGenerator&
  getDiscreteUniformGenerator() {
    return _discreteUniformGenerator;
  }

  //@}
};

END_NAMESPACE_STOCHASTIC

#define __stochastic_TauLeaping_ipp__
#include "TauLeaping.ipp"
#undef __stochastic_TauLeaping_ipp__

#endif
