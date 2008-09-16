// -*- C++ -*-

/*! 
  \file stochastic/Direct.h
  \brief The direct method for SSA.
*/

#if !defined(__stochastic_Direct_h__)
#define __stochastic_Direct_h__

#include "EssTerminationCondition.h"
#include "State.h"
#include "Propensities.h"

#include "../ads/array/StaticArrayOfArrays.h"
#include "../numerical/random/exponential/ExponentialGeneratorZiggurat.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_Direct)
#define DEBUG_stochastic_Direct
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Perform a stochastic simulation using Gillespie's direct method.
/*!
  \param _DiscreteFiniteGenerator Random deviate generator for the discrete,
  finite distribution with reaction propensities as scaled probabilities.
  \param _ExponentialGenerator Random deviate generator for the exponential
  distribution. By default the ziggurat algorithm is used.
  \param _PropensitiesFunctor Can calculate propensities as a function of the 
  reaction index and the populations.
  \param _State The state of the simulation: reactions, populations, and time.
*/
template<class _DiscreteFiniteGenerator,
	 class _ExponentialGenerator = 
	 numerical::ExponentialGeneratorZiggurat<>,
	 class _PropensitiesFunctor = PropensitiesSingle<>,
	 class _State = State<> >
class Direct {
  //
  // Public types.
  //
public:

  //! The state.
  typedef _State State;
  //! The propensities functor.
  typedef _PropensitiesFunctor PropensitiesFunctor;
  //! The exponential generator.
  typedef _ExponentialGenerator ExponentialGenerator;
  //! The discrete, finite generator.
  typedef _DiscreteFiniteGenerator DiscreteFiniteGenerator;
  //! The discrete, uniform generator.
  typedef typename ExponentialGenerator::DiscreteUniformGenerator
  DiscreteUniformGenerator;

  //! The number type.
  typedef typename State::Number Number;

  //
  // Enumerations.
  //
private:

  enum {ComputeIndividualPropensities = 
	Loki::IsSameType<typename PropensitiesFunctor::result_type, 
	Number>::value};

  //
  // Member data.
  //
private:

  State _state;
  PropensitiesFunctor _propensitiesFunctor;
  ads::StaticArrayOfArrays<int> _reactionInfluence;
  DiscreteUniformGenerator _discreteUniformGenerator;
  ExponentialGenerator _exponentialGenerator;
  DiscreteFiniteGenerator _discreteFiniteGenerator;
  std::vector<Number> _propensities;
  Number _tau;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  Direct();
  //! Copy constructor not implemented.
  Direct(const Direct&);
  //! Assignment operator not implemented.
  Direct&
  operator=(const Direct&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct.
  Direct(const State& state,
	 PropensitiesFunctor propensitiesFunctor,
	 const ads::StaticArrayOfArrays<int>& reactionInfluence) :
    // Copy.
    _state(state),
    _propensitiesFunctor(propensitiesFunctor),
    _reactionInfluence(reactionInfluence),
    // Construct.
    _discreteUniformGenerator(),
    _exponentialGenerator(&_discreteUniformGenerator),
    _discreteFiniteGenerator(&_discreteUniformGenerator),
    // Invalid value.
    _tau(-1) {
    // Allocate memory if needed.
    if (! ComputeIndividualPropensities) {
      _propensities.resize(state.getNumberOfReactions());
    }
  }

  //! Destruct.
  ~Direct() {
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
	     const Number time = 0) {
    // Initialize the state.
    _state.setPopulations(populations);
    _state.setTime(time);
    _state.resetReactionCounts();

    // Compute the propensities and initialize the discrete, finite generator.
    computePropensities();
    // Compute the initial time step.
    _tau = computeTau();
  }

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

  //! Simulate until the termination condition is reached.
  template<typename _TerminationCondition>
  void
  simulate(_TerminationCondition terminationCondition) {
    // Step until no more reactions can fire or we reach the termination
    // condition.
    while (step(terminationCondition)) {
    }
  }

  //! Try to take a step.  Return true if a step is taken.
  template<typename _TerminationCondition>
  bool
  step(_TerminationCondition terminationCondition);

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
  
  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Return a reference to the discrete, finite generator.
  DiscreteFiniteGenerator&
  getDiscreteFiniteGenerator() {
    return _discreteFiniteGenerator;
  }
  
  //! Return a reference to the discrete, uniform generator.
  DiscreteUniformGenerator&
  getDiscreteUniformGenerator() {
    return _discreteUniformGenerator;
  }
  
private:

  // Compute the time to the next reaction.
  Number
  computeTau() {
    const Number pmfSum = _discreteFiniteGenerator.getPmfSum();
    // If any reactions can fire.
    if (pmfSum != 0) {
      // Return the time to the next reaction.
      return _exponentialGenerator() / pmfSum;
    }
    // Otherwise return infinity.
    return std::numeric_limits<Number>::max();
  }

  void
  computePropensities() {
    _computePropensities(Loki::Int2Type<ComputeIndividualPropensities>());
  }

  void
  _computePropensities(Loki::Int2Type<false> /*dummy*/) {
    // Use the member function propensities array.
    _propensitiesFunctor(_state.getPopulations(), _propensities.begin());
    _discreteFiniteGenerator.initialize(_propensities.begin(), 
					_propensities.end());
  }

  void
  _computePropensities(Loki::Int2Type<true> /*dummy*/) {
    // Allocate a propensity array.
    std::vector<Number> propensities(_state.getNumberOfReactions());
    // Compute each propensity.
    for (std::size_t i = 0; i != propensities.size(); ++i) {
      propensities[i] = _propensitiesFunctor(i, _state.getPopulations());
    }
    _discreteFiniteGenerator.initialize(propensities.begin(), 
					propensities.end());
  }

  void
  updatePropensities(const int reactionIndex) {
    _updatePropensities
      (Loki::Int2Type<ComputeIndividualPropensities>(), reactionIndex);
  }

  void
  _updatePropensities(Loki::Int2Type<false> /*dummy*/,
		      const int reactionIndex) {
    // Recompute all of the propensity functions.
    _propensitiesFunctor(_state.getPopulations(), _propensities.begin());
    _discreteFiniteGenerator.setPmf(_propensities.begin());
    _discreteFiniteGenerator.updatePmf();
  }

  void
  _updatePropensities(Loki::Int2Type<true> /*dummy*/,
		      const int reactionIndex) {
    for (typename ads::StaticArrayOfArrays<int>::const_iterator 
	   i = _reactionInfluence.begin(reactionIndex); 
	 i != _reactionInfluence.end(reactionIndex); ++i) {
      _discreteFiniteGenerator.setPmf
	(*i, _propensitiesFunctor(*i, _state.getPopulations()));
    }
    _discreteFiniteGenerator.updatePmf();
  }

  //@}
};

//@}
  
END_NAMESPACE_STOCHASTIC

#define __stochastic_Direct_ipp__
#include "Direct.ipp"
#undef __stochastic_Direct_ipp__

#endif
