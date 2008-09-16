// -*- C++ -*-

/*! 
  \file stochastic/NextReaction.h
  \brief The next reaction method for SSA.
*/

#if !defined(__stochastic_NextReaction_h__)
#define __stochastic_NextReaction_h__

#include "EssTerminationCondition.h"
#include "State.h"
#include "Propensities.h"

#include "../ads/array/StaticArrayOfArrays.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_NextReaction)
#define DEBUG_stochastic_NextReaction
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Perform a stochastic simulation using Gibson and Bruck's next reaction method.
/*!
  \param _ReactionPriorityQueue The priority queue for the reactions.
  \param _PropensitiesFunctor Can calculate propensities as a function of the 
  reaction index and the populations.
  \param _State The state of the simulation: reactions, populations, and time.
*/
template<class _ReactionPriorityQueue,
	 class _PropensitiesFunctor = PropensitiesSingle<>,
	 class _State = State<> >
class NextReaction {
  //
  // Public types.
  //
public:

  //! The reaction priority queue.
  typedef _ReactionPriorityQueue ReactionPriorityQueue;
  //! The propensities functor.
  typedef _PropensitiesFunctor PropensitiesFunctor;
  //! The state.
  typedef _State State;

  //! The number type.
  typedef typename State::Number Number;
  //! The discrete uniform generator.
  typedef typename ReactionPriorityQueue::DiscreteUniformGenerator
  DiscreteUniformGenerator;

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
  // Store a pointer here because there are many options for this class.  
  // They require different constructor calls.
  ReactionPriorityQueue* _reactionPriorityQueue;
  ads::StaticArrayOfArrays<int> _reactionInfluence;
  std::vector<Number> _propensities;
  std::vector<Number> _oldPropensities;
  int _reactionIndex;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  NextReaction();
  //! Copy constructor not implemented.
  NextReaction(const NextReaction&);
  //! Assignment operator not implemented.
  NextReaction&
  operator=(const NextReaction&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct.
  NextReaction(const State& state,
	       PropensitiesFunctor propensitiesFunctor,
	       ReactionPriorityQueue* reactionPriorityQueue,
	       const ads::StaticArrayOfArrays<int>& reactionInfluence) :
    // Copy.
    _state(state),
    _propensitiesFunctor(propensitiesFunctor),
    _reactionPriorityQueue(reactionPriorityQueue),
    _reactionInfluence(reactionInfluence),
    // Allocate.
    _propensities(state.getNumberOfReactions()),
    _oldPropensities(ComputeIndividualPropensities ? 0 : _propensities.size()),
    // Invalid index.
    _reactionIndex(-1) {
    // Set the pointer to the propensities array if necessary.
    setPropensities(Loki::Int2Type<ReactionPriorityQueue::UsesPropensities>());
  }

  //! Destruct.
  ~NextReaction() {
  }

private:

  // The indexed priority queue does not use the propensities.
  void
  setPropensities(Loki::Int2Type<false> /*UsesPropensities*/) {
  }

  // The indexed priority queue uses the propensities.
  void
  setPropensities(Loki::Int2Type<true> /*UsesPropensities*/) {
    _reactionPriorityQueue->setPropensities(&_propensities);
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Simulation.
  //@{
public:

  //! Initialize the state with the initial populations and time.
  void
  initialize(const typename State::PopulationsContainer& populations,
	     const Number time = 0) {
    // Initialize the state.
    _state.setPopulations(populations);
    _state.setTime(time);
    _state.resetReactionCounts();

    _reactionIndex = -1;
    // Compute the initial propensities.
    computePropensities();

    // Initialize the reaction priority queue.
    _reactionPriorityQueue->clear();
    for (std::size_t i = 0; i != _propensities.size(); ++i) {
      if (_propensities[i] != 0) {
	_reactionPriorityQueue->push(i, _state.getTime(), _propensities[i]);
      }
    }
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

  //! Return the array of propensities.
  const std::vector<Number>&
  getPropensities() const {
    return _propensities;
  }

  //! Return a const reference to the discrete, uniform generator.
  const DiscreteUniformGenerator&
  getDiscreteUniformGenerator() const {
    return _reactionPriorityQueue->getDiscreteUniformGenerator();
  }
  
  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
private:

  void
  computePropensities() {
    _computePropensities(Loki::Int2Type<ComputeIndividualPropensities>());
  }

  void
  _computePropensities(Loki::Int2Type<false> /*dummy*/) {
    _propensitiesFunctor(_state.getPopulations(), _propensities.begin());
  }

  void
  _computePropensities(Loki::Int2Type<true> /*dummy*/) {
    // Compute each propensity.
    for (std::size_t i = 0; i != _propensities.size(); ++i) {
      _propensities[i] = _propensitiesFunctor(i, _state.getPopulations());
    }
  }

  void
  updatePropensitiesAndReactionTimes(const int reactionIndex) {
    _updatePropensitiesAndReactionTimes
      (Loki::Int2Type<ComputeIndividualPropensities>(), reactionIndex);
  }

  void
  _updatePropensitiesAndReactionTimes(Loki::Int2Type<false> /*dummy*/,
				      const int reactionIndex);

  void
  _updatePropensitiesAndReactionTimes(Loki::Int2Type<true> /*dummy*/,
				      const int reactionIndex);

  //@}
};

//@}
  
END_NAMESPACE_STOCHASTIC

#define __stochastic_NextReaction_ipp__
#include "NextReaction.ipp"
#undef __stochastic_NextReaction_ipp__

#endif
