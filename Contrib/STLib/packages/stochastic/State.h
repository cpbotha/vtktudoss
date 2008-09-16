// -*- C++ -*-

/*! 
  \file stochastic/State.h
  \brief The state of the stochastic simulation.
*/

#if !defined(__stochastic_State_h__)
#define __stochastic_State_h__

#include "StateVariables.h"

#include "../ads/array/SparseArray.h"

BEGIN_NAMESPACE_STOCHASTIC

// CONTINUE: Move scaleAdd here.

//! Build the state change vectors.
template<typename _ForwardIterator, class _Container>
inline
void
buildStateChangeVectors(const std::size_t numberOfSpecies,
			_ForwardIterator reactionsBeginning,
			_ForwardIterator reactionsEnd,
			_Container* stateChangeVectors) {
  typedef typename _Container::value_type::value_type PopulationType;
  // Resize the container using the number of reactions.
  stateChangeVectors->resize(std::distance(reactionsBeginning, reactionsEnd));
  // Construct the state change vectors from the reactions.
  std::vector<PopulationType> reactants(numberOfSpecies),
    products(numberOfSpecies);
  typename _Container::iterator scv = stateChangeVectors->begin();
  for ( ; reactionsBeginning != reactionsEnd; ++reactionsBeginning) {
    reactionsBeginning->getReactants(&reactants);
    reactionsBeginning->getProducts(&products);
    products -= reactants;
    *scv = products;
    ++scv;
  }
}


// CONTINUE: The data structure for the state change vectors could store
// pointers into the populations instead of indices.  This would be a 
// little faster as it would replace array indexing with dereferencing.
//----------------------------------------------------------------------------
//! The state of the stochastic simulation.
/*!
  Hold the time, the reaction count, the populations, and the state change
  vectors.

  \param _PopulationType is the number type for storing populations. By 
  default this is double.
  \param _ReactionCountType The number type for storing populations. By 
  default this is double.
  \param _StateChangeVector is the state change vector.
  By default this is a sparse array of the population type.
*/
template<typename _PopulationType = double,
	 typename _ReactionCountType = double,
	 class _StateChangeVector = ads::SparseArray<1, _PopulationType> >
class State : public StateVariables<_PopulationType, _ReactionCountType> {
  //
  // Private types.
  //
private:


  //
  // Public types.
  //
public:

  //! The base class.
  typedef StateVariables<_PopulationType, _ReactionCountType> Base;
  //! The number type for storing populations.
  typedef typename Base::PopulationType PopulationType;
  //! The number type for storing reaction counts.
  typedef typename Base::ReactionCountType ReactionCountType;
  //! The size type.
  typedef typename Base::SizeType SizeType;
  //! The floating-point number type.
  typedef typename Base::Number Number;
  //! The container for the populations.
  typedef typename Base::PopulationsContainer PopulationsContainer;
  //! The container for the state change vectors.
  typedef std::vector<_StateChangeVector> ScvContainer;

  //
  // Member data.
  //
private:

  using Base::_time;
  using Base::_totalReactionCount;
  using Base::_populations;
  using Base::_reactionCounts;
  //! The state change vectors.
  /*!
    I used to store a const reference to the state change vectors.  That
    enabled multiple states to use the same data structure.  I think that is 
    better to just store a copy.  OpenMP doesn't currently share very well.
  */
  ScvContainer _stateChangeVectors;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  State();
  //! Assignment operator not implemented.
  State&
  operator=(const State&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the number of species and the state change vectors.
  State(const SizeType numberOfSpecies,
	const ScvContainer& stateChangeVectors) :
    Base(numberOfSpecies, stateChangeVectors.size()),
    // Copy the state change vectors.
    _stateChangeVectors(stateChangeVectors) {
  }

  //! Construct from the populations and the state change vectors.
  State(const PopulationsContainer& populations,
	const ScvContainer& stateChangeVectors) :
    Base(populations, stateChangeVectors.size()),
    // Copy the state change vectors.
    _stateChangeVectors(stateChangeVectors) {
  }

  //! Copy constructor.
  State(const State& other) :
    Base(other),
    _stateChangeVectors(other._stateChangeVectors) {
  }

  //! Assignment operator for the base class.
  State&
  operator=(const Base& x) {
    Base::operator=(x);
    return *this;
  }

  //! Destructor.
  ~State() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  using Base::getTime;
  using Base::getReactionCount;
  using Base::getReactionCounts;
  using Base::getNumberOfSpecies;
  using Base::getNumberOfReactions;
  using Base::getPopulations;
  using Base::getPopulation;

  //! Get the state change vectors.
  const ScvContainer&
  getStateChangeVectors() const {
    return _stateChangeVectors;
  }

  //! Get the specified state change vector.
  const typename ScvContainer::value_type&
  getStateChangeVector(const SizeType reactionIndex) const {
    return _stateChangeVectors[reactionIndex];
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Validity.
  //@{
public:

  using Base::isValid;

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  using Base::setTime;
  using Base::advanceTime;
  using Base::resetReactionCounts;
  using Base::setPopulations;
  using Base::setPopulation;
  using Base::offsetPopulations;
  using Base::fixNegativePopulations;

  //! Fire the specified reaction.
  void
  fireReaction(const SizeType n) {
    ++_totalReactionCount;
    ++_reactionCounts[n];
    _populations += _stateChangeVectors[n];
  }

  //! Fire the specified reaction the specified number of times.
  void
  fireReaction(const SizeType reactionIndex, 
	       const ReactionCountType numberOfTimes) {
    _totalReactionCount += numberOfTimes;
    _reactionCounts[reactionIndex] += numberOfTimes;
    scaleAdd(&_populations, numberOfTimes, _stateChangeVectors[reactionIndex]);
  }

  //! Fire the specified reaction the specified number of times to update the specified populations.
  /*! Do not increment the reaction counts. */
  void
  fireReaction(PopulationsContainer* populations, const SizeType reactionIndex,
	       const ReactionCountType numberOfTimes) {
    scaleAdd(populations, numberOfTimes, _stateChangeVectors[reactionIndex]);
  }

  //@}
};

END_NAMESPACE_STOCHASTIC

#define __stochastic_State_ipp__
#include "State.ipp"
#undef __stochastic_State_ipp__

#endif
