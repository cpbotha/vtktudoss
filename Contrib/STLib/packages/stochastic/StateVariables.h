// -*- C++ -*-

/*! 
  \file stochastic/StateVariables.h
  \brief The state of the stochastic simulation.
*/

#if !defined(__stochastic_StateVariables_h__)
#define __stochastic_StateVariables_h__

#include "defs.h"

#include "../ext/vector.h"

BEGIN_NAMESPACE_STOCHASTIC

//----------------------------------------------------------------------------
//! The state variable for the stochastic simulation.
/*!
  Hold the time, the reaction counts, and the populations.

  \param _PopulationType is the number type for storing populations. By 
  default this is double.
  \param _ReactionCountType The number type for storing populations. By 
  default this is double.
*/
template<typename _PopulationType = double,
	 typename _ReactionCountType = double>
class StateVariables {
  //
  // Public types.
  //
public:

  //! The number type for storing populations.
  typedef _PopulationType PopulationType;
  //! The number type for storing reaction counts.
  typedef _ReactionCountType ReactionCountType;
  //! The size type.
  typedef std::size_t SizeType;
  //! The floating-point number type.
  typedef double Number;
  //! The container for the populations.
  typedef std::vector<PopulationType> PopulationsContainer;

  //
  // Member data.
  //
protected:

  //! The time.
  Number _time;
  //! The total number of reaction firings.
  ReactionCountType _totalReactionCount;
  //! The populations of the species.
  PopulationsContainer _populations;
  //! The number of reaction firings.
  std::vector<ReactionCountType> _reactionCounts;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Default constructor.
  StateVariables() :
    _time(0),
    _totalReactionCount(0),
    _populations(),
    _reactionCounts() {
  }

  //! Construct from the number of species and the state change vectors.
  StateVariables(const SizeType numberOfSpecies,
		 const SizeType numberOfReactions) :
    _time(0),
    _totalReactionCount(0),
    // Initialize the populations to zero.
    _populations(numberOfSpecies, 0),
    // Initialize the reaction counts.
    _reactionCounts(numberOfReactions, ReactionCountType(0)) {
  }

  //! Construct from the populations and the state change vectors.
  StateVariables(const PopulationsContainer& populations,
		 const SizeType numberOfReactions) :
    _time(0),
    _totalReactionCount(0),
    // Copy the populations.
    _populations(populations),
    // Initialize the reaction counts.
    _reactionCounts(numberOfReactions, ReactionCountType(0)) {
  }

  //! Copy constructor.
  StateVariables(const StateVariables& other) :
    _time(other._time),
    _totalReactionCount(other._totalReactionCount),
    _populations(other._populations),
    _reactionCounts(other._reactionCounts) {
  }

  //! Assignment operator.
  StateVariables&
  operator=(const StateVariables& other) {
    if (this != &other) {
      _time = other._time;
      _totalReactionCount = other._totalReactionCount;
      _populations = other._populations;
      _reactionCounts = other._reactionCounts;
    }
    return *this;
  }

  //! Destructor.
  ~StateVariables() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Get the time.
  Number
  getTime() const {
    return _time;
  }

  //! Get the total number of reaction firings.
  ReactionCountType
  getReactionCount() const {
    return _totalReactionCount;
  }

  //! Get the number of reaction firings for the specified reaction.
  ReactionCountType
  getReactionCount(const SizeType reactionIndex) const {
    return _reactionCounts[reactionIndex];
  }

  //! Get the vector of reaction counts.
  const std::vector<ReactionCountType>&
  getReactionCounts() const {
    return _reactionCounts;
  }

  //! Get the number of species.
  SizeType
  getNumberOfSpecies() const {
    return _populations.size();
  }

  //! Get the number of reactions.
  SizeType
  getNumberOfReactions() const {
    return _reactionCounts.size();
  }

  //! Get the populations.
  const PopulationsContainer&
  getPopulations() const {
    return _populations;
  }

  //! Get the specified population.
  PopulationType
  getPopulation(const SizeType speciesIndex) const {
    return _populations[speciesIndex];
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Validity.
  //@{
public:

  //! Return true if the state is valid.
  bool
  isValid() const {
    // Check that there are no negative populations.
    for (SizeType i = 0; i != _populations.size(); ++i) {
      if (_populations[i] < 0) {
	return false;
      }
    }
    return true;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Set the time.
  void
  setTime(const Number time) {
    _time = time;
  }

  //! Advance the time and return the new time.
  Number
  advanceTime(const Number increment) {
    return _time += increment;
  }

  //! Reset the reaction counts.
  void
  resetReactionCounts() {
    _totalReactionCount = 0;
    std::fill(_reactionCounts.begin(), _reactionCounts.end(), 0);
  }

  //! Set the populations.
  void
  setPopulations(const PopulationsContainer& populations) {
    assert(_populations.size() == populations.size());
    // Copy the populations.
    std::copy(populations.begin(), populations.end(), _populations.begin());
  }

  //! Set the specified population.
  void
  setPopulation(const SizeType index, const PopulationType population) {
    _populations[index] = population;
  }

  //! Offset the populations.
  void
  offsetPopulations(const PopulationsContainer& change) {
    _populations += change;
  }

#if 0
  //! Increment the specified reaction counts.
  void
  incrementReactionCounts(const SizeType reactionIndex,
			  const ReactionCountType numberOfTimes) {
    _totalReactionCount += numberOfTimes;
    _reactionCounts[reactionIndex] += numberOfTimes;
  }
#endif

  //! Fix any negative populations (by making them zero).
  void
  fixNegativePopulations() {
    // Fix the negative populations.
    const SizeType end = _populations.size();
    for (SizeType n = 0; n != end; ++n) {
      if (_populations[n] < 0) {
	_populations[n] = 0;
      }
    }
  }

  //@}
};

//--------------------------------------------------------------------------
//! \defgroup stochastic_StateVariablesFunctions Free functions for StateVariables.
//@{

//! Return true if the states are equal.
/*! \relates StateVariables */
template<typename _PType, typename _RCType>
bool
operator==(const StateVariables<_PType, _RCType>& x,
	   const StateVariables<_PType, _RCType>& y);

//! Return true if the states are not equal.
/*! \relates StateVariables */
template<typename _PType, typename _RCType>
inline
bool
operator!=(const StateVariables<_PType, _RCType>& x,
	   const StateVariables<_PType, _RCType>& y) {
  return !(x == y);
}

//! Write the populations in ascii format.
/*! \relates StateVariables */
template<typename _PType, typename _RCType>
void
writePopulationsAscii(std::ostream& out, 
		      const StateVariables<_PType, _RCType>& x);

//! Read the populations in ascii format.
/*! \relates StateVariables */
template<typename _PType, typename _RCType>
void
readPopulationsAscii(std::istream& in, StateVariables<_PType, _RCType>* x);

//! Write the state in ascii format.
/*! \relates StateVariables */
template<typename _PType, typename _RCType>
void
writeAscii(std::ostream& out, const StateVariables<_PType, _RCType>& x);

//@}

END_NAMESPACE_STOCHASTIC

#define __stochastic_StateVariables_ipp__
#include "StateVariables.ipp"
#undef __stochastic_StateVariables_ipp__

#endif
