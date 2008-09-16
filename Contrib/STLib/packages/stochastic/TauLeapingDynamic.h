// -*- C++ -*-

/*!
  \file stochastic/TauLeapingDynamic.h
  \brief The tau-leaping method for SSA.
*/

#if !defined(__stochastic_TauLeapingDynamic_h__)
#define __stochastic_TauLeapingDynamic_h__

#include "State.h"
#include "ReactionSet.h"

#include "../numerical/random/poisson/PoissonGeneratorInvAcNormSure.h"
#include "../ads/array/StaticArrayOfArrays.h"

#include <set>

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_TauLeapingDynamic)
#define DEBUG_stochastic_TauLeapingDynamic
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Perform a stochastic simulation using the tau-leaping method.
/*!
  \param _State The state of the simulation: reactions, populations, and time.
*/
template<class _State>
class TauLeapingDynamic {
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
  // Member data.  
  //
private:

  //State _state;
  //ReactionSet _reactionSet;
  //std::vector<Number> _propensities;
  //DiscreteUniformGenerator _discreteUniformGenerator;
  NormalGenerator _normalGenerator;
  PoissonGenerator _poissonGenerator;

  // The active reactions.
  std::vector<int> _activeReactions;
  // The reactants for each reaction.
  ads::StaticArrayOfArrays<int> _reactants;
  // The active species (reactants of the active reactions).
  std::vector<int> _activeSpecies;
  // The mean population change.
  std::vector<Number> _mu;
  // The variance in the population change.
  std::vector<Number> _sigmaSquared;
  std::vector<int> _highestOrder;
  std::vector<int> _highestIndividualOrder;
  Number _epsilon;

  //
  // Not implemented.
  //
private:

  //! Default constructor not implemented.
  TauLeapingDynamic();
  //! Copy constructor not implemented.
  TauLeapingDynamic(const TauLeapingDynamic&);
  //! Assignment operator not implemented.
  TauLeapingDynamic&
  operator=(const TauLeapingDynamic&);
  
  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct.
  TauLeapingDynamic(const State& state, const ReactionSet& reactionSet,
		    DiscreteUniformGenerator* discreteUniformGenerator,
		    const Number epsilon) :
    _normalGenerator(discreteUniformGenerator),
    // CONTINUE: normal threshhold.
    _poissonGenerator(&_normalGenerator, 1000),
    // Initially there are no active reactions.
    _activeReactions(),
    _reactants(),
    _activeSpecies(),
    _mu(state.getNumberOfSpecies()),
    _sigmaSquared(state.getNumberOfSpecies()),
    _highestOrder(state.getNumberOfSpecies()),
    _highestIndividualOrder(state.getNumberOfSpecies()),
    _epsilon() {
    initialize(reactionSet);
    // Set epsilon, the normal threshhold, and the sure threshhold.
    setEpsilon(epsilon);
  }

  //! Destruct.
  ~TauLeapingDynamic()
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Simulation.
  //@{
public:

  //! Compute the maximum allowed time step.
  Number
  computeStep(const std::vector<ads::SparseArray<1, PopulationType> >& 
	      listOfStateChangeVectors,
	      const std::vector<Number>& propensities,
	      const std::vector<PopulationType>& populations);

  //! Generate a Poisson deviate with the specified mean.
  int
  generatePoisson(const Number mean) {
    return _poissonGenerator(mean);
  }

private:

  //! Compute the orders for the species.
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
  initialize(const ReactionSet& reactionSet);

  //! Compute mu and sigma squared.
  void
  computeMuAndSigmaSquared
  (const std::vector<ads::SparseArray<1, PopulationType> >& 
   listOfStateChangeVectors,
   const std::vector<Number>& propensities);

  //! Compute the g described in "Efficient step size selection for the tau-leaping simulation method".
  Number
  computeG(int speciesIndex, PopulationType population) const;

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Get the value of epsilon.
  Number
  getEpsilon() const {
    return _epsilon;
  }

  //! Get the active reactions.
  const std::vector<int>&
  getActiveReactions() const {
    return _activeReactions;
  }

  //! Return true if there are no active reactions.
  bool
  empty() const {
    return _activeReactions.empty();
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Set the value of epsilon.
  void
  setEpsilon(const Number epsilon) {
    _epsilon = epsilon;
    // The relative error in the mean is less than 0.1 * error.
    // continuityError / mean < 0.1 * error
    // 1 / mean < 0.1 * error
    // mean > 10 / error
    const Number t = 10. / epsilon;
    _poissonGenerator.setNormalThreshhold(t);
    // The relative error in neglecting the standard deviation is less 
    // than 0.1 * error.
    // sqrt(mean) / mean < 0.1 * error
    // mean > 100 / error^2
    _poissonGenerator.setSureThreshhold(t * t);
  }

  //! Clear the active set.
  void
  clear() {
    _activeReactions.clear();
    _activeSpecies.clear();
  }

  //! Insert a reaction to the active set.
  void
  insert(const int index) {
    _activeReactions.push_back(index);
    computeActiveSpecies();
  }

  //! Remove a reaction from the active set.
  void
  erase(const int index) {
    const std::ptrdiff_t i = std::find(_activeReactions.begin(), 
				       _activeReactions.end(), index)
      - _activeReactions.begin();
#ifdef DEBUG_stlib
    assert(std::size_t(i) != _activeReactions.size());
#endif
    _activeReactions[i] = _activeReactions.back();
    _activeReactions.pop_back();
    computeActiveSpecies();
  }

private:

  //! Compute the active species from the active reactions.
  void
  computeActiveSpecies() {
    typedef ads::StaticArrayOfArrays<int>::const_iterator const_iterator;
    _activeSpecies.clear();
    std::set<int> active;
    for (std::vector<int>::const_iterator i = _activeReactions.begin(); 
	 i != _activeReactions.end(); ++i) {
      for (const_iterator j = _reactants.begin(*i); j != _reactants.end(*i); 
	   ++j) {
	active.insert(*j);
      }
    }
    std::copy(active.begin(), active.end(), 
	      std::back_inserter(_activeSpecies));
  }

  //@}
};

END_NAMESPACE_STOCHASTIC

#define __stochastic_TauLeapingDynamic_ipp__
#include "TauLeapingDynamic.ipp"
#undef __stochastic_TauLeapingDynamic_ipp__

#endif
