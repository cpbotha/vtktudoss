// -*- C++ -*-

/*! 
  \file stochastic/ReactionSet.h
  \brief The state of the stochastic simulation.
*/

#if !defined(__stochastic_ReactionSet_h__)
#define __stochastic_ReactionSet_h__

#include "Reaction.h"

#include "../ads/array/SparseArray.h"
#include "../ext/vector.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_ReactionSet)
#define DEBUG_stochastic_ReactionSet
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! The state of the stochastic simulation.
/*!
  \param T The number type.
*/
template<typename T = double>
class ReactionSet {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef T Number;
  //! The reaction type.
  typedef Reaction<Number> Reaction;

  //
  // Private types.
  //
private:

  typedef std::vector<Reaction> ReactionContainer;

  //
  // More public types.
  //
public:

  //! A const iterator on reactions.
  typedef typename ReactionContainer::const_iterator ReactionConstIterator;

  //
  // Protected types.
  //
protected:

  //! A sparse array of integers.
  typedef ads::SparseArray<1,int> SparseArrayInt;

  //
  // Member data.
  //
private:

  //! The reactions.
  ReactionContainer _reactions;

  //
  // Not implemented.
  //
private:  

  //! Assignment operator not implemented.
  ReactionSet&
  operator=(const ReactionSet&);


  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Default constructor.
  ReactionSet() :
    _reactions()
  {}

  //! Construct from a range of reactions.
  template<typename _InputIterator>
  ReactionSet(_InputIterator reactionsBeginning, _InputIterator reactionsEnd) :
    _reactions(reactionsBeginning, reactionsEnd) {
  }

  //! Copy constructor.
  ReactionSet(const ReactionSet& other) :
    _reactions(other._reactions)
  {}

  //! Destructor.
  ~ReactionSet() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Get the number of reactions.
  std::size_t
  getSize() const {
    return _reactions.size();
  }

  //! Get the specified reaction.
  const Reaction&
  getReaction(const int n) const {
    return _reactions[n];
  }

  //! Get the beginning of the range of reactions.
  ReactionConstIterator
  getBeginning() const {
    return _reactions.begin();
  }

  //! Get the end of the range of reactions.
  ReactionConstIterator
  getEnd() const {
    return _reactions.end();
  }

  //! Return the specified propensity function.
  template<typename Container>
  Number
  computePropensity(const int n, const Container& populations) const {
    return _reactions[n].computePropensityFunction(populations);
  }

  // CONTINUE
#if 0
  //! Compute the propensity functions.
  template<typename _Container, typename _RandomAccessIterator>
  void
  computePropensities(const _Container& populations,
		      _RandomAccessIterator propensities) const {
    for (std::size_t i = 0; i != _reactions.size(); ++i) {
      propensities[i] = _reactions[i].computePropensityFunction(populations);
    }
  }
#endif

  //! Compute the number of species from the reactions.
  int
  computeNumberOfSpecies() const {
    int maxIndex = -1;
    for (std::size_t i = 0; i != _reactions.size(); ++i) {
      // computeMaximum checks for empty arrays.
      maxIndex = std::max(maxIndex, ads::computeMaximum
			  (_reactions[i].getReactants().getIndices()));
      maxIndex = std::max(maxIndex, ads::computeMaximum
			  (_reactions[i].getProducts().getIndices()));
    }
    return maxIndex + 1;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Set the rate constants.
  template<typename _ForwardIterator>
  void
  setRateConstants(_ForwardIterator beginning, _ForwardIterator end) {
#ifdef DEBUG_stochastic
    assert(std::distance(beginning, end) == getSize());
#endif
    for (std::size_t i = 0; i != getSize(); ++i, ++beginning) {
      _reactions[i].setRateConstant(*beginning);
    }
  }

  //! Set all of the rate constants to a value.
  void
  setRateConstants(const Number value) {
    for (std::size_t i = 0; i != getSize(); ++i) {
      _reactions[i].setRateConstant(value);
    }
  }

  //! Set the specified rate constant.
  void
  setRateConstant(const std::size_t n, const Number value) {
    _reactions[n].setRateConstant(value);
  }

  //! Build from the reactions.
  template<typename _InputIterator>
  void
  rebuild(_InputIterator reactionsBeginning, _InputIterator reactionsEnd);

  //@}
};

  
//--------------------------------------------------------------------------
//! \defgroup stochastic_ReactionSetFunctions Free functions for ReactionSet.
//@{

//! Return true if the states are equal.
/*! \relates ReactionSet */
template<typename T>
bool
operator==(const ReactionSet<T>& x, const ReactionSet<T>& y);

//! Return true if the states are not equal.
/*! \relates ReactionSet */
template<typename T>
inline
bool
operator!=(const ReactionSet<T>& x, const ReactionSet<T>& y) {
  return !(x == y);
}

//! Write the reactions in ascii format.
/*! \relates ReactionSet */
template<typename T>
void
writeAscii(std::ostream& out, const ReactionSet<T>& x);

//! Read the reactions in ascii format.
/*! \relates ReactionSet */
template<typename T>
void
readAscii(std::istream& in, ReactionSet<T>* x);

//! Read the reactants and products in ascii format.
/*! \relates ReactionSet */
template<typename T>
void
readReactantsAndProductsAscii(std::istream& in, ReactionSet<T>* x);

//! Read the reactants and products in ascii format.
/*! \relates ReactionSet */
template<typename T>
void
readReactantsAndProductsAscii(std::istream& in, 
			      const std::size_t numberOfReactions, 
			      ReactionSet<T>* x);

//! Read the rate constants in ascii format.
/*! \relates ReactionSet */
template<typename T>
void
readRateConstantsAscii(std::istream& in, ReactionSet<T>* x);

//@}

END_NAMESPACE_STOCHASTIC

#define __stochastic_ReactionSet_ipp__
#include "ReactionSet.ipp"
#undef __stochastic_ReactionSet_ipp__

#endif
