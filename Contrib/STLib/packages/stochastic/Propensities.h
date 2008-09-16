// -*- C++ -*-

/*! 
  \file stochastic/Propensities.h
  \brief Functors for computing propensities.
*/

#if !defined(__stochastic_Propensities_h__)
#define __stochastic_Propensities_h__

#include "ReactionSet.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_Propensities)
#define DEBUG_stochastic_Propensities
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Functor for computing a single propensity.
/*!
  \param _Number The number type.

  \note This functor is expensive to copy.
*/
template<typename _Number = double>
class PropensitiesSingle :
  private ReactionSet<_Number> {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _Number Number;
  //! A set of reactions.
  typedef ReactionSet<Number> ReactionSet;
  //! The result type.
  typedef Number result_type;

  //
  // Private types.
  //
private:

  typedef ReactionSet Base;

  //
  // Not implemented.
  //
private:  

  // Default constructor not implemented.
  PropensitiesSingle();
  //! Assignment operator not implemented.
  PropensitiesSingle&
  operator=(const PropensitiesSingle&);


  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the set of reactions.
  PropensitiesSingle(const ReactionSet& reactionSet) :
    Base(reactionSet) {
  }

  //! Copy constructor.
  PropensitiesSingle(const PropensitiesSingle& other) :
    Base(other) {
  }

  //! Destructor.
  ~PropensitiesSingle() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Functor.
  //@{
public:

  //! Return the specified propensity function.
  template<typename Container>
  result_type
  operator()(const int n, const Container& populations) const {
    return Base::computePropensity(n, populations);
  }

  //@}
};

  
//! Functor for computing all of the propensities.
/*!
  \param _Number The number type.

  \note This functor is expensive to copy.
*/
template<typename _Number = double>
class PropensitiesAll :
  private ReactionSet<_Number> {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _Number Number;
  //! A set of reactions.
  typedef ReactionSet<Number> ReactionSet;
  //! The result type.
  typedef void result_type;

  //
  // Private types.
  //
private:

  typedef ReactionSet Base;

  //
  // Not implemented.
  //
private:  

  // Default constructor not implemented.
  PropensitiesAll();
  //! Assignment operator not implemented.
  PropensitiesAll&
  operator=(const PropensitiesAll&);


  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the set of reactions.
  PropensitiesAll(const ReactionSet& reactionSet) :
    Base(reactionSet) {
  }

  //! Copy constructor.
  PropensitiesAll(const PropensitiesAll& other) :
    Base(other) {
  }

  //! Destructor.
  ~PropensitiesAll() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Functor.
  //@{
public:

  //! Compute the propensity functions.
  template<typename _Container, typename _RandomAccessIterator>
  result_type
  operator()(const _Container& populations,
	     _RandomAccessIterator propensities) const {
    for (int i = 0; i != Base::getSize(); ++i) {
      propensities[i] = Base::computePropensity(i, populations);
    }
  }

  //@}
};

  
END_NAMESPACE_STOCHASTIC

#endif
