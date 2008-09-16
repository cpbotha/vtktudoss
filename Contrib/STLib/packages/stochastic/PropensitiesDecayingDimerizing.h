// -*- C++ -*-

/*! 
  \file stochastic/PropensitiesDecayingDimerizing.h
  \brief The state of the stochastic simulation.
*/

#if !defined(__stochastic_PropensitiesDecayingDimerizing_h__)
#define __stochastic_PropensitiesDecayingDimerizing_h__

#include "defs.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_PropensitiesDecayingDimerizing)
#define DEBUG_stochastic_PropensitiesDecayingDimerizing
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Propensities for the decaying-dimerizing problem.
template<typename _Number = double>
class PropensitiesSingleDecayingDimerizing {
  enum {NumberOfReactions = 4};

  //
  // Public types.
  //
public:

  //! The number type.
  typedef _Number Number;
  //! The result type.
  typedef Number result_type;

  //
  // Private types.
  //
private:

  //! A pointer to a member function that computes a single propensity.
  typedef Number (PropensitiesSingleDecayingDimerizing::* PropensityMember)
  (const int*) const;

  //
  // Member data.
  //
private:

  PropensityMember _propensityFunctions[NumberOfReactions];
  
  //
  // Not implemented.
  //
private:  

  //! Assignment operator not implemented.
  PropensitiesSingleDecayingDimerizing&
  operator=(const PropensitiesSingleDecayingDimerizing&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Default constructor.
  PropensitiesSingleDecayingDimerizing() {
    _propensityFunctions[0] = &PropensitiesSingleDecayingDimerizing::f0;
    _propensityFunctions[1] = &PropensitiesSingleDecayingDimerizing::f1;
    _propensityFunctions[2] = &PropensitiesSingleDecayingDimerizing::f2;
    _propensityFunctions[3] = &PropensitiesSingleDecayingDimerizing::f3;
  }

  //! Copy constructor.
  PropensitiesSingleDecayingDimerizing
  (const PropensitiesSingleDecayingDimerizing& other) {
    for (int i = 0; i != NumberOfReactions; ++i) {
      _propensityFunctions[i] = other._propensityFunctions[i];
    }
  }

  //! Destructor.
  ~PropensitiesSingleDecayingDimerizing() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Functor.
  //@{
public:

#if 0
  //! Get the number of reactions.
  int
  getSize() const {
    return NumberOfReactions;
  }
#endif

  //! Return the specified propensity function.
  template<typename _Container>
  result_type
  operator()(const int n, const _Container& populations) const {
    return (this->*_propensityFunctions[n])(populations.begin());
  }

  //! Compute the propensity functions.
  template<typename _Container, typename _RandomAccessIterator>
  result_type
  operator()(const _Container& populations,
	     _RandomAccessIterator propensities) const {
    propensities[0] = populations[0];
    propensities[1] = (0.002 / 2.0) * populations[0] * (populations[0] - 1);
    propensities[2] = 0.5 * populations[1];
    propensities[3] = 0.04 * populations[1];
  }

  //@}
  //--------------------------------------------------------------------------
  // Compute propensities.
private:

  Number
  f0(const int* x) const {
    return x[0];
  }
  
  Number
  f1(const int* x) const {
    return (0.002 / 2.0) * x[0] * (x[0] - 1);
  }
  
  Number
  f2(const int* x) const {
    return 0.5 * x[1];
  }
  
  Number
  f3(const int* x) const {
    return 0.04 * x[1];
  }
  
};


//! Propensities for the decaying-dimerizing problem.
template<typename _Number = double>
class PropensitiesAllDecayingDimerizing {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _Number Number;
  //! The result type.
  typedef void result_type;

  //--------------------------------------------------------------------------
  //! \name Functor.
  //@{
public:

  //! Compute the propensity functions.
  template<typename _Container, typename _RandomAccessIterator>
  result_type
  operator()(const _Container& populations,
	     _RandomAccessIterator propensities) const {
    propensities[0] = populations[0];
    propensities[1] = (0.002 / 2.0) * populations[0] * (populations[0] - 1);
    propensities[2] = 0.5 * populations[1];
    propensities[3] = 0.04 * populations[1];
  }

  //@}
};

END_NAMESPACE_STOCHASTIC

#endif
