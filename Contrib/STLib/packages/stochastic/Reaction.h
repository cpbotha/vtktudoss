// -*- C++ -*-

/*! 
  \file stochastic/Reaction.h
  \brief A reaction.
*/

#if !defined(__stochastic_Reaction_h__)
#define __stochastic_Reaction_h__

#include "defs.h"

#include "../ads/array/SparseArray.h"

#include "../numerical/specialFunctions/BinomialCoefficient.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_Reaction)
#define DEBUG_stochastic_Reaction
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! A reaction.
/*!
  \param T The number type.
*/
template<typename T = double>
class Reaction {
public:

  //
  // Public types.
  //

  //! The number type.
  typedef T Number;
  //! An array of integers.
  typedef ads::SparseArray<1, int> ArrayInt;

private:

  //
  // Member data.
  //

  //! The reactants.
  ArrayInt _reactants;
  //! The products.
  ArrayInt _products;
  //! Specific reaction probability rate constant.
  Number _rateConstant;

public:

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{

  //! Default constructor.
  Reaction() :
    _reactants(),
    _products(),
    _rateConstant(0)
  {}

  //! Construct from the reactants, the products, and the rate constant.
  template<typename IntegerForwardIterator>
  Reaction(IntegerForwardIterator reactantIndicesBeginning,
	   IntegerForwardIterator reactantIndicesEnd,
	   IntegerForwardIterator reactantCoefficientsBeginning,
	   IntegerForwardIterator reactantCoefficientsEnd,
	   IntegerForwardIterator productIndicesBeginning,
	   IntegerForwardIterator productIndicesEnd,
	   IntegerForwardIterator productCoefficientsBeginning,
	   IntegerForwardIterator productCoefficientsEnd,
	   const Number rateConstant) :
    _reactants(reactantIndicesBeginning, reactantIndicesEnd,
	       reactantCoefficientsBeginning, reactantCoefficientsEnd, 0),
    _products(productIndicesBeginning, productIndicesEnd,
	      productCoefficientsBeginning, productCoefficientsEnd, 0),
    _rateConstant(rateConstant)
  {}

  //! Rebuild from the reactants, the products, and the rate constant.
  template<typename IntegerForwardIterator>
  void
  rebuild(IntegerForwardIterator reactantIndicesBeginning,
	  IntegerForwardIterator reactantIndicesEnd,
	  IntegerForwardIterator reactantCoefficientsBeginning,
	  IntegerForwardIterator reactantCoefficientsEnd,
	  IntegerForwardIterator productIndicesBeginning,
	  IntegerForwardIterator productIndicesEnd,
	  IntegerForwardIterator productCoefficientsBeginning,
	  IntegerForwardIterator productCoefficientsEnd,
	  const Number rateConstant);

  //! Construct from the reactants, the products, and the rate constant.
  Reaction(const ArrayInt& reactants,
	   const ArrayInt& products,
	   const Number rateConstant) :
    _reactants(reactants),
    _products(products),
    _rateConstant(rateConstant) {
    assert(_reactants.getNull() == 0 && _products.getNull() == 0);
  }

  //! Rebuild from the reactants, the products, and the rate constant.
  void
  rebuild(const ArrayInt& reactants,
	  const ArrayInt& products,
	  const Number rateConstant) {
    _reactants = reactants;
    _products = products;
    assert(_reactants.getNull() == 0 && _products.getNull() == 0);
    _rateConstant = rateConstant;
  }

  //! Copy constructor.
  Reaction(const Reaction& other) :
    _reactants(other._reactants),
    _products(other._products),
    _rateConstant(other._rateConstant)
  {}
  
  //! Assignment operator.
  Reaction&
  operator=(const Reaction& other) {
    if (&other != this) {
      _reactants = other._reactants;
      _products = other._products;
      _rateConstant = other._rateConstant;
    }
    return *this;
  }
  
  //! Destructor.
  ~Reaction()
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{

  //! Get the reactants as a sparse array.
  const ArrayInt&
  getReactants() const {
    return _reactants;
  }

  //! Get the products as a sparse array.
  const ArrayInt&
  getProducts() const {
    return _products;
  }

  //! Get the reactants as a dense array.
  template<typename _Integer>
  void
  getReactants(std::vector<_Integer>* reactants) const {
    std::fill(reactants->begin(), reactants->end(), 0);
    // For each non-zero element.
    for (int i = 0; i != _reactants.size(); ++i) {
#ifdef DEBUG_stochastic_Reaction
      assert(_reactants.getIndex(i) < int(reactants->size()));
#endif
      (*reactants)[_reactants.getIndex(i)] = _reactants[i];
    }
  }

  //! Get the products as a dense array.
  template<typename _Integer>
  void
  getProducts(std::vector<_Integer>* products) const {
    std::fill(products->begin(), products->end(), 0);
    // For each non-zero element.
    for (int i = 0; i != _products.size(); ++i) {
#ifdef DEBUG_stochastic_Reaction
      assert(_products.getIndex(i) < int(products->size()));
#endif
      (*products)[_products.getIndex(i)] = _products[i];
    }
  }

  //! Get the specific reaction probability rate constant.
  Number 
  getRateConstant() const {
    return _rateConstant;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{

  //! Set the reactants with a sparse array.
  void
  setReactants(const ArrayInt& reactants) {
    _reactants = reactants;
  }

  //! Set the products with a sparse array.
  void
  setProducts(const ArrayInt& products) {
    _products = products;
  }

  //! Set the reactants with indices and coefficients.
  template<typename _ForwardIterator>
  void
  setReactants(_ForwardIterator indicesBeginning, 
	       _ForwardIterator indicesEnd,
	       _ForwardIterator coefficientsBeginning, 
	       _ForwardIterator coefficientsEnd) {
    _reactants.rebuild(indicesBeginning, indicesEnd,
		       coefficientsBeginning, coefficientsEnd);
  }

  //! Set the products with indices and coefficients.
  template<typename _ForwardIterator>
  void
  setProducts(_ForwardIterator indicesBeginning, 
	      _ForwardIterator indicesEnd,
	      _ForwardIterator coefficientsBeginning, 
	      _ForwardIterator coefficientsEnd) {
    _products.rebuild(indicesBeginning, indicesEnd,
		      coefficientsBeginning, coefficientsEnd);
  }

  //! Set the specific reaction probability rate constant.
  void
  setRateConstant(const Number rateConstant) {
    _rateConstant = rateConstant;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Functions.
  //@{

  //! Compute the propensity function given the species populations.
  /*!
    Let the reactants be
    \f[
    \sum_{i \in R} b_i x_i
    \f]
    where \f$R\f$ is the set of reactant indices.
    The propensity function is 
    \f[
    a(\mathbf{x}) = c \prod_{i \in R} 
    \left(\begin{array}{c} x_i \\ b_i \end{array} \right)
    \f]
    where c is the rate constant.
  */
  template<typename Container>
  Number
  computePropensityFunction(const Container& populations) const;

  //! Compute the propensity function derivatives given the species populations.
  /*!
    Let the reactants be
    \f[
    \sum_{i \in R} b_i x_i
    \f]
    where \f$R\f$ is the set of reactant indices.
    Note that the derivative of the binomial function is
    \f[
    \frac{\mathrm{d}}{\mathrm{d}x} 
    \left(\begin{array}{c} x \\ b \end{array} \right) = 
    (H_{x} - H_{x - b}) \left(\begin{array}{c} x \\ b \end{array} \right).
    \f]
    Here \f$H_n\f$ is the n_th harmonic number.  
    Thus the derivative of the propensity function with respect to \f$x_i\f$ is
    \f[
    \frac{\partial a}{\partial x_i} = 
    (H_{x_i} - H_{x_i - b_i}) a(\mathbf(x)).
    \f]
  */
  template<typename Container>
  void
  computePropensityFunctionDerivatives(const Container& populations,
				       ArrayInt* derivatives) const;

  //! Compute the state change vector for the reaction.
  void
  computeStateChangeVector(ArrayInt* stateChangeVector) const {
    ads::computeDifference(_products, _reactants, stateChangeVector);
  }

  //@}
};


//--------------------------------------------------------------------------
//! \defgroup stochastic_ReactionFunctions Free functions for Reaction.
//@{


//! Return true if the reaction is valid.
/*! 
  \relates Reaction 

  The species indices must be in the range [0..numberOfSpecies).  The 
  species coefficients must be positive.  The rate constant must be 
  non-negative.
*/
template<typename T>
bool
isValid(const Reaction<T>& reaction, const int numberOfSpecies);

//! Write the reaction in ascii format.
/*! \relates Reaction */
template<typename T>
void
writeAscii(std::ostream& out, const Reaction<T>& x);


//! Read the reaction in ascii format.
/*! \relates Reaction */
template<typename T>
void
readAscii(std::istream& in, Reaction<T>* x);


//! Read the reaction in ascii format.
/*! \relates Reaction */
template<typename T>
void
readReactantsAscii(std::istream& in, Reaction<T>* x);

//! Read the products in ascii format.
/*! \relates Reaction */
template<typename T>
void
readReactantsAscii(std::istream& in, Reaction<T>* x);

//! Write the reaction in ascii format.
/*! \relates Reaction */
template<typename T>
inline
std::ostream&
operator<<(std::ostream& out, const Reaction<T>& x) {
  writeAscii(out, x);
  return out;
}


//! Read the reaction in ascii format.
/*! \relates Reaction */
template<typename T>
inline
std::istream&
operator>>(std::istream& in, Reaction<T>& x) {
  readAscii(in, &x);
  return in;
}


//! Return true if the two reactions are equal.
/*! \relates Reaction */
template<typename T>
inline
bool
operator==(const Reaction<T>& x, const Reaction<T>& y) {
  return x.getReactants() == y.getReactants() &&
    x.getProducts() == y.getProducts() &&
    x.getRateConstant() == y.getRateConstant();
}


//! Return true if the two reactions are equal.
/*! \relates Reaction */
template<typename T>
inline
bool
operator!=(const Reaction<T>& x, const Reaction<T>& y) {
  return !(x == y);
}


//@}

  
END_NAMESPACE_STOCHASTIC

#define __stochastic_Reaction_ipp__
#include "Reaction.ipp"
#undef __stochastic_Reaction_ipp__

#endif
