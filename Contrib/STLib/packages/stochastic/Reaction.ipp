// -*- C++ -*-

#if !defined(__stochastic_Reaction_ipp__)
#error This file is an implementation detail of Reaction.
#endif

BEGIN_NAMESPACE_STOCHASTIC


// Rebuild from the reactants, the products, and the rate constant.
template<typename T>
template<typename IntegerForwardIterator>
inline
void
Reaction<T>::
rebuild(IntegerForwardIterator reactantIndicesBeginning,
	IntegerForwardIterator reactantIndicesEnd,
	IntegerForwardIterator reactantCoefficientsBeginning,
	IntegerForwardIterator reactantCoefficientsEnd,
	IntegerForwardIterator productIndicesBeginning,
	IntegerForwardIterator productIndicesEnd,
	IntegerForwardIterator productCoefficientsBeginning,
	IntegerForwardIterator productCoefficientsEnd,
	const Number rateConstant) {
  _reactants.rebuild(reactantIndicesBeginning, reactantIndicesEnd,
		     reactantCoefficientsBeginning, reactantCoefficientsEnd,
		     0);
  _products.rebuild(productIndicesBeginning, productIndicesEnd,
		    productCoefficientsBeginning, productCoefficientsEnd, 0);
  _rateConstant = rateConstant;
}



// 3.7
// Compute the binomial coefficient for the case that n >= 0 and k > 0.
// This is a little faster than the version in the numerical package.
template<typename _Result, typename _Argument>
inline
_Result
_computeBinomialCoefficient(const _Argument n, const int k) {
#ifdef DEBUG_stochastic_Reaction
  assert(n >= 0 && k > 0);
#endif
  _Result result = n;
  for (int i = 1; i != k; ++i) {
    result *= (n - i);
    result /= (i + 1);
  }
  return result;
}


// Compute the propensity function given the the species populations.
template<typename T>
template<typename Container>
inline
typename Reaction<T>::Number
Reaction<T>::
computePropensityFunction(const Container& populations) const {
  Number result = _rateConstant;

  int speciesIndex;
  for (int i = 0; i != _reactants.size(); ++i) {
    speciesIndex = _reactants.getIndices()[i];
#ifdef DEBUG_stochastic_Reaction
    assert(0 <= speciesIndex && speciesIndex < int(populations.size()));
#endif
    // If a single reaction firing would make the population negative.
    if (populations[speciesIndex] < _reactants[i]) {
      // The reaction cannot occur.
      return 0;
    }
    result *= _computeBinomialCoefficient<Number>
      (populations[speciesIndex], _reactants[i]);
  }
#ifdef DEBUG_stochastic_Reaction
  assert(result >= 0);
#endif
  return result;
}


// Compute the propensity function derivatives given the the species 
// populations.
template<typename T>
template<typename Container>
inline
void
Reaction<T>::
computePropensityFunctionDerivatives(const Container& populations,
				     ArrayInt* derivatives) const {
  // The value of the propensity function is used the evaluation of its 
  // derivatives.
  const T propensityFunction = computePropensityFunction(populations);

  // Compute the derivatives.
  std::vector<int> indices;
  std::vector<T> values;
  T value;
  int speciesIndex;
  typename Container::value_type population;
  for (int i = 0; i != _reactants.size(); ++i) {
    speciesIndex = _reactants.getIndices()[i];
    population = populations[speciesIndex];
    value = numerical::computeDifferenceOfHarmonicNumbers<T>
      (population, population - _reactants[i]) * propensityFunction;
    // Record the non-zero values.
    if (value != 0) {
      indices.push_back(speciesIndex);
      values.push_back(value);
    }
  }

  // CONTINUE: This is not efficient.
  // Rebuild the sparse array from the vectors.
  derivatives->rebuild(indices.begin(), indices.end(), 
		       values.begin(), values.end(), 0);
}

//--------------------------------------------------------------------------
// Accessors.
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Free functions.
//--------------------------------------------------------------------------

// Return true if the reaction is valid.
template<typename T>
inline
bool
isValid(const Reaction<T>& reaction, const int numberOfSpecies) {
  // The reactants and products may not both be empty.
  if (reaction.getReactants().size() == 0 &&
      reaction.getProducts().size() == 0) {
    return false;
  }

  // Check the reactant indices and coefficients.
  for (int i = 0; i != reaction.getReactants().size(); ++i) {
    if (reaction.getReactants().getIndices()[i] < 0 || 
	reaction.getReactants().getIndices()[i] >= numberOfSpecies) {
      return false;
    }
    if (reaction.getReactants()[i] <= 0) {
      return false;
    }
  }

  // Check the product indices and coefficients.
  for (int i = 0; i != reaction.getProducts().size(); ++i) {
    if (reaction.getProducts().getIndices()[i] < 0 || 
	reaction.getProducts().getIndices()[i] >= numberOfSpecies) {
      return false;
    }
    if (reaction.getProducts()[i] <= 0) {
      return false;
    }
  }

  // Check the rate constant.
  if (reaction.getRateConstant() < 0) {
    return false;
  }

  return true;
}
  
// Write the reaction in ascii format.
template<typename T>
inline
void
writeAscii(std::ostream& out, const Reaction<T>& x) {
  out << x.getReactants() << x.getProducts() << x.getRateConstant() << "\n";
}


// Read the reaction in ascii format.
template<typename T>
inline
void
readAscii(std::istream& in, Reaction<T>* x) {
  typename Reaction<T>::ArrayInt reactants, products;
  T rateConstant;
  in >> reactants >> products >> rateConstant;
  x->rebuild(reactants, products, rateConstant);
}


// Read the reaction in ascii format.
template<typename T>
inline
void
readReactantsAscii(std::istream& in, Reaction<T>* x) {
  int size = -1;
  in >> size;
#ifdef DEBUG_stochastic
  assert(size >= 0);
#endif
  std::vector<int> indices(size), coefficients(size);
  for (int i = 0; i != size; ++i) {
    in >> indices[i] >> coefficients[i];
  }
  x->setReactants(indices.begin(), indices.end(), coefficients.begin(),
		  coefficients.end());
}


// Read the products in ascii format.
template<typename T>
inline
void
readProductsAscii(std::istream& in, Reaction<T>* x) {
  int size = -1;
  in >> size;
#ifdef DEBUG_stochastic
  assert(size >= 0);
#endif
  std::vector<int> indices(size), coefficients(size);
  for (int i = 0; i != size; ++i) {
    in >> indices[i] >> coefficients[i];
  }
  x->setProducts(indices.begin(), indices.end(), coefficients.begin(),
		 coefficients.end());
}


END_NAMESPACE_STOCHASTIC

// End of file.
