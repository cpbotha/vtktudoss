// -*- C++ -*-

/*! 
  \file numerical/random/normal/NormalGeneratorZigguratGsl.h
  \brief Normal random deviate with zero mean and unit variance.
*/

#if !defined(__numerical_NormalGeneratorZigguratGsl_h__)
#define __numerical_NormalGeneratorZigguratGsl_h__

#include "../uniform/DiscreteUniformGeneratorMt19937Gsl.h"

#include <gsl/gsl_randist.h>

#include <cmath>

// If we are debugging the whole numerical package.
#if defined(DEBUG_numerical) && !defined(DEBUG_NormalGeneratorZigguratGsl)
#define DEBUG_NormalGeneratorZigguratGsl
#endif

BEGIN_NAMESPACE_NUMERICAL

template<class _Generator = DiscreteUniformGeneratorMt19937Gsl>
class NormalGeneratorZigguratGsl;

//! Normal random deviate with zero mean and unit variance.
/*!
  This is a template specialization for the case that the uniform deviate
  generator is UniformGeneratorMt19937Gsl .
  This generator can be initialized in the constructor or with seed().

  This functor wraps the function call to the Marsaglia-Tsang ziggurat 
  algorithm for computing Gaussian deviates in the
  <a href="http://www.gnu.org/software/gsl/">GNU Scientific Library</a>.
  It returns a floating point value that is a random deviate drawn from a 
  normal (Gaussian) distribution with zero mean and unit variance.
*/
template<>
class NormalGeneratorZigguratGsl<DiscreteUniformGeneratorMt19937Gsl> {
public:

  //! The discrete uniform generator.
  typedef DiscreteUniformGeneratorMt19937Gsl DiscreteUniformGenerator;
  //! The number type.
  typedef double Number;
  //! The argument type.
  typedef void argument_type;
  //! The result type.
  typedef Number result_type;

  //
  // Member data.
  //

private:

  //! The discrete uniform generator.
  DiscreteUniformGenerator* _discreteUniformGenerator;

  //
  // Not implemented.
  //

  //! Default constructor not implemented.
  NormalGeneratorZigguratGsl();

public:

  //! Construct using the uniform generator.
  explicit
  NormalGeneratorZigguratGsl(DiscreteUniformGenerator* generator) :
    _discreteUniformGenerator(generator) 
  {}

  //! Copy constructor.
  /*!
    \note The discrete, uniform generator is not copied.  Only the pointer
    to it is copied.
  */
  NormalGeneratorZigguratGsl(const NormalGeneratorZigguratGsl& other) :
    _discreteUniformGenerator(other._discreteUniformGenerator)
  {}

  //! Assignment operator.
  /*!
    \note The discrete,uniform generator is not copied.  Only the pointer
    to it is copied.
  */
  NormalGeneratorZigguratGsl&
  operator=(const NormalGeneratorZigguratGsl& other) {
    if (this != &other) {
      _discreteUniformGenerator = other._discreteUniformGenerator;
    }
    return *this;
  }

  //! Destructor.
  /*!
    The memory for the discrete, uniform generator is not freed.
  */
  ~NormalGeneratorZigguratGsl()
  {}

  //! Seed the uniform random number generator.
  void
  seed(const DiscreteUniformGenerator::result_type seedValue) {
    _discreteUniformGenerator->seed(seedValue);
  }

  //! Return a standard normal deviate.
  result_type
  operator()() {
    return gsl_ran_gaussian_ziggurat(_discreteUniformGenerator->getGenerator(),
				     1);
  }

  //! Return a normal deviate with specified mean and variance.
  result_type
  operator()(const Number mean, const Number variance) {
    return std::sqrt(variance) * operator()() + mean;
  }

  //! Get the discrete uniform generator.
  DiscreteUniformGenerator*
  getDiscreteUniformGenerator() {
    return _discreteUniformGenerator;
  }
};


END_NAMESPACE_NUMERICAL

#endif
