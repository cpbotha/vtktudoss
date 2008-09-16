// -*- C++ -*-

/*! 
  \file numerical/random/poisson/PoissonGeneratorStochKit.h
  \brief Uniform random deviates.
*/

#if !defined(__numerical_PoissonGeneratorStochKit_h__)
#define __numerical_PoissonGeneratorStochKit_h__

#include "../normal/Default.h"

#include <algorithm>

// If we are debugging the whole numerical package.
#if defined(DEBUG_numerical) && !defined(DEBUG_PoissonGeneratorStochKit)
#define DEBUG_PoissonGeneratorStochKit
#endif

BEGIN_NAMESPACE_NUMERICAL

//! Generator for Poisson deviates.
/*!
  \image html random/poisson/same/sameStochKitSmallArgument.jpg "Execution times for the same means."
  \image latex random/poisson/same/sameStochKitSmallArgument.pdf "Execution times for the same means." width=0.5\textwidth

  \image html random/poisson/same/sameStochKitLargeArgument.jpg "Execution times for the same means."
  \image latex random/poisson/same/sameStochKitLargeArgument.pdf "Execution times for the same means." width=0.5\textwidth


  \image html random/poisson/different/differentStochKitSmallArgument.jpg "Execution times for different means."
  \image latex random/poisson/different/differentStochKitSmallArgument.pdf "Execution times for different means." width=0.5\textwidth

  \image html random/poisson/different/differentStochKitLargeArgument.jpg "Execution times for different means."
  \image latex random/poisson/different/differentStochKitLargeArgument.pdf "Execution times for different means." width=0.5\textwidth


  \image html random/poisson/distribution/distributionStochKitSmallArgument.jpg "Execution times for a distribution of means."
  \image latex random/poisson/distribution/distributionStochKitSmallArgument.pdf "Execution times for a distribution of means." width=0.5\textwidth

  \image html random/poisson/distribution/distributionStochKitLargeArgument.jpg "Execution times for a distribution of means."
  \image latex random/poisson/distribution/distributionStochKitLargeArgument.pdf "Execution times for a distribution of means." width=0.5\textwidth
*/
template<class _Uniform = DISCRETE_UNIFORM_GENERATOR_DEFAULT,
	 template<class> class _Normal = NORMAL_GENERATOR_DEFAULT,
	 typename _Result = double>
class PoissonGeneratorStochKit {
public:

  //! The number type.
  typedef double Number;
  //! The argument type.
  typedef Number argument_type;
  //! The result type.
  typedef _Result result_type;
  //! The discrete uniform generator.
  typedef _Uniform DiscreteUniformGenerator;
  //! The normal generator.
  typedef _Normal<DiscreteUniformGenerator> NormalGenerator;

private:

  //
  // Member data.
  //

  //! The normal generator.
  NormalGenerator* _normalGenerator;

  //
  // Not implemented.
  //

  //! Default constructor not implemented.
  PoissonGeneratorStochKit();

public:

  //! Construct using the normal generator.
  explicit
  PoissonGeneratorStochKit(NormalGenerator* normalGenerator) :
    _normalGenerator(normalGenerator)
  {}

  //! Copy constructor.
  PoissonGeneratorStochKit(const PoissonGeneratorStochKit& other) :
    _normalGenerator(other._normalGenerator)
  {}

  //! Assignment operator.
  PoissonGeneratorStochKit&
  operator=(const PoissonGeneratorStochKit& other) {
    if (this != &other) {
      _normalGenerator = other._normalGenerator;
    }
    return *this;
  }

  //! Destructor.
  ~PoissonGeneratorStochKit()
  {}

  //! Seed the uniform random number generator.
  void
  seed(const typename DiscreteUniformGenerator::result_type seedValue) {
    _normalGenerator->seed(seedValue);
  }

  //! Return a Poisson deviate with the specifed mean.
  result_type
  operator()(argument_type mean);

private:

  int
  ignpoi(float mu);

  float
  snorm();

  float
  sexpo();

  float
  fsign(float num, float sign);

};

END_NAMESPACE_NUMERICAL

#define __numerical_random_PoissonGeneratorStochKit_ipp__
#include "PoissonGeneratorStochKit.ipp"
#undef __numerical_random_PoissonGeneratorStochKit_ipp__

#endif
