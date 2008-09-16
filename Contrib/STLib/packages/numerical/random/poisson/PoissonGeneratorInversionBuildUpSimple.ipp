// -*- C++ -*-

#if !defined(__numerical_random_PoissonGeneratorInversionBuildUpSimple_ipp__)
#error This file is an implementation detail of PoissonGeneratorInversionBuildUpSimple.
#endif

#define __numerical_random_PoissonGeneratorInversionMaximumMean_ipp__
#include "PoissonGeneratorInversionMaximumMean.ipp"
#undef __numerical_random_PoissonGeneratorInversionMaximumMean_ipp__

BEGIN_NAMESPACE_NUMERICAL

template<class _Uniform, typename _Result>
inline
typename PoissonGeneratorInversionBuildUpSimple<_Uniform, _Result>::result_type
PoissonGeneratorInversionBuildUpSimple<_Uniform, _Result>::
operator()(const argument_type mean) {
#ifdef DEBUG_numerical_PoissonGeneratorInversionBuildUpSimple
  // If the mean is too large, we will get underflow in computing p.
  // The algorithm will give incorrect results.
  assert(mean < Number(PoissonGeneratorInversionMaximumMean<Number>::Value));
#endif

  // CONTINUE
  const int NumericalFailureBound = 
    2 * PoissonGeneratorInversionMaximumMean<Number>::Value;
  const Number exponential = std::exp(-mean);
  while (true) {
    // Poisson random deviate.
    int deviate = 0;
    // Probability density function.
    Number pdf = exponential;
    // Cumulative distribution function.
    Number cdf = pdf;
    // Uniform random deviate.
    const Number u = transformDiscreteDeviateToContinuousDeviateClosed<Number>
      ((*_discreteUniformGenerator)());
    do {
      if (cdf >= u) {
	return deviate;
      }
      ++deviate;
      pdf *= mean / deviate;
      cdf += pdf;
    } while (deviate != NumericalFailureBound);
  }
}

END_NAMESPACE_NUMERICAL
