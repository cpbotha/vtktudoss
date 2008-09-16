// -*- C++ -*-

#if !defined(__numerical_random_PoissonGeneratorExponentialInterArrival_ipp__)
#error This file is an implementation detail of PoissonGeneratorExponentialInterArrival.
#endif

BEGIN_NAMESPACE_NUMERICAL

template<class _Uniform, template<class> class _Exponential, typename _Result>
inline
typename PoissonGeneratorExponentialInterArrival<_Uniform, _Exponential, _Result>::result_type
PoissonGeneratorExponentialInterArrival<_Uniform, _Exponential, _Result>::
operator()(argument_type mean) {
  result_type deviate = -1;
  do {
    ++deviate;
    mean -= (*_exponentialGenerator)();
  } while (mean > 0);
  return deviate;
}

END_NAMESPACE_NUMERICAL
