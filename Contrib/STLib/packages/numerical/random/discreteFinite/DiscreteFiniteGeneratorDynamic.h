// -*- C++ -*-

/*! 
  \file numerical/random/discreteFinite/DiscreteFiniteGeneratorDynamic.h
  \brief Discrete, finite deviate. Dynamic data structure.
*/

#if !defined(__numerical_DiscreteFiniteGeneratorDynamic_h__)
#define __numerical_DiscreteFiniteGeneratorDynamic_h__

#include "linearSearch.h"

#include "../uniform/ContinuousUniformGenerator.h"
#include "../../../ext/vector.h"

#include <numeric>

BEGIN_NAMESPACE_NUMERICAL

//! Discrete, finite deviate. Dynamic data structure.
/*!
  \param Generator is the discrete, uniform generator.
  \param Number is the floating point number type. By default it is double.

  The probabilities vary in time. They are integrated to form cumulative
  probabilities for each event. Also, new events may be added. The deviate
  is determined with a linear search on the probabilities. When a deviate is 
  drawn, the probabilities are reset to zero. Events may be deleted only 
  after drawing a deviate and before integrating any probabilities for the 
  next event. 
*/
template<class _Generator = DISCRETE_UNIFORM_GENERATOR_DEFAULT,
	 typename _Number = double>
class DiscreteFiniteGeneratorDynamic {
  //
  // Public types.
  //
public:

  //! The discrete uniform generator.
  typedef _Generator DiscreteUniformGenerator;
  //! The continuous uniform generator.
  typedef ContinuousUniformGeneratorClosed<DiscreteUniformGenerator>
  ContinuousUniformGenerator;
  //! The number type.
  typedef _Number Number;
  //! The argument type.
  typedef void argument_type;
  //! The result type.
  typedef int result_type;

  //
  // Member data.
  //
private:

  //! The continuous uniform generator.
  ContinuousUniformGenerator _continuousUniformGenerator;
  //! The event indices.
  std::vector<int> _events;
  //! Probability mass function.  (This is scaled and may not sum to unity.)
  std::vector<Number> _pmf;
  //! The sum of the PMF.
  Number _pmfSum;

  //
  // Not implemented.
  //
private:

  //! Default constructor not implemented.
  DiscreteFiniteGeneratorDynamic();

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct using the uniform generator.
  explicit
  DiscreteFiniteGeneratorDynamic(DiscreteUniformGenerator* generator) :
    // Make a continuous uniform generator using the discrete uniform generator.
    _continuousUniformGenerator(generator),
    _events(),
    _pmf(),
    _pmfSum(0)
  {}

  //! Copy constructor.
  /*!
    \note The discrete, uniform generator is not copied.  Only the pointer
    to it is copied.
  */
  DiscreteFiniteGeneratorDynamic(const DiscreteFiniteGeneratorDynamic& other) :
    _continuousUniformGenerator(other._continuousUniformGenerator),
    _events(other._events),
    _pmf(other._pmf),
    _pmfSum(other._pmfSum)
  {}

  //! Assignment operator.
  /*!
    \note The discrete, uniform generator is not copied.  Only the pointer
    to it is copied.
  */
  DiscreteFiniteGeneratorDynamic&
  operator=(const DiscreteFiniteGeneratorDynamic& other) {
    if (this != &other) {
      _continuousUniformGenerator = other._continuousUniformGenerator;
      _events = other._events;
      _pmf = other._pmf;
      _pmfSum = other._pmfSum;
    }
    return *this;
  }

  //! Destructor.
  /*! The memory for the discrete, uniform generator is not freed. */
  ~DiscreteFiniteGeneratorDynamic()
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Random number generation.
  //@{
public:

  //! Seed the uniform random number generator.
  void
  seed(const typename DiscreteUniformGenerator::result_type seedValue) {
    _continuousUniformGenerator.seed(seedValue);
  }

#if 0
  //! Return true if any events have non-zero probability.
  bool
  isEvent() {
    _pmfSum = sum(_pmf);
    return _pmfSum != 0;
  }
#endif

  //! Return a discrete, finite deviate. Reset the probabilities to zero.
  result_type
  operator()() {
    _pmfSum = sum(_pmf);
#ifdef DEBUG_stlib
    assert(_pmfSum != 0);
#endif
    // Draw the deviate.
    int index;
    do {
      index = linearSearchChopDownUnguarded
	(_pmf.begin(), _pmf.end(), _continuousUniformGenerator() * _pmfSum);
    } while (_pmf[index] == 0);
    // Reset the probabilities.
    std::fill(_pmf.begin(), _pmf.end(), 0);
    // Return the event.
    return _events[index];
  }

  //! Clear the events.
  void
  clear() {
    _events.clear();
    _pmf.clear();
  }

  //! Insert a new event.
  void
  insert(const int index) {
    _events.push_back(index);
    _pmf.push_back(0);
  }

  //! Erase an event.
  /*! \pre All of the probabilities must be zero. */
  void
  erase(const int index) {
    const std::ptrdiff_t i = std::find(_events.begin(), _events.end(), index)
      - _events.begin();
#ifdef DEBUG_stlib
    assert(std::size_t(i) != _events.size());
#endif
    _events[i] = _events.back();
#ifdef DEBUG_stlib
    assert(_pmf[i] == 0);
#endif
    _events.pop_back();
    _pmf.pop_back();
  }

  //! Add the probabilities contribution for the specified time interval.
  void
  addPmf(const std::vector<Number>& probabilities, const Number delta) {
    for (std::size_t n = 0; n != _events.size(); ++n) {
      _pmf[n] += delta * probabilities[_events[n]];
    }
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Return the container of events.
  const std::vector<int>&
  getEvents() const {
    return _events;
  }

  //! Get the PMF for the specified event.
  Number
  getPmf(const int event) const {
    return _pmf[event];
  }

  //! Return true if the specified event has the maximum probability.
  bool
  isMaximum(const int event) const {
    const Number p = _pmf[event];
    for (std::size_t i = 0; i != _pmf.size(); ++i) {
      if (_pmf[i] > p) {
	return false;
      }
    }
    return true;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
public:

  bool
  operator==(const DiscreteFiniteGeneratorDynamic& other) const {
    return _events == other._events && _pmf == other._pmf &&
      _pmfSum == other._pmfSum;
  }

  //@}
};

END_NAMESPACE_NUMERICAL

#endif
