// -*- C++ -*-

/*! 
  \file stochastic/ReactionPriorityQueue.h
  \brief A reaction.
*/

#if !defined(__stochastic_ReactionPriorityQueue_h__)
#define __stochastic_ReactionPriorityQueue_h__

#include "defs.h"

#include "../ads/indexedPriorityQueue/IndexedPriorityQueueLinearSearch.h"
#include "../numerical/random/exponential/Default.h"
#include "../third-party/loki/TypeManip.h"

#include <vector>

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_ReactionPriorityQueue)
#define DEBUG_stochastic_ReactionPriorityQueue
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! The reaction priority queue.
/*!
  CONTINUE.
*/
template<class _IndexedPriorityQueue =
	 ads::IndexedPriorityQueueLinearSearch<>,
	 class _ExponentialGenerator = numerical::EXPONENTIAL_GENERATOR_DEFAULT<>,
	 bool _UseInfinity = false>
class ReactionPriorityQueue :
  private _IndexedPriorityQueue {
  //
  // Private types.
  //
private:

  typedef _IndexedPriorityQueue Base;

  //
  // Enumerations.
  //
public:

  //! Whether to use the infinity value defined by IEEE 754.
  enum {UseInfinity = _UseInfinity, UsesPropensities = Base::UsesPropensities};

  //
  // Public types.
  //
public:

  //! The number type.
  typedef typename Base::Key Number;
  //! The discrete uniform random number generator.
  typedef typename _ExponentialGenerator::DiscreteUniformGenerator 
  DiscreteUniformGenerator;

  //
  // Member data.
  //
private:

  _ExponentialGenerator _exponential;

  //
  // Not implemented.
  //
private:

  //! Default constructor not implemented.
  ReactionPriorityQueue();
  //! Copy constructor not implemented.
  ReactionPriorityQueue(const ReactionPriorityQueue&);
  //! Assignment operator not implemented.
  const ReactionPriorityQueue&
  operator=(const ReactionPriorityQueue&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the size and the uniform generator.
  ReactionPriorityQueue(const int size, DiscreteUniformGenerator* uniform) :
    Base(size),
    _exponential(uniform)
  {}

  //! Construct.  The base class uses hashing.
  ReactionPriorityQueue(const int size,
			DiscreteUniformGenerator* uniform,
			const int hashTableSize,
			const Number targetLoad) :
    Base(size, hashTableSize, targetLoad),
    _exponential(uniform)
  {}

  //! Store a pointer to the propensities in the indexed priority queue.
  void
  setPropensities(const std::vector<Number>* propensities) {
    Base::setPropensities(propensities);
  }

  //! Destructor.
  ~ReactionPriorityQueue() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Return the key of the specified element.
  using Base::get;

  //! Return the index of the top element.
  using Base::top;

  //! Return the value that signifies infinity.
  static
  Number
  infinity() {
    return infinity(Loki::Int2Type<UseInfinity>());
  }

  //! Return a const reference to the discrete, uniform generator.
  const DiscreteUniformGenerator&
  getDiscreteUniformGenerator() const {
    return *_exponential.getDiscreteUniformGenerator();
  }

private:

  //! Return the value that signifies infinity.
  static
  Number
  infinity(Loki::Int2Type<false> /*dummy*/) {
    return std::numeric_limits<Number>::max();
  }

  //! Return the value that signifies infinity.
  static
  Number
  infinity(Loki::Int2Type<true> /*dummy*/) {
    return std::numeric_limits<Number>::infinity();
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Push the value into the queue.
  void
  push(const int index, const Number time, const Number propensity) {
    push(index, time, propensity, Loki::Int2Type<UseInfinity>());
  }

  //! Push the top value into the queue.
  void
  pushTop(const Number time, const Number propensity) {
    pushTop(time, propensity, Loki::Int2Type<UseInfinity>());
  }

  //! Update the value in the queue.
  void
  update(const int index, const Number time, const Number oldPropensity,
	 const Number newPropensity) {
    update(index, time, oldPropensity, newPropensity, 
	   Loki::Int2Type<UseInfinity>());
  }

  //! Clear the queue.
  using Base::clear;

  //! Set the constant used to balance costs.
  void
  setCostConstant(const Number costConstant) {
    Base::setCostConstant(costConstant);
  }

private:

  //! Push the value into the queue.
  void
  push(const int index, const Number time, const Number propensity,
       Loki::Int2Type<false> /*dummy*/) {
#ifdef DEBUG_stochastic
    assert(propensity != 0);
#endif
    Base::push(index, time + _exponential() / propensity);
  }

  //! Push the value into the queue.
  void
  push(const int index, const Number time, const Number propensity,
       Loki::Int2Type<true> /*dummy*/) {
    Base::push(index, time + _exponential() / propensity);
  }

  //! Push the top value into the queue.
  void
  pushTop(const Number time, const Number propensity,
	  Loki::Int2Type<false> /*dummy*/) {
    if (propensity != 0) {
      Base::pushTop(time + _exponential() / propensity);
    }
    else {
      Base::popTop();
    }
  }

  //! Push the top value into the queue.
  void
  pushTop(const Number time, const Number propensity,
	  Loki::Int2Type<true> /*dummy*/) {
    Base::pushTop(time + _exponential() / propensity);
  }

  //! Update the value in the queue.
  void
  update(const int index, const Number time, const Number oldPropensity,
	 const Number newPropensity, Loki::Int2Type<false> /*dummy*/) {
    if (oldPropensity != newPropensity) {
      if (oldPropensity == 0) {
	push(index, time, newPropensity);
      }
      else if (newPropensity == 0) {
	Base::pop(index);
      }
      else {
	Base::set(index, time + (oldPropensity / newPropensity) * 
		  (Base::get(index) - time));
      }
    }
  }

  //! Update the value in the queue.
  void
  update(const int index, const Number time, const Number oldPropensity,
	 const Number newPropensity, Loki::Int2Type<true> /*dummy*/) {
    if (oldPropensity == 0) {
      push(index, time, newPropensity);
    }
    else {
      Base::set(index, time + (oldPropensity / newPropensity) * 
		(Base::get(index) - time));
    }
  }

  //@}
};

END_NAMESPACE_STOCHASTIC

#endif
