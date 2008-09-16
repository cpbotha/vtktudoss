// -*- C++ -*-

/*! 
  \file shortest_paths/UniformRandom.h
  \brief Implements a class for a random number generator.

  The random number generators return numbers in a range with a given ratio
  between the max and min.
*/

#if !defined(__UniformRandom_h__)
#define __UniformRandom_h__

#include "../ads/functor/UniformRandom.h"

#include <cassert>

//! A random number generator.
template <typename NumberType>
class UniformRandom
{
private:

  //
  // Member data.
  //

  //! The random number generator.
  ads::SubtractiveRNG _random;

  //! The min number.
  NumberType _min;

  //! The range.
  NumberType _range;

  //
  // Not implemented.
  //
  
  //! Default constructor not implemented.
  UniformRandom();

  //! Copy constructor not implemented.
  UniformRandom( const UniformRandom& );

  //! Assignment operator not implemented.
  UniformRandom& 
  operator=( const UniformRandom& );

public:

  //
  // Constructors, Destructor.
  //
  
  //! Construct from a ratio.
  /*!
    If \param ratio != 0, the random numbers will be in the range 
    [1 / ratio .. 1].  \param ratio = 0 signifies an infinite ratio.  
    The random numbers will be in the range [0 .. 1].
  */
  explicit 
  UniformRandom( NumberType ratio ) :
    _random()
  {
    assert( ratio >= 0 );
    if ( ratio == 0 ) {
      _min = 0;
      _range = 1;
    }
    else {
      _min = 1 / ratio;
      _range = (1 - _min);
    }
  }

  //! Trivial destructor.
  virtual 
  ~UniformRandom()
  {}

  //
  // Functional
  //

  //! Return a random number.
  NumberType operator()()
  { 
    return _min + _range * (_random(1024) / 1023.0);
  }

};

//! A random integer generator.
template <>
class UniformRandom<int>
{
private:

  //
  // Member data.
  //

  //! The random number generator.
  ads::SubtractiveRNG _random;

  //! The min number.
  int _min;

  //! The range.
  int _range;

  //
  // Not implemented.
  //
  
  //! Default constructor not implemented.
  UniformRandom();

  //! Copy constructor not implemented.
  UniformRandom( const UniformRandom& );

  //! Assignment operator not implemented.
  const UniformRandom& operator=( const UniformRandom& );

public:

  //
  // Constructors, Destructor.
  //
  
  //! Construct from a ratio.
  /*!
    If \param ratio != 0, the random numbers will be in the range [1 .. ratio].
    \param ratio = 0 signifies an infinite ratio.  The random numbers will 
    be in the range [0 .. 1000].
  */
  explicit UniformRandom( int ratio ) :
    _random()
  {
    assert( ratio >= 0 );
    if ( ratio == 0 ) {
      _min = 0;
      _range = 1001;
    }
    else {
      _min = 1;
      _range = ratio;
    }
  }

  //! Trivial destructor.
  virtual ~UniformRandom()
  {}

  //
  // Functional
  //

  //! Return a random number.
  int operator()()
  { 
    return _min + _random( _range );
  }

};

#endif
