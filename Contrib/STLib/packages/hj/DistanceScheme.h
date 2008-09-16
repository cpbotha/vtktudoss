// -*- C++ -*-

/*! 
  \file DistanceScheme.h
  \brief Base class for finite difference schemes for computing distance.
*/

#if !defined(__hj_DistanceScheme_h__)
#define __hj_DistanceScheme_h__

#include "status.h"

#include "../ads/array/Array.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_DistanceScheme)
#define DEBUG_DistanceScheme
#endif

BEGIN_NAMESPACE_HJ

//! Base class for finite differences for computing distance.
/*!
  \param N is the space dimension.
  \param T is the number type.

  The \c Distance class defined the finite difference operations.
  Here we hold copies of the solution array and the status array.
  Classes which derive from \c DistanceScheme will use these operations
  and arrays to implement finite difference schemes.
*/
template<int N, typename T>
class DistanceScheme {
protected:

  //
  // Member data
  //

  //! A reference for the solution array.
  ads::Array<N,T,false> _solution;

  //! A reference for the status array.
  ads::Array<N,Status,false> _status;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  DistanceScheme();
  //! Copy constructor not implemented.
  DistanceScheme(const DistanceScheme&);
  //! Assignment operator not implemented.
  DistanceScheme& 
  operator=(const DistanceScheme&);

public:

  //
  // Constructors
  //

  //! Constructor.
  /*!
    \param solution is the solution array.
    \param status is the status array.
  */
  template<bool A1, bool A2>
  DistanceScheme(ads::Array<N,T,A1>& solution, 
		  ads::Array<N,Status,A2>& status) :
    _solution(solution),
    _status(status)
  {}

  ~DistanceScheme()
  {}

  //
  // File I/O
  //

  //! Write that distance is being computed.
  void 
  put(std::ostream& out) const {
    out << "This is an equation for computing distance.\n";
  }

};

//
// File I/O
//

//! Write to a file stream.
template<int N, typename T>
inline
std::ostream& 
operator<<(std::ostream& out, const DistanceScheme<N,T>& x) {
  x.put(out);
  return out;
}

END_NAMESPACE_HJ

#endif
