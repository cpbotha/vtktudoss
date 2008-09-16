// -*- C++ -*-

/*! 
  \file EikonalScheme.h
  \brief Base class for finite difference schemes for the eikonal equation.
*/

#if !defined(__hj_EikonalScheme_h__)
#define __hj_EikonalScheme_h__

#include "status.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_EikonalScheme)
#define DEBUG_EikonalScheme
#endif

BEGIN_NAMESPACE_HJ

//! Base class for finite differences for the eikonal equation.
/*!
  \param N is the space dimension.
  \param T is the number type.

  The \c Eikonal class defined the finite difference operations.
  Here we hold copies of the solution array and the status array and
  store the speed array.  Classes which derive from \c EikonalScheme
  will use these operations and arrays to implement finite
  difference schemes.
*/
template<int N, typename T>
class EikonalScheme {
protected:

  //
  // Member data
  //

  //! A reference for the solution array.
  ads::Array<N,T,false> _solution;

  //! A reference for the status array.
  ads::Array<N,Status,false> _status;

  //! The inverse speeed array.
  ads::Array<N,T> _inverseSpeed;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  EikonalScheme();
  //! Copy constructor not implemented.
  EikonalScheme(const EikonalScheme&);
  //! Assignment operator not implemented.
  EikonalScheme& 
  operator=(const EikonalScheme&);

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
  EikonalScheme(ads::Array<N,T,A1>& solution, 
		ads::Array<N,Status,A2>& status) :
    _solution(solution),
    _status(status),
    _inverseSpeed(solution.ranges())
  {}

  ~EikonalScheme()
  {}

  //
  // Manipulators.
  //

  //! Return a reference to the inverse speed array.
  ads::Array<N,T>&
  getInverseSpeed() {
    return _inverseSpeed;
  }

  //
  // File I/O
  //

  //! Write the inverse speed array.
  void 
  put(std::ostream& out) const;
};

//
// File I/O
//

//! Write to a file stream.
template<int N, typename T>
inline
std::ostream& 
operator<<(std::ostream& out, const EikonalScheme<N,T>& x) {
  x.put(out);
  return out;
}

END_NAMESPACE_HJ

#define __hj_EikonalScheme_ipp__
#include "EikonalScheme.ipp"
#undef __hj_EikonalScheme_ipp__

#endif
