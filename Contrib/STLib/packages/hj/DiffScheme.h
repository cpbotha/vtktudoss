// -*- C++ -*-

/*! 
  \file DiffScheme.h
  \brief Finite difference operations for an N-D grid.
*/

#if !defined(__hj_DiffScheme_h__)
#define __hj_DiffScheme_h__

#include "status.h"

#include "../ads/priority_queue/PriorityQueueBinaryHeapArray.h"
#include "../ads/priority_queue/PriorityQueueBinaryHeapStoreKeys.h"
#include "../ads/priority_queue/PriorityQueueCellArray.h"

#include <deque>

// If we are debugging the whole hj namespace.
#if defined(DEBUG_hj) && !defined(DEBUG_DiffScheme)
#define DEBUG_DiffScheme
#endif

BEGIN_NAMESPACE_HJ

//! Finite difference scheme.
/*!
  \param N The space dimension.
  \param T The number type.
  \param Equation Represents the equation to be solved.

  The DiffScheme class implements the common functionality of difference
  schemes.  It stores the status array and the equation to be solved.
  Here we define the \c label() functions for the MCC algorithm and the 
  FM method.
*/
template<int N, typename T, class Equation>
class DiffScheme {
 public:

  //! The number type.
  typedef T Number;

  //! An array reference of numbers.
  typedef ads::Array<N,Number,false> NumberArrayRef;
  //! An array of status variables.
  typedef ads::Array<N,Status> StatusArray;

  //! A multi-index.
  typedef typename StatusArray::index_type Index;
  //! An index range.
  typedef typename StatusArray::range_type Range;

 protected:

  // 
  // Member data
  //

  //! The index ranges of the grid, not counting ghost points.
  Range _index_ranges;

  //! A reference for the solution array.
  NumberArrayRef _solution;

  //! The status array.
  StatusArray _status;

  //! The equation.
  Equation _equation;

 private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  DiffScheme();
  //! Copy constructor not implemented.
  DiffScheme(const DiffScheme&);
  //! Assignment operator not implemented.
  DiffScheme& 
  operator=(const DiffScheme&);

 public:

  //
  // Constructors
  //

  //! Constructor.
  /*!
    \param index_ranges are the index ranges of the real grid points 
    (not counting the ghost grid points, if there are any).
    \param solution is the solution array.
    \param dx is the grid spacing.
    \param is_concurrent indicates whether we are using the concurrent 
    algorithm.
  */
  template<bool A>
  DiffScheme(const Range& index_ranges, ads::Array<N,Number,A>& solution, 
	      Number dx, bool is_concurrent);

  ~DiffScheme()
  {}

  //
  // Mathematical member functions.
  //

  //! Initialize the status array.
  /*!
    Set the interior status points to \c UNLABELED.
    Set the boundary status points to \c VOID.
  */
  void 
  initialize();

  //! Set the status at the specified grid point to \c INITIAL.
  void 
  set_initial(const Index& i) {
    _status(i) = INITIAL;
  }

  //! Set the status at the specified grid point to \c KNOWN.
  void 
  set_known(const Index& i) {
    _status(i) = KNOWN;
  }

  //! Label the specified grid point with the \c value.
  /*!
    If previously unlabeled, add to \c unlabeled_neighbors.
  */
  template<class Container>
  void 
  label(Container& unlabeled_neighbors, const Index& i,
	 Number value);

  //! Label the specified grid point with the \c value.
  /*!
    If previously unlabeled, add to labeled.  If the distance decreases, 
    adjust the position in the \c labeled heap.
  */
  template<typename HandleType>
  void 
  label(ads::PriorityQueueBinaryHeapArray<HandleType>& labeled, 
	 const Index& i, Number value);

  //! Label the specified grid point with the \c value.
  /*!
    If previously unlabeled, or if the solution decreases, add to \c labeled.
  */
  template<typename HandleType>
  void 
  label(ads::PriorityQueueBinaryHeapStoreKeys<HandleType>& labeled, 
	 const Index& i, Number value);

  //! Label the specified grid point with the \c value.
  /*!
    If previously unlabeled, or if the solution decreases, add to \c labeled.
  */
  template<typename HandleType>
  void 
  label(ads::PriorityQueueCellArray<HandleType>& labeled, 
	 const Index& i, Number value);

  //! Label the specified grid point with the \c value.
  /*!
    The sorted grid points method uses this function.
  */
  void 
  label(int, const Index& i, Number value);

  //! Return a lower bound on correct values.
  /*!
    This is just a wrapper for the \c lower_bound() function defined in
    the equation.
  */
  Number
  lower_bound(const Index& i, const Number min_unknown) const { 
    return _equation.lower_bound(i, min_unknown);
  }

  //
  // Accesors
  //

  //! Return the radius of the stencil.
  static
  int
  radius() {
    return Equation::radius();
  }

  //! Return true if the status of the grid point is \c KNOWN.
  bool 
  is_known(const Index& i) const { 
    return _status(i) == KNOWN; 
  }

  //! Return true if the status of the grid point is \c INITIAL.
  bool 
  is_initial(const Index& i) const { 
    return _status(i) == INITIAL; 
  }

  //! Return true if the status of the grid point is \c LABELED.
  bool 
  is_labeled(const Index& i) const { 
    return _status(i) == LABELED;
  }

  //! Return true if the status of the grid point is \c UNLABELED.
  bool 
  is_unlabeled(const Index& i) const { 
    return _status(i) == UNLABELED;
  }

  //! Return the status array.
  const StatusArray& 
  status() const { 
    return _status; 
  }
  
  //! Return a const reference to the equation.
  const Equation&
  equation() const {
    return _equation;
  }

  //
  // Manipulators
  //

  //! Return the status array.
  StatusArray& 
  status() { 
    return _status; 
  }
  
  //! Return a reference to the equation.
  Equation&
  equation() {
    return _equation;
  }

  //
  // File I/O
  //

  //! Write statistics for the status array.
  void
  print_statistics(std::ostream& out) const;

  //! Write the status array and data associated with the equation.
  void 
  put(std::ostream& out) const;

 protected:

  //! Return true if the grid point is labeled or unlabeled.
  bool
  is_labeled_or_unlabeled(const Index& i) const {
    return _status(i) == LABELED || _status(i) == UNLABELED;
  }

};

//
// File I/O
//

//! Write to a file stream.
template<int N, typename T, class Equation>
std::ostream& 
operator<<(std::ostream& out, const DiffScheme<N,T,Equation>& x);

END_NAMESPACE_HJ

#define __hj_DiffScheme_ipp__
#include "DiffScheme.ipp"
#undef __hj_DiffScheme_ipp__

#endif
