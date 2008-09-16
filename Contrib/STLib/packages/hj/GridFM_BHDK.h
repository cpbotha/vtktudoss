// -*- C++ -*-

/*! 
  \file GridFM_BHDK.h
  \brief An N-D grid.  Fast Marching method.  Binary Heap with Dynamic Keys.
*/

#if !defined(__hj_GridFM_BHDK_h__)
#define __hj_GridFM_BHDK_h__

#include "Grid.h"

#include "../ads/priority_queue/PriorityQueueBinaryHeapArray.h"

#include <vector>
#include <iterator>
#include <algorithm>

// If we are debugging the whole hj namespace.
#if defined(DEBUG_hj) && !defined(DEBUG_GridFM_BHDK)
#define DEBUG_GridFM_BHDK
#endif

BEGIN_NAMESPACE_HJ

//! The fast marching method for solving static H-J equations.
/*!
  This class implements the fast marching method for solving static 
  Hamilton-Jacobi equations in the \c solve() member function.  It uses
  a binary heap with dynamic keys.
*/
template <int N, typename T, class DifferenceScheme>
class GridFM_BHDK :
  public Grid<N,T,DifferenceScheme>
{
private:

  typedef Grid<N,T,DifferenceScheme> Base;
  //! The index type.
  typedef typename Base::Index Index;

public:

  //! The number type.
  typedef typename Base::Number Number;
  //! A handle to a value type stored in an array.
  typedef typename Base::handle handle;
  //! A const handle to a value type stored in an array.
  typedef typename Base::const_handle const_handle;

private:

  //! The solution array.
  using Base::_solution;
  //! The finite difference scheme.
  using Base::_scheme;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  GridFM_BHDK();
  //! Copy constructor not implemented.
  GridFM_BHDK( const GridFM_BHDK& );
  //! Assignment operator not implemented.
  GridFM_BHDK& 
  operator=( const GridFM_BHDK& );

public:

  //
  // Constructors
  //

  //! Construct from the solution array and the grid spacing
  /*!
    \param solution is the solution grid.
    \param dx the spacing between adjacent grid points.
  */
  template <bool A>
  GridFM_BHDK( ads::Array<N,Number,A>& solution, const Number dx ) :
    Base( solution, dx )
  {}

  //! Trivial destructor.
  virtual 
  ~GridFM_BHDK() 
  {}

  //
  // Mathematical functions.
  //

  //! Solve the Hamilton-Jacobi equation with the FM method.
  /*!
    Find the solution for all grid points with solution less than or
    equal to \c max_solution.  If \c max_solution is zero, then the solution
    is determined for the entire grid.  The default value of \c max_solution
    is zero.

    This function uses the calls the \c label_neighbors() function 
    defined in the finite difference scheme.
  */
  void 
  solve( Number max_solution = 0 );
};

END_NAMESPACE_HJ

#define __hj_GridFM_BHDK_ipp__
#include "GridFM_BHDK.ipp"
#undef __hj_GridFM_BHDK_ipp__

#endif
