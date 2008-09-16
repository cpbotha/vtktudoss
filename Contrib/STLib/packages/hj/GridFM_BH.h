// -*- C++ -*-

/*! 
  \file GridFM_BH.h
  \brief N-D grid.  Fast Marching method.  Binary Heap with static keys.
*/

#if !defined(__hj_GridFM_BH_h__)
#define __hj_GridFM_BH_h__

#include "Grid.h"
#include "status.h"

#include "../ads/priority_queue/PriorityQueueBinaryHeapStoreKeys.h"

// If we are debugging the whole hj namespace.
#if defined(DEBUG_hj) && !defined(DEBUG_GridFM_BH)
#define DEBUG_GridFM_BH
#endif

BEGIN_NAMESPACE_HJ

//! The fast marching method for solving static H-J equations.
/*!
  This class implements the fast marching method for solving static 
  Hamilton-Jacobi equations in the \c solve() member function.  It uses
  a simple priority queue that does not have a decrease_key() operation.
*/
template <int N, typename T, class DifferenceScheme>
class GridFM_BH :
  public Grid<N,T,DifferenceScheme>
{
private:

  typedef Grid<N,T,DifferenceScheme> Base;

protected:

  //! The index type.
  typedef typename Base::Index Index;

public:

  //! The number type.
  typedef typename Base::Number Number;
  //! A handle to a value type stored in an array.
  typedef typename Base::handle handle;
  //! A const handle to a value type stored in an array.
  typedef typename Base::const_handle const_handle;

protected:

  //! The solution array.
  using Base::_solution;
  //! The finite difference scheme.
  using Base::_scheme;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  GridFM_BH();
  //! Copy constructor not implemented.
  GridFM_BH( const GridFM_BH& );
  //! Assignment operator not implemented.
  GridFM_BH& 
  operator=( const GridFM_BH& );

public:

  //
  // Constructors
  //

  //! Construct from the solution array and the grid spacing
  /*!
    \param solution The solution grid.
    \param dx The spacing between adjacent grid points.
  */
  template <bool A>
  GridFM_BH( ads::Array<N,Number,A>& solution, const Number dx ) :
    Base( solution, dx )
  {}

  //! Trivial destructor.
  virtual 
  ~GridFM_BH() 
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

#define __hj_GridFM_BH_ipp__
#include "GridFM_BH.ipp"
#undef __hj_GridFM_BH_ipp__

#endif
