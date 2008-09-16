// -*- C++ -*-

/*! 
  \file GridBFS.h
  \brief Breadth First Search method.  Does not produce the correct solution.
*/

#if !defined(__hj_GridBFS_h__)
#define __hj_GridBFS_h__

#include "Grid.h"

#include <deque>

// If we are debugging the whole hj namespace.
#if defined(DEBUG_hj) && !defined(DEBUG_GridBFS)
#define DEBUG_GridBFS
#endif

BEGIN_NAMESPACE_HJ

//! The breadth first search algorithm for "solving" static H-J equations.
/*!
  This class implements a placebo method for "solving" static 
  Hamilton-Jacobi equations in the \c solve() member function.
*/
template <int N, typename T, class DifferenceScheme>
class GridBFS :
  public Grid<N,T,DifferenceScheme>
{
private:

  typedef Grid<N,T,DifferenceScheme> Base;
  //! The index type.
  typedef typename Base::Index Index;

public:

  //! The number type.
  typedef T Number;

private:

  //! The solution array.
  using Base::_solution;
  //! The index ranges of the solution, not counting ghost points.
  //using Base::_ranges;
  //! The grid spacing.
  //using Base::_dx;
  //! The finite difference scheme.
  using Base::_scheme;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  GridBFS();
  //! Copy constructor not implemented.
  GridBFS( const GridBFS& );
  //! Assignment operator not implemented.
  GridBFS& 
  operator=( const GridBFS& );

public:

  //
  // Constructors
  //

  //! Construct from the solution array and the grid spacing.
  /*!
    \param solution is the solution grid.
    \param dx the spacing between adjacent grid points.
  */
  template <bool A>
  GridBFS( ads::Array<N,Number,A>& solution, const Number dx ) :
    Base( solution, dx )
  {}

  //! Trivial destructor.
  virtual 
  ~GridBFS() 
  {}

  //
  // Mathematical functions.
  //

  //! "Solve" the Hamilton-Jacobi equation with the placebo method.
  /*!
    Order the finite difference labeling operations with a breadth-first
    search of the grid points.
    
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

#define __hj_GridBFS_ipp__
#include "GridBFS.ipp"
#undef __hj_GridBFS_ipp__

#endif
