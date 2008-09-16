// -*- C++ -*-

/*! 
  \file DistanceAdjDiag1st2.h
  \brief Distance equation.  First-order, adjacent-diagonal scheme.
*/

#if !defined(__hj_DistanceAdjDiag1st2_h__)
#error This file is an implementation detail of the class DistanceAdjDiag1st.
#endif

BEGIN_NAMESPACE_HJ

//! Distance equation.  Adjacent-diagonal difference scheme.  1st order.
/*!
  \param T is the number type.
*/
template <typename T>
class DistanceAdjDiag1st<2,T> :
  public Distance<2,T>,
  public DistanceScheme<2,T>
{
private:

  typedef Distance<2,T> EquationBase;
  typedef DistanceScheme<2,T> SchemeBase;

  using EquationBase::diff_a1;
  using EquationBase::diff_d1;
  using EquationBase::diff_a1_d1;

  using EquationBase::_dx_t_sqrt2;
  using EquationBase::_dx_o_sqrt2;

  using SchemeBase::_status;
  using SchemeBase::_solution;

public:

  //! The number type.
  typedef T Number;
  //! A multi-index.
  typedef ads::FixedArray<2,int> Index;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  DistanceAdjDiag1st();
  //! Copy constructor not implemented.
  DistanceAdjDiag1st( const DistanceAdjDiag1st& );
  //! Assignment operator not implemented.
  DistanceAdjDiag1st& 
  operator=( const DistanceAdjDiag1st& );

public:

  //
  // Constructors
  //

  //! Constructor.
  /*!
    \param solution is the solution array.
    \param status is the status array.
    \param dx is the grid spacing.
  */
  template <bool A1, bool A2>
  DistanceAdjDiag1st( ads::Array<2,T,A1>& solution, 
		      ads::Array<2,Status,A2>& status,
		      const Number dx ) :
    EquationBase( dx ),
    SchemeBase( solution, status )
  {}

  ~DistanceAdjDiag1st()
  {}

  //
  // Accessors
  //

  //! Return the radius of the stencil.
  static
  int
  radius()
  {
    return 1;
  }

  //! Return the minimum change from a known solution to a labeled solution when labeling a neighbor.
  Number
  min_delta() const
  {
    return _dx_o_sqrt2;
  }

  //! Return the maximum change from a known solution to a labeled solution when labeling a neighbor.
  Number
  max_delta() const
  {
    return _dx_t_sqrt2;
  }

  //
  // Mathematical member functions.
  //

  //! Return a lower bound on correct values.
  Number 
  lower_bound( const Index& i, const Number min_unknown ) const
  { 
    return min_unknown + _dx_o_sqrt2;
  }

  //! Use the grid point i+d to compute the solution at i.
  Number
  diff_adj( const Index& i, const Index& d ) const;

  //! Use the two specified neighbors to compute the solution at i.
  Number
  diff_adj_diag( const Index& i, const Index& a, const Index& b ) const;

  //! Use the grid point i+d to compute the solution at i.
  Number
  diff_diag( const Index& i, const Index& d ) const;

  //! Use the two specified neighbors to compute the solution at i.
  Number
  diff_diag_adj( const Index& i, const Index& a, const Index& b ) const;
};

END_NAMESPACE_HJ

#define __hj_DistanceAdjDiag1st2_ipp__
#include "DistanceAdjDiag1st2.ipp"
#undef __hj_DistanceAdjDiag1st2_ipp__
