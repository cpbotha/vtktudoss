// -*- C++ -*-

/*! 
  \file DistanceAdjDiag2nd2.h
  \brief Distance equation.  Second-order, adjacent-diagonal scheme.
*/

#if !defined(__hj_DistanceAdjDiag2nd2_h__)
#error This file is an implementation detail of the class DistanceAdjDiag2nd.
#endif

BEGIN_NAMESPACE_HJ

//! Distance equation.  Adjacent-diagonal difference scheme.  2nd order.
/*!
  \param T is the number type.
*/
template <typename T>
class DistanceAdjDiag2nd<2,T> :
  public Distance<2,T>,
  public DistanceScheme<2,T>
{
private:

  typedef Distance<2,T> EquationBase;
  typedef DistanceScheme<2,T> SchemeBase;

  using EquationBase::diff_a1;
  using EquationBase::diff_a2;
  using EquationBase::diff_d1;
  using EquationBase::diff_d2;
  using EquationBase::diff_a1_a1;
  using EquationBase::diff_a1_d1;
  using EquationBase::diff_a1_d2;
  using EquationBase::diff_a2_d1;
  using EquationBase::diff_a2_d2;

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
  DistanceAdjDiag2nd();
  //! Copy constructor not implemented.
  DistanceAdjDiag2nd( const DistanceAdjDiag2nd& );
  //! Assignment operator not implemented.
  DistanceAdjDiag2nd& 
  operator=( const DistanceAdjDiag2nd& );

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
  DistanceAdjDiag2nd( ads::Array<2,T,A1>& solution, 
		      ads::Array<2,Status,A2>& status,
		      const Number dx ) :
    EquationBase( dx ),
    SchemeBase( solution, status )
  {}

  ~DistanceAdjDiag2nd()
  {}

  //
  // Accessors
  //

  //! Return the radius of the stencil.
  static
  int
  radius()
  {
    return 2;
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

  //! Use the two specified directions to compute the solution at i.
  Number
  diff_adj_diag( const Index& i, const Index& a, const Index& b ) const;

  //! Use the direction d to compute the solution at i.
  Number
  diff_diag( const Index& i, const Index& d ) const;

  //! Use the two specified directions to compute the solution at i.
  Number
  diff_diag_adj( const Index& i, const Index& a, const Index& b ) const;
};

END_NAMESPACE_HJ

#define __hj_DistanceAdjDiag2nd2_ipp__
#include "DistanceAdjDiag2nd2.ipp"
#undef __hj_DistanceAdjDiag2nd2_ipp__
