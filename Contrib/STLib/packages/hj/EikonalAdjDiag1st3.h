// -*- C++ -*-

/*! 
  \file EikonalAdjDiag1st3.h
  \brief Eikonal equation.  First-order, adjacent-diagonal scheme.
*/

#if !defined(__hj_EikonalAdjDiag1st3_h__)
#error This file is an implementation detail of the class EikonalAdjDiag1st.
#endif

BEGIN_NAMESPACE_HJ

//! Eikonal equation.  Adjacent-diagonal difference scheme.  1st order.
/*!
  \param T is the number type.
*/
template<typename T>
class EikonalAdjDiag1st<3,T> :
  public Eikonal<3,T>,
  public EikonalScheme<3,T> {
private:

  typedef Eikonal<3,T> EquationBase;
  typedef EikonalScheme<3,T> SchemeBase;

  using EquationBase::diff_a1;
  using EquationBase::diff_d1;
  using EquationBase::diff_a1_d1;
  using EquationBase::diff_d1_d1;
  using EquationBase::diff_a1_d1_d1;
  using EquationBase::diff_d1_d1_d1;

  using EquationBase::_dx_t_sqrt2;
  using EquationBase::_dx_o_sqrt2;

  using SchemeBase::_status;
  using SchemeBase::_solution;
  using SchemeBase::_inverseSpeed;

public:

  //! The number type.
  typedef T Number;
  //! A multi-index.
  typedef ads::FixedArray<3,int> Index;

private:

  // 
  // Not implemented.
  //

  //! Default constructor not implemented.
  EikonalAdjDiag1st();
  //! Copy constructor not implemented.
  EikonalAdjDiag1st(const EikonalAdjDiag1st&);
  //! Assignment operator not implemented.
  EikonalAdjDiag1st& 
  operator=(const EikonalAdjDiag1st&);

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
  template<bool A1, bool A2>
  EikonalAdjDiag1st(ads::Array<3,T,A1>& solution, 
		    ads::Array<3,Status,A2>& status,
		    const Number dx) :
    EquationBase(dx),
    SchemeBase(solution, status)
  {}

  ~EikonalAdjDiag1st()
  {}

  //
  // Accessors
  //

  //! Return the radius of the stencil.
  static
  int
  radius() {
    return 1;
  }

  //! Return the minimum change from a known solution to a labeled solution when labeling a neighbor.
  /*!
    This is an expensive function.  It loops over the speed array.
  */
  Number
  min_delta() const {
    return _dx_o_sqrt2 * *std::max_element(_inverseSpeed.begin(), 
					   _inverseSpeed.end());
  }

  //! Return the maximum change from a known solution to a labeled solution when labeling a neighbor.
  /*!
    This is an expensive function.  It loops over the speed array.
  */
  Number
  max_delta() const {
    return _dx_t_sqrt2 * *std::min_element(_inverseSpeed.begin(), 
					   _inverseSpeed.end());
  }

  //
  // Manipulators.
  //

  using SchemeBase::getInverseSpeed;

  //
  // Mathematical member functions.
  //

  //! Return a lower bound on correct values.
  Number 
  lower_bound(const Index& i, const Number min_unknown) const { 
    return min_unknown + _dx_o_sqrt2 * _inverseSpeed(i);
  }

  //! Use the grid point i+d to compute the solution at i.
  Number
  diff_adj(const Index& i, const Index& d) const;

  //! Use the two specified neighbors to compute the solution at i.
  Number
  diff_adj_diag(const Index& i, const Index& a, const Index& b) const;

  //! Use the three specified neighbors to compute the solution at i.
  Number
  diff_adj_diag_diag(const Index& i, const Index& a, const Index& b, 
		     const Index& c) const;

  //! Use the grid point i+d to compute the solution at i.
  Number
  diff_diag(const Index& i, const Index& d) const;

  //! Use the two specified neighbors to compute the solution at i.
  Number
  diff_diag_adj(const Index& i, const Index& a, const Index& b) const;

  //! Use the two specified neighbors to compute the solution at i.
  Number
  diff_diag_diag(const Index& i, const Index& a, const Index& b) const;

  //! Use the three specified neighbors to compute the solution at i.
  Number
  diff_diag_adj_diag(const Index& i, const Index& a, const Index& b, 
		     const Index& c) const;

  //! Use the three specified neighbors to compute the solution at i.
  Number
  diff_diag_diag_diag(const Index& i, const Index& a, const Index& b, 
		      const Index& c) const;
};

END_NAMESPACE_HJ

#define __hj_EikonalAdjDiag1st3_ipp__
#include "EikonalAdjDiag1st3.ipp"
#undef __hj_EikonalAdjDiag1st3_ipp__
