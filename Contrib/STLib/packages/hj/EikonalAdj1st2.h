// -*- C++ -*-

/*! 
  \file EikonalAdj1st2.h
  \brief Eikonal equation.  First-order, adjacent scheme.
*/

#if !defined(__hj_EikonalAdj1st2_h__)
#error This file is an implementation detail of the class EikonalAdj1st.
#endif

BEGIN_NAMESPACE_HJ

//! Eikonal equation.  Adjacent difference scheme.  1st order.
/*!
  \param T is the number type.
*/
template<typename T>
class EikonalAdj1st<2,T> :
  public Eikonal<2,T>,
  public EikonalScheme<2,T> {
private:

  typedef Eikonal<2,T> EquationBase;
  typedef EikonalScheme<2,T> SchemeBase;

  using EquationBase::diff_a1;
  using EquationBase::diff_a1_a1;

  using SchemeBase::_status;
  using SchemeBase::_solution;
  using SchemeBase::_inverseSpeed;

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
  EikonalAdj1st();
  //! Copy constructor not implemented.
  EikonalAdj1st(const EikonalAdj1st&);
  //! Assignment operator not implemented.
  EikonalAdj1st& 
  operator=(const EikonalAdj1st&);

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
  EikonalAdj1st(ads::Array<2,T,A1>& solution, 
		ads::Array<2,Status,A2>& status,
		const Number dx) :
    EquationBase(dx),
    SchemeBase(solution, status)
  {}

  ~EikonalAdj1st()
  {}

  //
  // Accessors.
  //

  //! Return the radius of the stencil.
  static
  int
  radius() {
    return 1;
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
    return min_unknown; 
  }

  //! Use the grid point i+d to compute the solution at i.
  Number
  diff_adj(const Index& i, const Index& d) const;

  //! Use the two specified neighbors to compute the solution at i.
  Number
  diff_adj_adj(const Index& i,	const Index& a, const Index& b) const;

};

END_NAMESPACE_HJ

#define __hj_EikonalAdj1st2_ipp__
#include "EikonalAdj1st2.ipp"
#undef __hj_EikonalAdj1st2_ipp__
