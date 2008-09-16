// -*- C++ -*-

/*! 
  \file EikonalAdjDiag1st.h
  \brief Eikonal equation.  First-order, adjacent-diagonal scheme.
*/

#if !defined(__EikonalAdjDiag1st_h__)
#define __EikonalAdjDiag1st_h__

#include "Eikonal.h"
#include "EikonalScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_EikonalAdjDiag1st)
#define DEBUG_EikonalAdjDiag1st
#endif

#ifdef DEBUG_EikonalAdjDiag1st
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Eikonal equation.  Adjacent difference scheme.  1st order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template<int N, typename T>
class EikonalAdjDiag1st;

END_NAMESPACE_HJ

#define __hj_EikonalAdjDiag1st2_h__
#include "EikonalAdjDiag1st2.h"
#undef __hj_EikonalAdjDiag1st2_h__

#define __hj_EikonalAdjDiag1st3_h__
#include "EikonalAdjDiag1st3.h"
#undef __hj_EikonalAdjDiag1st3_h__

#endif
