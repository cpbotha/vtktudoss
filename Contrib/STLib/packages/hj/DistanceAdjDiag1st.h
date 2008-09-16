// -*- C++ -*-

/*! 
  \file DistanceAdjDiag1st.h
  \brief Distance equation.  First-order, adjacent-diagonal scheme.
*/

#if !defined(__DistanceAdjDiag1st_h__)
#define __DistanceAdjDiag1st_h__

#include "Distance.h"
#include "DistanceScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_DistanceAdjDiag1st)
#define DEBUG_DistanceAdjDiag1st
#endif

#ifdef DEBUG_DistanceAdjDiag1st
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Distance equation.  Adjacent difference scheme.  1st order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template <int N, typename T>
class DistanceAdjDiag1st;

END_NAMESPACE_HJ

#define __hj_DistanceAdjDiag1st2_h__
#include "DistanceAdjDiag1st2.h"
#undef __hj_DistanceAdjDiag1st2_h__

#define __hj_DistanceAdjDiag1st3_h__
#include "DistanceAdjDiag1st3.h"
#undef __hj_DistanceAdjDiag1st3_h__

#endif
