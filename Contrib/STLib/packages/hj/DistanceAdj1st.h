// -*- C++ -*-

/*! 
  \file DistanceAdj1st.h
  \brief Distance equation.  First-order, adjacent scheme.
*/

#if !defined(__hj_DistanceAdj1st_h__)
#define __hj_DistanceAdj1st_h__

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_DistanceAdj1st)
#define DEBUG_DistanceAdj1st
#endif

#ifdef DEBUG_DistanceAdj1st
// Include the debugging code.
#include "debug.h"
#endif

#include "Distance.h"
#include "DistanceScheme.h"

BEGIN_NAMESPACE_HJ

//! Distance equation.  Adjacent difference scheme.  1st order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template<int N, typename T>
class DistanceAdj1st;

END_NAMESPACE_HJ

#define __hj_DistanceAdj1st2_h__
#include "DistanceAdj1st2.h"
#undef __hj_DistanceAdj1st2_h__

#define __hj_DistanceAdj1st3_h__
#include "DistanceAdj1st3.h"
#undef __hj_DistanceAdj1st3_h__

#endif
