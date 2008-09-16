// -*- C++ -*-

/*! 
  \file DistanceAdjDiag2nd.h
  \brief Distance equation.  Second-order, adjacent-diagonal scheme.
*/

#if !defined(__DistanceAdjDiag2nd_h__)
#define __DistanceAdjDiag2nd_h__

#include "Distance.h"
#include "DistanceScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_DistanceAdjDiag2nd)
#define DEBUG_DistanceAdjDiag2nd
#endif

#ifdef DEBUG_DistanceAdjDiag2nd
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Distance equation.  Adjacent-diagonal difference scheme.  2nd order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template <int N, typename T>
class DistanceAdjDiag2nd;

END_NAMESPACE_HJ

#define __hj_DistanceAdjDiag2nd2_h__
#include "DistanceAdjDiag2nd2.h"
#undef __hj_DistanceAdjDiag2nd2_h__

//#define __hj_DistanceAdjDiag2nd3_h__
//#include "DistanceAdjDiag2nd3.h"
//#undef __hj_DistanceAdjDiag2nd3_h__

#endif
