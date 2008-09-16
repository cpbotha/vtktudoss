// -*- C++ -*-

/*! 
  \file DistanceAdj2nd.h
  \brief Distance equation.  Second-order, adjacent scheme.
*/

#if !defined(__hj_DistanceAdj2nd_h__)
#define __hj_DistanceAdj2nd_h__

#include "Distance.h"
#include "DistanceScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_DistanceAdj2nd)
#define DEBUG_DistanceAdj2nd
#endif

#ifdef DEBUG_DistanceAdj2nd
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Distance equation.  Adjacent difference scheme.  2nd order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template <int N, typename T>
class DistanceAdj2nd;

END_NAMESPACE_HJ

#define __hj_DistanceAdj2nd2_h__
#include "DistanceAdj2nd2.h"
#undef __hj_DistanceAdj2nd2_h__

//#define __hj_DistanceAdj2nd3_h__
//#include "DistanceAdj2nd3.h"
//#undef __hj_DistanceAdj2nd3_h__

#endif
