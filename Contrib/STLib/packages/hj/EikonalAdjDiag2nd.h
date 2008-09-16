// -*- C++ -*-

/*! 
  \file EikonalAdjDiag2nd.h
  \brief Eikonal equation.  Second-order, adjacent-diagonal scheme.
*/

#if !defined(__EikonalAdjDiag2nd_h__)
#define __EikonalAdjDiag2nd_h__

#include "Eikonal.h"
#include "EikonalScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_EikonalAdjDiag2nd)
#define DEBUG_EikonalAdjDiag2nd
#endif

#ifdef DEBUG_EikonalAdjDiag2nd
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Eikonal equation.  Adjacent-diagonal difference scheme.  2nd order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template<int N, typename T>
class EikonalAdjDiag2nd;

END_NAMESPACE_HJ

#define __hj_EikonalAdjDiag2nd2_h__
#include "EikonalAdjDiag2nd2.h"
#undef __hj_EikonalAdjDiag2nd2_h__

//#define __hj_EikonalAdjDiag2nd3_h__
//#include "EikonalAdjDiag2nd3.h"
//#undef __hj_EikonalAdjDiag2nd3_h__

#endif
