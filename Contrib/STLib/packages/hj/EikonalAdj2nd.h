// -*- C++ -*-

/*! 
  \file EikonalAdj2nd.h
  \brief Eikonal equation.  Second-order, adjacent scheme.
*/

#if !defined(__hj_EikonalAdj2nd_h__)
#define __hj_EikonalAdj2nd_h__

#include "Eikonal.h"
#include "EikonalScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_EikonalAdj2nd)
#define DEBUG_EikonalAdj2nd
#endif

#ifdef DEBUG_EikonalAdj2nd
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Eikonal equation.  Adjacent difference scheme.  2nd order.
/*!
  \param N is the space dimension.
  \param T is the number type.
*/
template<int N, typename T>
class EikonalAdj2nd;

END_NAMESPACE_HJ

#define __hj_EikonalAdj2nd2_h__
#include "EikonalAdj2nd2.h"
#undef __hj_EikonalAdj2nd2_h__

//#define __hj_EikonalAdj2nd3_h__
//#include "EikonalAdj2nd3.h"
//#undef __hj_EikonalAdj2nd3_h__

#endif
