// -*- C++ -*-

/*! 
  \file EikonalAdj1st.h
  \brief Eikonal equation.  First-order, adjacent scheme.
*/

#if !defined(__hj_EikonalAdj1st_h__)
#define __hj_EikonalAdj1st_h__

#include "Eikonal.h"
#include "EikonalScheme.h"

// If we are debugging the whole hj package.
#if defined(DEBUG_hj) && !defined(DEBUG_EikonalAdj1st)
#define DEBUG_EikonalAdj1st
#endif

#ifdef DEBUG_EikonalAdj1st
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
class EikonalAdj1st;

END_NAMESPACE_HJ

#define __hj_EikonalAdj1st2_h__
#include "EikonalAdj1st2.h"
#undef __hj_EikonalAdj1st2_h__

#define __hj_EikonalAdj1st3_h__
#include "EikonalAdj1st3.h"
#undef __hj_EikonalAdj1st3_h__

#endif
