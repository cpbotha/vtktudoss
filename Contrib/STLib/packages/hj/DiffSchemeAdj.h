// -*- C++ -*-

/*! 
  \file DiffSchemeAdj.h
  \brief A class that supports finite difference operations for an N-D grid.
  
  Scheme with an adjacent stencil.
*/

#if !defined(__hj_DiffSchemeAdj_h__)
#define __hj_DiffSchemeAdj_h__

#include "DiffScheme.h"

#include "../ads/algorithm/min_max.h"

#include <limits>

#include <cmath>

// If we are debugging the whole hj namespace.
#if defined(DEBUG_hj) && !defined(DEBUG_DiffSchemeAdj)
#define DEBUG_DiffSchemeAdj
#endif

#ifdef DEBUG_DiffSchemeAdj
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Adjacent difference scheme.
/*!
  \param N is the space dimension.
  \param T is the number type.
  \param Equation represents the equation to be solved.  The equation must
  supply functions that perform the finite differencing in up to N
  adjacent directions.

  This class implements the labeling operations for adjacent difference 
  schemes in the \c label_neighbors() member function.
*/
template <int N, typename T, class Equation>
class DiffSchemeAdj;

END_NAMESPACE_HJ

#define __hj_DiffSchemeAdj2_h__
#include "DiffSchemeAdj2.h"
#undef __hj_DiffSchemeAdj2_h__

#define __hj_DiffSchemeAdj3_h__
#include "DiffSchemeAdj3.h"
#undef __hj_DiffSchemeAdj3_h__

#endif
