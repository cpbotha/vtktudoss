// -*- C++ -*-

/*! 
  \file DiffSchemeAdjDiag.h
  \brief A class that supports finite difference operations for an N-D grid.
  
  Scheme with an adjacent-diagonal stencil.
*/

#if !defined(__hj_DiffSchemeAdjDiag_h__)
#define __hj_DiffSchemeAdjDiag_h__

#include "DiffScheme.h"

#include <limits>

#include <cmath>

// If we are debugging the whole hj namespace.
#if defined(DEBUG_hj) && !defined(DEBUG_DiffSchemeAdjDiag)
#define DEBUG_DiffSchemeAdjDiag
#endif

#ifdef DEBUG_DiffSchemeAdjDiag
// Include the debugging code.
#include "debug.h"
#endif

BEGIN_NAMESPACE_HJ

//! Adjacent-diagonal difference scheme.
/*!
  \param N is the space dimension.
  \param T is the number type.
  \param Equation represents the equation to be solved.  The equation must
  supply functions that perform the finite differencing in adjacent 
  and diagonal directions.

  This class implements the labeling operations for adjacent-diagonal 
  difference schemes in the \c label_neighbors() member function.
*/
template <int N, typename T, class Equation>
class DiffSchemeAdjDiag;

END_NAMESPACE_HJ

#define __hj_DiffSchemeAdjDiag2_h__
#include "DiffSchemeAdjDiag2.h"
#undef __hj_DiffSchemeAdjDiag2_h__

#define __hj_DiffSchemeAdjDiag3_h__
#include "DiffSchemeAdjDiag3.h"
#undef __hj_DiffSchemeAdjDiag3_h__

#endif
