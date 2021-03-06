// -*- C++ -*-

#if !defined(__numerical_h__)
#define __numerical_h__

#include "numerical/constants.h"
#include "numerical/derivative.h"
#include "numerical/grid_interp_extrap.h"
#include "numerical/interpolation.h"
#include "numerical/optimization.h"
#include "numerical/partition.h"
#include "numerical/polynomial.h"
#include "numerical/random.h"
#include "numerical/specialFunctions.h"

BEGIN_NAMESPACE_NUMERICAL

/*!
\mainpage Numerical Algorithms Package
\anchor numerical

\section numerical_introduction Introduction

This is a numerical algorithms package that I use in various
projects.  It's not a general purpose library.  I just add
functionality as I need it.  

This is a templated C++ class library.  All the functionality is
implemented in header files.  Thus there is no library to compile or
link with.  Just include the appropriate header files in your
application code when you compile.

This package is composed of a number of sub-packages.  All classes
and functions are in the \c numerical namespace.
- The \ref numerical_constants "Mathematical constants" micro-package has 
  constants for \f$\pi\f$, Euler's constant e, and the like.
- The \ref derivative "derivative" package has functions and functors for 
  evaluating derivatives.
- The \ref grid_interp_extrap "grid interpolation/extrapolation" package
  is useful for interpolating field values in level-set applications.
- The \ref interpolation has functions for performing linear 
  interpolation.
- The \ref optimization "optimization" package implements a quasi-Newton
  method, a downhill simplex method and a coordinate descent method.  In 
  addition, it implements the penalty method for constrained optimization
  with equality constraints.
- The \ref partition "partition" micro-package has a function for 
  fair partitioning of an integer.
- The \ref numerical_polynomial "polynomial" package has a function for 
  evaluating polynomials.
- The \ref numerical_random "random number" package has functors for
  generating random numbers.
- The \ref numerical_specialFunctions "special functions" package has 
  functors for \f$\Gamma\f$ and the like.
*/

END_NAMESPACE_NUMERICAL

#endif
