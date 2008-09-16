// -*- C++ -*-

/*! 
  \file tesselate_sphere_unindexed.h
  \brief Function for tesselating a unit sphere.
*/

#if !defined(__tesselate_sphere_unindexed_h__)
#define __tesselate_sphere_unindexed_h__

#include "defs.h"

#include "../ads/array/Array.h"
#include "../geom/kernel/Point.h"

#include <list>

#include <cmath>
#include <cassert>

BEGIN_NAMESPACE_MST


END_NAMESPACE_MST

#define __tesselate_sphere_unindexed_ipp__
#include "tesselate_sphere_unindexed.ipp"
#undef __tesselate_sphere_unindexed_ipp__

#endif
