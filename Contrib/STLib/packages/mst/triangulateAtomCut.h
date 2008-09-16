// -*- C++ -*-

/*! 
  \file triangulateAtomCut.h
  \brief Function for tesselating an atom using cut clipping.
*/

#if !defined(__mst_triangulateAtomCut_h__)
#define __mst_triangulateAtomCut_h__

#include "triangulateAtom.h"

#include "../ads/functor/constant.h"
#include "../geom/kernel/Circle3.h"
#include "../geom/mesh/simplicial/coarsen.h"
#include "../numerical/constants/Constants.h"

BEGIN_NAMESPACE_MST


//! Triangulate the visible surface.
/*!
  \param atom The atom to triangulate.
  \param clippingIdentifiers The identifiers of the atoms that might 
  clip the surface.
  \param clippingAtoms The atoms that might clip the surface.
  \param edgeLengthSlope The slope of the maximum edge length function.
  \param edgeLengthOffset The offset of the maximum edge length function.
  \param refinementLevel The refinement level for the initial mesh.  
  (Used only if it is non-negative.)
  \param mesh The output triangle mesh.
  \param actuallyClip The identifiers of the atoms that actually clip the mesh.
  \param epsilon If a vertex is within epsilon of the clipping surface it 
  is considered to be on the clipping surface.

  \return The target edge length.
*/
template<typename T, typename V, typename IS, typename IntOutputIterator>
T
triangulateVisibleSurfaceWithCutClipping
(const Atom<T>& atom,
 std::vector<int>& clippingIdentifiers,
 std::vector< Atom<T> >& clippingAtoms,
 T edgeLengthSlope,
 T edgeLengthOffset,
 int refinementLevel,
 geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh,
 IntOutputIterator actuallyClip,
 T epsilon);


END_NAMESPACE_MST

#define __mst_triangulateAtomCut_ipp__
#include "triangulateAtomCut.ipp"
#undef __mst_triangulateAtomCut_ipp__

#endif
