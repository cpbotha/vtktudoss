// -*- C++ -*-

#if !defined(__mst_triangulateAtom_ipp__)
#error This file is an implementation detail of triangulateAtom.
#endif

BEGIN_NAMESPACE_MST


template<typename T, typename V, typename IS, typename IntOutputIterator>
inline
T
triangulateVisibleSurface(const Atom<T>& atom,
			  std::vector<int>& clippingIdentifiers,
			  std::vector< Atom<T> >& clippingAtoms,
			  const T edgeLengthSlope,
			  const T edgeLengthOffset,
			  const int refinementLevel,
			  geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh,
			  IntOutputIterator actuallyClip,
			  const T maximumStretchFactor,
			  const T epsilon,
			  const bool areUsingCircularEdges) {
  // Compute the clipping atoms and form the initial tesselation.
  const T targetEdgeLength = 
    computeClippingAtomsAndTesselate(atom, clippingIdentifiers,
				     clippingAtoms, edgeLengthSlope,
				     edgeLengthOffset, refinementLevel,
				     mesh, actuallyClip);

  // If the atom is completely erased by the clipping.
  if (mesh->getSimplicesSize() == 0) {
    return targetEdgeLength;
  }

  // If we allow any stretching.
  if (maximumStretchFactor != 1) {
    clipWithRubberClipping(atom, clippingIdentifiers, clippingAtoms, mesh,
			   actuallyClip, maximumStretchFactor,
			   areUsingCircularEdges);
  }
  clipWithCutClipping(atom, clippingIdentifiers, clippingAtoms,
		      mesh, actuallyClip, epsilon);

  return targetEdgeLength;
}


END_NAMESPACE_MST

// End of file.
