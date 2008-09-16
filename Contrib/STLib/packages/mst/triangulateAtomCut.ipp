// -*- C++ -*-

#if !defined(__mst_triangulateAtomCut_ipp__)
#error This file is an implementation detail of triangulateAtomCut.
#endif

BEGIN_NAMESPACE_MST


// CONTINUE: Deal with special cases, like only zero/one exterior vertex or 
// only zero/one interior vertex.
template<typename T, typename V, typename IS>
inline
void
clipWithCuts(const Atom<T>& atom, const Atom<T>& clippingAtom, 
	     geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh,
	     const ads::Array<1,int>& signOfDistance) {
  // A Cartesian point.
  typedef typename Atom<T>::Point Point;
  // The mesh type.
  typedef geom::IndSimpSetIncAdj<3,2,true,T,V,IS> Mesh;
  // An iterator over 1-faces (edges) in the mesh.
  typedef typename Mesh::FaceIterator FaceIterator;
  // An indexed simplex.
  typedef typename Mesh::IndexedSimplex IndexedSimplex;
  // An indexed simplex const iterator.
  typedef typename Mesh::IndexedSimplexConstIterator 
    IndexedSimplexConstIterator;
  // A pair of vertex indices define an edge.
  typedef std::pair<int,int> Edge;

  // The simplex dimension.
  const int M = 2;

  // Make the circle that is the intersection of the two spheres.
  // We will use this in finding the intersection points below.
  geom::Circle3<T> circle;
  makeIntersectionCircle(atom, clippingAtom, &circle);

  //
  // Determine the edges that cross the clipping surface.  Compute the 
  // intersection points.
  //
  // This will keep track of the new vertex index as they are computed.
  int vertexIndex = mesh->getVerticesSize();
  // Each simplex will determine if each of its faces (edges) will have 
  // a new intersection vertex.  These are the indices of these vertices.
  // We initialize the array with the null value -1.
  ads::Array<2,int> 
    newVertexIndices(ads::FixedArray<2,int>(mesh->getSimplicesSize(), M + 1),
		     -1);
  std::vector<Point> newVertices;
  // For each edge.
  int simplexIndex, localIndex, vertexIndex1, vertexIndex2, adjacent;
  Point p;
  for (FaceIterator i = mesh->getFacesBeginning(); i != mesh->getFacesEnd();
       ++i) {
    simplexIndex = i->first;
    localIndex = i->second;
    vertexIndex1 = 
      mesh->getIndexedSimplex(simplexIndex)[(localIndex + 1) % (M + 1)];
    vertexIndex2 = 
      mesh->getIndexedSimplex(simplexIndex)[(localIndex + 2) % (M + 1)];
    // If the edge crosses the clipping curve.
    if (signOfDistance[vertexIndex1] * signOfDistance[vertexIndex2] == -1) {
      // Compute the intersection point.
      geom::computeClosestPoint(circle, mesh->getVertex(vertexIndex1),
				mesh->getVertex(vertexIndex2), &p);
      newVertices.push_back(p);
      // Record the new vertex index for this simplex.
      newVertexIndices(simplexIndex, localIndex) = vertexIndex;
      // Record the new vertex index for the adjacent simplex.
      adjacent = mesh->getAdjacent(simplexIndex, localIndex);
      // If there is an adjacent simplex.
      if (adjacent != -1) {
	newVertexIndices(adjacent, mesh->getMirrorIndex(simplexIndex, 
							localIndex)) 
	  = vertexIndex;
      }
      ++vertexIndex;
    }
  }

  //
  // Build a vector of the vertices for the clipped mesh.
  //
  std::vector<Point> vertices;
  // Reserve room for the old and new vertices.
  vertices.reserve(vertexIndex);
  // Add the old vertices.
  std::copy(mesh->getVerticesBeginning(), mesh->getVerticesEnd(), 
	    std::back_inserter(vertices));
  // Add the new vertices.
  std::copy(newVertices.begin(), newVertices.end(), 
	    std::back_inserter(vertices));

  //
  // Build a vector of the indexed simplices for the clipped mesh.
  //
  std::vector<IndexedSimplex> indexedSimplices;
  IndexedSimplex s;
  ads::FixedArray<M + 1,int> indices;
  typename Mesh::Vertex centroid;
  // For each triangle.
  for (simplexIndex = 0; simplexIndex != mesh->getSimplicesSize(); 
       ++simplexIndex) {
    // The indexed simplex.
    s = mesh->getIndexedSimplex(simplexIndex);

    // Make an array of the sign of the distance for this simplex.
    ads::FixedArray<M + 1,int> simplexSignOfDistance;
    for (int i = 0; i != M + 1; ++i) {
      simplexSignOfDistance[i] = signOfDistance[s[i]];
    }

    // -1 -1 -1 = -3  Keep.
    // -1 -1  0 = -2  Keep.
    // -1 -1  1 = -1  Make 2 triangles.
    // -1  0  0 = -1  Keep.
    // -1  0  1 =  0  Make 1 triangle, use 1 new vertex.
    // -1  1  1 =  1  Make 1 triangle, use 2 new vertices.
    //  0  0  0 =  0  Use the centroid to determine whether to keep or discard.
    //  0  0  1 =  1  Discard.
    //  0  1  1 =  2  Discard.
    //  1  1  1 =  3  Discard.

    // The sign of the distance in sorted order.
    ads::FixedArray<M + 1,int> sd(simplexSignOfDistance);
    sd.sort();

    // If all three vertices are visible or 
    // two are visible and one is on the clipping surface or
    // one is visible and two are on the clipping surface.
    if (sd[0] == -1 && sd[1] <= 0 && sd[2] <= 0) {
      // Add the old simplex.
      indexedSimplices.push_back(s);
    }
    // If two vertices are visible and one is not visible.
    else if (sd[0] == -1 && sd[1] == -1 && sd[2] == 1) {
      //
      // Move the non-visible vertex to the third position.
      //
      indices[0] = newVertexIndices(simplexIndex, 0);
      indices[1] = newVertexIndices(simplexIndex, 1);
      indices[2] = newVertexIndices(simplexIndex, 2);
      if (simplexSignOfDistance[0] == 1) {
	std::swap(s[0], s[1]);
	std::swap(s[1], s[2]);
	std::swap(indices[0], indices[1]);
	std::swap(indices[1], indices[2]);
      }
      else if (simplexSignOfDistance[1] == 1) {
	std::swap(s[1], s[2]);
	// Do the second swap to preserve the orientation.
	std::swap(s[0], s[1]);
	std::swap(indices[1], indices[2]);
	std::swap(indices[0], indices[1]);
      }

      // The two triangles that compose the visible quadrilateral.
      s[2] = indices[0];
      indexedSimplices.push_back(s);
      s[1] = indices[0];
      s[2] = indices[1];
      indexedSimplices.push_back(s);
    }
    // If one vertex is visible, one is on the clipping curve,
    // and one is not visible.
    else if (sd[0] == -1 && sd[1] == 0 && sd[2] == 1) {
      //
      // Move the visible vertex to the first position.
      //
      indices[0] = newVertexIndices(simplexIndex, 0);
      indices[1] = newVertexIndices(simplexIndex, 1);
      indices[2] = newVertexIndices(simplexIndex, 2);
      if (simplexSignOfDistance[1] == -1) {
	std::swap(s[0], s[1]);
	// Do the second swap to preserve the orientation.
	std::swap(s[1], s[2]);
	std::swap(indices[0], indices[1]);
	std::swap(indices[1], indices[2]);
	std::swap(simplexSignOfDistance[0], simplexSignOfDistance[1]);
	std::swap(simplexSignOfDistance[1], simplexSignOfDistance[2]);
      }
      else if (simplexSignOfDistance[2] == -1) {
	std::swap(s[1], s[2]);
	std::swap(s[0], s[1]);
	std::swap(indices[1], indices[2]);
	std::swap(indices[0], indices[1]);
	std::swap(simplexSignOfDistance[1], simplexSignOfDistance[2]);
	std::swap(simplexSignOfDistance[0], simplexSignOfDistance[1]);
      }
      // If the second vertex is on the clipping surface.
      if (simplexSignOfDistance[1] == 0) {
	s[2] = indices[1];
      }
      // Otherwise, the third vertex is on the clipping surface.
      else {
	s[1] = indices[2];
      }
      indexedSimplices.push_back(s);
    }
    // If one vertex is visible and two are not visible.
    else if (sd[0] == -1 && sd[1] == 1 && sd[2] == 1) {
      //
      // Move the visible vertex to the first position.
      //
      indices[0] = newVertexIndices(simplexIndex, 0);
      indices[1] = newVertexIndices(simplexIndex, 1);
      indices[2] = newVertexIndices(simplexIndex, 2);
      if (simplexSignOfDistance[1] == -1) {
	std::swap(s[0], s[1]);
	// Do the second swap to preserve the orientation.
	std::swap(s[1], s[2]);
	std::swap(indices[0], indices[1]);
	std::swap(indices[1], indices[2]);
      }
      else if (simplexSignOfDistance[2] == -1) {
	std::swap(s[1], s[2]);
	std::swap(s[0], s[1]);
	std::swap(indices[1], indices[2]);
	std::swap(indices[0], indices[1]);
      }

      // The visible portion of the triangle is also a triangle.
      s[1] = indices[2];
      s[2] = indices[1];
      indexedSimplices.push_back(s);
    }
    // If three vertices are on the clipping surface.
    else if (sd[0] == 0 && sd[1] == 0 && sd[2] == 0) {
      // Use the centroid to determine whether to keep or discard the triangle.
      geom::getCentroid(*mesh, simplexIndex, &centroid);
      // If the centroid has negative distance.
      if (clippingAtom.getRadius() - 
	  geom::computeDistance(clippingAtom.getCenter(), centroid) < 0) {
	// Keep the triangle.
	indexedSimplices.push_back(s);
      }
    }
    // Otherwise no vertices are visible.
    else {
      assert(sd[0] >= 0 && sd[1] >= 0 && sd[2] == 1);
      // The triangle is removed by the clipping.
    }
  }

  // Make the new mesh.
  mesh->build(vertices.size(), &vertices[0], indexedSimplices.size(),
	      &indexedSimplices[0]);
  // Pack the mesh to get rid of unused vertices.
  geom::pack(mesh);
}



// Clip the mesh.  Return true if any clipping was done.
template<typename T, typename V, typename IS>
inline
bool
clipWithCuts(const Atom<T>& atom, const Atom<T>& clippingAtom, 
	     geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh,
	     const T epsilon) {
  typedef geom::IndSimpSetIncAdj<3,2,true,T,V,IS> Mesh;

  assert(epsilon >= 0);

  // Compute the negative of the signed distance from the vertices to the 
  // clipping atom.  Then the vertices with negative distance are visible.
  ads::Array<1,T> distance(mesh->getVerticesSize());
  for (int i = 0; i != distance.size(); ++i) {
    distance[i] = clippingAtom.getRadius() 
      - geom::computeDistance(clippingAtom.getCenter(), mesh->getVertex(i));
  }

  // Make an array of the sign of the distance.
  ads::Array<1,int> signOfDistance(distance.size());
  for (int i = 0; i != distance.size(); ++i) {
    // If the distance is within epsilon of zero, we consider it to be zero.
    if (distance[i] > epsilon) {
      signOfDistance[i] = 1;
    }
    else if (distance[i] < -epsilon) {
      signOfDistance[i] = -1;
    }
    else {
      signOfDistance[i] = 0;
    }
  }

  // The number of vertices outside the visible domain.
  const int numOutside = std::count(signOfDistance.begin(), 
				    signOfDistance.end(), 1);
  
  // If all of the vertices are visible.
  if (numOutside == 0) {
    // Do not alter the mesh.
    return false;
  }
  // If none of the vertices are visible.
  if (numOutside == 0) {
    // Clear the mesh.
    *mesh = Mesh();
  }
  else {
    // Otherwise, clip the mesh.
    clipWithCuts(atom, clippingAtom, mesh, signOfDistance);
  }
  return true;
}



template<typename T, typename V, typename IS>
inline
void
optimize(const T maxEdgeLength,
	 geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh) {
  typedef geom::SimpMeshRed<3,2,T> SimplicialMesh;
  typedef typename SimplicialMesh::Vertex Vertex;

#if 0
  // CONTINUE
  std::cerr << "begin optimize " << mesh->getVerticesSize() << " "
	    << mesh->getSimplicesSize() << " "
	    << maxEdgeLength << "\n";
#endif

  // Make a mesh for mesh optimization.
  SimplicialMesh optimizationMesh(*mesh);

  // Apply coarsening to get rid of very short edges.
  {
    const T minimumAllowedQuality = 0.0;
    const T qualityFactor = 0.0;
    geom::coarsen<geom::SimplexModCondNum<2> >
      (&optimizationMesh, 
       ads::unary_constant<Vertex,T>(0.01 * maxEdgeLength), 
       minimumAllowedQuality, qualityFactor);
  }

#if 0
  // CONTINUE: This crashes.
  // Apply coarsening to get rid of short edges and improve the quality of 
  // the mesh.
  const T minimumAllowedQuality = 0.0;
  const T qualityFactor = 1.01;
  const T edgeDeviation = numerical::Constants<T>::Pi() / 6.0;
  const T cornerDeviation = numerical::Constants<T>::Pi() / 6.0;
  geom::coarsen<geom::SimplexModCondNum<2> >
    (&optimizationMesh, 
     ads::unary_constant<Vertex,T>(0.1 * maxEdgeLength), 
     minimumAllowedQuality,qualityFactor,
     edgeDeviation,cornerDeviation);
#endif

  // Transfer the result back to the ISS.
  // CONTINUE
  //std::cerr << "transfer\n";
  buildIndSimpSetFromSimpMeshRed(optimizationMesh, mesh);
  // CONTINUE
  //std::cerr << "done " << mesh->getVerticesSize() << " "
  //<< mesh->getSimplicesSize() << "\n";

#if 0
  // CONTINUE
  std::cerr << "end optimize " << mesh->getVerticesSize() << " "
	    << mesh->getSimplicesSize() << " "
	    << maxEdgeLength << "\n";
#endif
}



template<typename T, typename V, typename IS, typename IntOutputIterator>
inline
void
clipWithCutClipping(const Atom<T>& atom,
		    std::vector<int>& clippingIdentifiers,
		    std::vector< Atom<T> >& clippingAtoms,
		    geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh,
		    IntOutputIterator actuallyClip,
		    const T epsilon) {
  // We will need the edge length when optimizing the mesh below.
  const T edgeLength = geom::computeMaximumEdgeLength(*mesh);

  // For each clipping atom, and while we have not erased the mesh 
  // with clipping.
  for (int i = 0; i != int(clippingAtoms.size()) && 
	 mesh->getSimplicesSize() != 0; ++i) {
    // See if we clip the mesh using that atom.
    if (clipWithCuts(atom, clippingAtoms[i], mesh, epsilon)) {
      // Record the clipping atom.
      *actuallyClip++ = clippingIdentifiers[i];
    }
  }

  // Improve the quality of the mesh.
  optimize(edgeLength, mesh);
}




template<typename T, typename V, typename IS, typename IntOutputIterator>
inline
T
triangulateVisibleSurfaceWithCutClipping
(const Atom<T>& atom,
 std::vector<int>& clippingIdentifiers,
 std::vector< Atom<T> >& clippingAtoms,
 const T edgeLengthSlope,
 const T edgeLengthOffset,
 const int refinementLevel,
 geom::IndSimpSetIncAdj<3,2,true,T,V,IS>* mesh,
 IntOutputIterator actuallyClip,
 const T epsilon) {
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

  clipWithCutClipping(atom, clippingIdentifiers, clippingAtoms,
		      mesh, actuallyClip, epsilon);

  return targetEdgeLength;
}


END_NAMESPACE_MST

// End of file.
