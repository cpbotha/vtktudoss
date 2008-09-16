// -*- C++ -*-

#if !defined(__elc_EulerianComm_ipp__)
#error This file is an implementation detail of the class EulerianComm.
#endif

BEGIN_NAMESPACE_ELC


//
// Communication
//


template<int N, typename T>
inline
void
EulerianComm<N,T>::
receiveMesh(const BBox& domain) {
#ifdef ELC_USE_CPP_INTERFACE
  const int commWorldRank = _comm.Get_rank();
#else
  int commWorldRank;
  MPI_Comm_rank(_comm, &commWorldRank);
#endif

  //
  // Determine the point-to-point communication pattern.
  //

  _lagrangianData.clear();
  // We don't use the overlapping Lagrangian domains, so we use a 
  // trivial output iterator.
  _pointToPoint.solve(domain, commWorldRank,
		      ads::constructTrivialOutputIterator(),
		      std::back_inserter(_lagrangianData));

  //
  // Resize the vectors of boundary information.
  //

  const int numComm = _lagrangianData.size();
  _identifiers.resize(numComm);
  _positions.resize(numComm);
  _velocities.resize(numComm);
  _connectivities.resize(numComm);

  //
  // Receive the boundary from each relevant Lagrangian processor.
  //

  _identifierRequests.resize(numComm);
  _positionRequests.resize(numComm);
  _velocityRequests.resize(numComm);
  _connectivityRequests.resize(numComm);

  for (int i = 0; i != numComm; ++i) {
    const int lagrangianProc = _lagrangianData[i][0];
    const int numNodesToReceive = _lagrangianData[i][1];
    const int numFacesToReceive = _lagrangianData[i][2];

    _identifiers[i].resize(numNodesToReceive);
    _positions[i].resize(numNodesToReceive);
    _velocities[i].resize(numNodesToReceive);
    _connectivities[i].resize(numFacesToReceive);
#ifdef ELC_USE_CPP_INTERFACE
    _identifierRequests[i] =
      _comm.Irecv(_identifiers[i].data(), numNodesToReceive, MPI::INT,
		   lagrangianProc, TagIdentifiers);

    _positionRequests[i] =
      _comm.Irecv(_positions[i].data(), numNodesToReceive * N, _mpiNumber, 
		   lagrangianProc, TagPositions);

    _velocityRequests[i] = 
      _comm.Irecv(_velocities[i].data(), numNodesToReceive * N, _mpiNumber,
		   lagrangianProc, TagVelocities);

    _connectivityRequests[i] =
      _comm.Irecv(_connectivities[i].data(), numFacesToReceive * N, MPI::INT, 
		   lagrangianProc, TagFaceData);
#else
    MPI_Irecv(_identifiers[i].data(), numNodesToReceive, MPI_INT, 
	       lagrangianProc, TagIdentifiers, _comm, 
	       &_identifierRequests[i]);

    MPI_Irecv(_positions[i].data(), numNodesToReceive * N, _mpiNumber, 
	       lagrangianProc, TagPositions, _comm, &_positionRequests[i]);

    MPI_Irecv(_velocities[i].data(), numNodesToReceive * N, _mpiNumber, 
	       lagrangianProc, TagVelocities, _comm, 
	       &_velocityRequests[i]);

    MPI_Irecv(_connectivities[i].data(), numFacesToReceive * N, MPI_INT, 
	       lagrangianProc, TagFaceData, _comm, 
	       &_connectivityRequests[i]);
#endif
  }
}



template<int N, typename T>
inline
void
EulerianComm<N,T>::
waitForMesh() {
  //
  // Receive the identifiers and the positions.
  // If the solid does not provide global identifiers, the positions are 
  // needed to construct these identifiers.
  //
  bool shouldGenerateIdentifiers = true;
  MpiStatus status;
  const int numComm = _lagrangianData.size();
  for (int i = 0; i != numComm; ++i) {
    const int numNodesToReceive = _lagrangianData[i][1];
    int count;
    // Wait for the identifier receives to complete.
#ifdef ELC_USE_CPP_INTERFACE
    _identifierRequests[i].Wait(status);
    count = status.Get_count(MPI::INT);
#else
    MPI_Wait(&_identifierRequests[i], &status);
    MPI_Get_count(&status, MPI_INT, &count);
#endif
    if (count != 0) {
      shouldGenerateIdentifiers = false;
    }
    if (shouldGenerateIdentifiers) {
      assert(count == 0);
    }
    else {
      assert(count == numNodesToReceive);
    }

    // Wait for the position receives to complete.
#ifdef ELC_USE_CPP_INTERFACE
    _positionRequests[i].Wait(status);
    assert(status.Get_count(_mpiNumber) == numNodesToReceive * N);
#else
    MPI_Wait(&_positionRequests[i], &status);
    MPI_Get_count(&status, _mpiNumber, &count);
    assert(count == numNodesToReceive * N);
#endif

    // Wait for the velocity receives to complete.
#ifdef ELC_USE_CPP_INTERFACE
    _velocityRequests[i].Wait(status);
    assert(status.Get_count(_mpiNumber) == numNodesToReceive * N);
#else
    MPI_Wait(&_velocityRequests[i], &status);
    MPI_Get_count(&status, _mpiNumber, &count);
    assert(count == numNodesToReceive * N);
#endif

    // Wait for the connectivity receives to complete.
    const int numFacesSent = _lagrangianData[i][2];
#ifdef ELC_USE_CPP_INTERFACE
    _connectivityRequests[i].Wait(status);
    assert(status.Get_count(MPI::INT) == numFacesSent * N);
#else
    MPI_Wait(&_connectivityRequests[i], &status);
    MPI_Get_count(&status, MPI_INT, &count);
    assert(count == numFacesSent * N);
#endif
  }

#if 0
  // Disabled.
  if (shouldGenerateIdentifiers) {
    generateIdentifiers();
  }
#endif

  //
  // Build the mapping from node identifiers to node indices in the assembled
  // boundary.
  //

  _identifierToIndex.clear();
  std::pair<std::map<int,int>::iterator, bool> insertResult;
  // Loop over the boundary portions.
  for (int i = 0; i != numComm; ++i) {
    // Loop over nodes.
    for (int n = 0; n != _identifiers[i].size(); ++n) {
      // If the identifier is not in the mapping.
      if (! _identifierToIndex.count(_identifiers[i][n])) {
	// Add it to the mapping.
	std::map<int,int>::value_type 
	  insertValue(_identifiers[i][n], _identifierToIndex.size());
	insertResult = _identifierToIndex.insert(insertValue);
	// Make sure that it was inserted.
	assert(insertResult.second);
      }
    }
  }

  //
  // Build the assembled boundary.
  //

  // Determine the number of nodes and faces.
  const int numAssembledNodes = _identifierToIndex.size();
  int numAssembledFaces = 0;
  for (std::size_t i = 0; i != _connectivities.size(); ++i) {
    numAssembledFaces += _connectivities[i].size();
  }

  // Allocate memory.
  _assembledPositions.resize(numAssembledNodes);
  _assembledVelocities.resize(numAssembledNodes);
  _assembledConnectivities.resize(numAssembledFaces);

  // Set the positions and velocities.
  {
    std::map<int,int>::const_iterator pairIterator;
    int index;
    for (int i = 0; i != numComm; ++i) {
      for (int n = 0; n != _positions[i].size(); ++n) {
	// Find the identifier.
	pairIterator = _identifierToIndex.find(_identifiers[i][n]);
	// Make sure that we found it.
	assert(pairIterator != _identifierToIndex.end());
	// Extract the index.
	index = pairIterator->second;
	// Set the position and velocity.
	_assembledPositions[index] = _positions[i][n];
	_assembledVelocities[index] = _velocities[i][n];
      }

      // Free the position and velocity memory for the i_th processor.
      {
	// I do this to placate the xlC compiler on frost.
	int size = 0;
	_positions[i].resize(size);
	_velocities[i].resize(size);
      }
      // Don't free the identifier memory.  We'll need that for sending
      // the pressure.
    }
  }

  // Set the connectivities.
  {
    std::map<int,int>::const_iterator pairIterator;
    int cellIndex = 0;
    int localIndex, globalIdentifier;
    for (int i = 0; i != numComm; ++i) {
      for (int j = 0; j != _connectivities[i].size(); ++j) {
	for (int n = 0; n != N; ++n) {
	  if (_vertexIdentifierStyle == LocalIndices) {
	    // The local node index in the i_th Lagrangrian processors 
	    // with which we communicate.
	    localIndex = _connectivities[i][j][n];
	    // Make sure the local index is in the right range.
	    assert(0 <= localIndex && localIndex < _identifiers[i].size());
	    // Switch from the local Lagrangian identifiers to the global
	    // Lagrangian identifiers.
	    globalIdentifier = _identifiers[i][localIndex];
	  }
	  else {
	    assert(_vertexIdentifierStyle == GlobalIdentifiers);
	    globalIdentifier = _connectivities[i][j][n];
	  }
	  // Find the identifier.
	  pairIterator = _identifierToIndex.find(globalIdentifier);
	  // Make sure that we found it.
	  assert(pairIterator != _identifierToIndex.end());
	  // Extract the node index.
	  _assembledConnectivities[cellIndex][n] = pairIterator->second;
	}
	++cellIndex;
      }
      // Free the connectivities memory for the i_th processor.
      _connectivities[i].resize(0);
    }
    // Sanity check.
    assert(cellIndex == _assembledConnectivities.size());
  }
  
  initializePressure();
}



template<int N, typename T>
inline
void
EulerianComm<N,T>::
generateIdentifiers() {
  //
  // Initially assign each node an identifier.
  //
  const int numComm = _lagrangianData.size();
  int globalIdentifier = 0;
  // Loop over the boundary portions.
  for (int i = 0; i != numComm; ++i) {
    // Loop over nodes in this patch.
    for (int n = 0; n != _identifiers[i].size(); ++n) {
      _identifiers[i][n] = globalIdentifier++;
    }
  }
    
  //
  // Now find the nodes that are on the boundary of each patch.
  //
  std::vector<std::pair<int*, const Point*> > boundaryNodes;
  std::pair<int*, const Point*> value;
  std::set<int> boundaryIndices;
  // Loop over the boundary portions.
  for (int i = 0; i != numComm; ++i) {
    // Get the boundary nodes from this patch.
    geom::IndSimpSetIncAdj<N,N-1,false,Number> 
      patch(_positions[i].size(), _positions[i].data(),
	    _connectivities[i].size(), _connectivities[i].data());
    
    geom::determineBoundaryVertices
      (patch, std::inserter(boundaryIndices, boundaryIndices.end()));
    for (std::set<int>::const_iterator iter = boundaryIndices.begin();
	 iter != boundaryIndices.end(); ++iter) {
      value.first = &_identifiers[i][*iter];
      value.second = &_positions[i][*iter];
      boundaryNodes.push_back(value);
    }
    boundaryIndices.clear();
  }

  //
  // Sort the boundary nodes.  Use the position as a composite number.
  //
  // CONTINUE: An ORQ would have better computational complexity.
  {
    // Functor that selects the second of the pair and then dereferences.
    typedef 
      ads::unary_compose_unary_unary<ads::Dereference<const Point*>, 
      ads::Select2nd< std::pair<int*, const Point*> > > GetPoint;
    // Functor that compares the std::pairs.
    typedef
      ads::binary_compose_binary_unary<ads::less_composite<N,Point>,
      GetPoint, GetPoint> Comp;

    // Set the x coordinate as the first in the composite number.
    ads::less_composite<N,Point> lc;
    lc.set(0);

    // Make the comparison functor.
    GetPoint gp;
    Comp comp(lc, gp, gp);
      
    // Sort the nodes.
    std::sort(boundaryNodes.begin(), boundaryNodes.end(), comp);
  }

  //
  // Detect which nodes coincide and fix those identifiers.
  //
    
  // Any nodes closer than epsilon to each other will be considered the 
  // same node.
  // CONTINUE: Use the bounding box of all of the nodes.
  const Number epsilon = 
    std::sqrt(std::numeric_limits<Number>::epsilon());
  const Number epsilonSquared = epsilon * epsilon;

  // Loop over the boundary nodes.
  const int size = boundaryNodes.size();
  Number x;
  int m;
  for (int n = 1; n < size; ++n) {
    // Look for nodes that coincide.
    x = (*boundaryNodes[n].second)[0];
    m = n - 1;
    while (m >= 0 && x - (*boundaryNodes[m].second)[0] < epsilon) {
      if (geom::computeSquaredDistance(*boundaryNodes[m].second,
				       *boundaryNodes[n].second) <
	  epsilonSquared) {
	// Fix the identifier.
	*boundaryNodes[n].first = *boundaryNodes[m].first;
	break;
      }
      --m;
    }
  }
}


namespace internal {

  template<typename T>
  inline
  void
  computeFaceNormal(const ads::Array<1, ads::FixedArray<2,T> >& pos,
		    const ads::Array<1, ads::FixedArray<2,int> >& conn,
		    const int n, ads::FixedArray<2,T>* normal,
		    Loki::Int2Type<2> /* dummy */) {
    ads::FixedArray<2,T> tangent(pos[conn[n][1]]);
    tangent -= pos[conn[n][0]];
    geom::normalize(&tangent);
    (*normal)[0] = tangent[1];
    (*normal)[1] = - tangent[0];
  }

  template<typename T>
  inline
  void
  computeFaceNormal(const ads::Array<1, ads::FixedArray<3,T> >& pos,
		    const ads::Array<1, ads::FixedArray<3,int> >& conn,
		    const int n, ads::FixedArray<3,T>* normal,
		    Loki::Int2Type<3> /* dummy */) {
    const ads::FixedArray<3,int>& face = conn[n];
    geom::computeCrossProduct(pos[face[2]] - pos[face[1]],
			      pos[face[0]] - pos[face[1]], normal);
    geom::normalize(normal);
  }

  template<int N, typename T>
  inline
  void
  computeFaceNormal(const ads::Array<1, ads::FixedArray<N,T> >& pos,
		    const ads::Array<1, ads::FixedArray<N,int> >& conn,
		    const int n, ads::FixedArray<N,T>* normal) {
    computeFaceNormal(pos, conn, n, normal, Loki::Int2Type<N>());
  }

}


// Compute the face normals.
template<int N, typename T>
inline
void
EulerianComm<N,T>::
computeFaceNormals() {
  // Resize the array.
  _faceNormals.resize(getNumberOfFaces());

  // Loop over the faces.
  for (int i = 0; i != getNumberOfFaces(); ++i) {
    // Compute the face normal.
    computeFaceNormal(i, &_faceNormals[i]);
  }
}

  
// Compute the face centroids.
template<int N, typename T>
inline
void
EulerianComm<N,T>::
computeFaceCentroids() {
  // Resize the array.
  _faceCentroids.resize(getNumberOfFaces());

  // Loop over the faces.
  for (int i = 0; i != getNumberOfFaces(); ++i) {
    // Compute the face centroid.
    computeFaceCentroid(i, &_faceCentroids[i]);
  }
}

  
template<int N, typename T>
inline
void
EulerianComm<N,T>::
computeFaceNormal(const int n, Point* normal) const {
  // Get the face normal.
  internal::computeFaceNormal<N>(_assembledPositions, 
				 _assembledConnectivities, n, normal);
}


template<int N, typename T>
inline
void
EulerianComm<N,T>::
computeFaceCentroid(const int n, Point* centroid) const {
  // The centroid is the arithmetic mean of the face node positions.

  *centroid = 0.0;
  // CONTINUE
  //assert(0 <= n && n < _assembledConnectivities.size());
  const IndexedFace& face = _assembledConnectivities[n];
  // Loop over the nodes of the face.
  for (int j = 0; j != N; ++j) {
    //centroid += _assembledPositions[ _assembledConnectivities[n][j] ];
    *centroid += _assembledPositions[face[j]];
  }
  *centroid /= N;
}

END_NAMESPACE_ELC
