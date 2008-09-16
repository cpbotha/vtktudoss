// -*- C++ -*-

#if !defined(__elc_LagrangianComm_ipp__)
#error This file is an implementation detail of the class LagrangianComm.
#endif

BEGIN_NAMESPACE_ELC

//
// Communication
//


template<int N, typename T>
inline
void
LagrangianComm<N,T>::
sendMesh(const int numNodes, const void* identifiers, const void* positions,
	 const void* velocities, 
	 const int numFaces, const void* connectivities) {
#ifdef ELC_USE_CPP_INTERFACE
  const int commWorldRank = _comm.Get_rank();
#else
  int commWorldRank;
  MPI_Comm_rank(_comm, &commWorldRank);
#endif

  // Record these for waitMesh().
  _numNodes = numNodes;
  _numFaces = numFaces;

  // See if the solid has given us the global node identifiers.
  const int numIdentifiers = (identifiers == 0 ? 0 : numNodes);

  //
  // Compute a bounding box around the boundary nodes.
  //
  BBox domain;
  {
    ads::Array<1,const Point,false> positionsRef(numNodes, positions);
    domain.bound(positionsRef.begin(), positionsRef.end());
  }

  //
  // Determine the point-to-point communication pattern.
  //

  _eulerianProcessors.clear();
  // We don't use the overlapping Eulerian domains, so we use a 
  // trivial output iterator.
  _pointToPoint.solve(domain, 
		      ads::FixedArray<3,int>(commWorldRank,
					     numNodes, numFaces),
		      ads::constructTrivialOutputIterator(),
		      std::back_inserter(_eulerianProcessors));

  //
  // Send the boundary to each relevant Eulerian processor.
  //
  const int numComm = _eulerianProcessors.size();
  _identifierRequests.resize(numComm);
  _positionRequests.resize(numComm);
  _velocityRequests.resize(numComm);
  _connectivityRequests.resize(numComm);
  int proc = 0;
  for (int i = 0; i != numComm; ++i) {
    proc = _eulerianProcessors[i];
#ifdef ELC_USE_CPP_INTERFACE
    _identifierRequests[i] = 
      _comm.Isend(identifiers, numIdentifiers, MPI::INT, proc, TagIdentifiers);
    _positionRequests[i] = 
      _comm.Isend(positions, numNodes * N, _mpiNumber, proc, TagPositions);
    _velocityRequests[i] = 
      _comm.Isend(velocities, numNodes * N, _mpiNumber, proc, TagVelocities);
    _connectivityRequests[i] = 
      _comm.Isend(connectivities, numFaces * N, MPI::INT, proc, TagFaceData);
#else
    MPI_Isend(const_cast<void*>(identifiers), numIdentifiers, MPI_INT, 
	      proc, TagIdentifiers, _comm, &_identifierRequests[i]);
    MPI_Isend(const_cast<void*>(positions), numNodes * N, 
	      _mpiNumber, proc, TagPositions, _comm, &_positionRequests[i]);
    MPI_Isend(const_cast<void*>(velocities), numNodes * N, 
	      _mpiNumber, proc, TagVelocities, _comm, &_velocityRequests[i]);
    MPI_Isend(const_cast<void*>(connectivities), numFaces * N, 
	      MPI_INT, proc, TagFaceData, _comm, &_connectivityRequests[i]);
#endif
  }
}



template<int N, typename T>
inline
void
LagrangianComm<N,T>::
waitForMesh() {
  const int numComm = _eulerianProcessors.size();
  // For each Eulerian processor with which we are communicating.
  for (int i = 0; i != numComm; ++i) {
#ifdef ELC_USE_CPP_INTERFACE
    _identifierRequests[i].Wait();
    _positionRequests[i].Wait();
    _velocityRequests[i].Wait();
    _connectivityRequests[i].Wait();
#else
    MPI_Wait(&_identifierRequests[i], MPI_STATUS_IGNORE);
    MPI_Wait(&_positionRequests[i], MPI_STATUS_IGNORE);
    MPI_Wait(&_velocityRequests[i], MPI_STATUS_IGNORE);
    MPI_Wait(&_connectivityRequests[i], MPI_STATUS_IGNORE);
#endif
  }
}



template<int N, typename T>
inline
void
LagrangianComm<N,T>::
receivePressure(const int numPoints, void* pressures) {
  // Record for wait_pressure().
  _numPoints = numPoints;
  _compositePressures = ads::Array<1,Number,false>(numPoints, pressures);
  _compositePressures = std::numeric_limits<Number>::max();

  //
  // Receive the pressure.
  //
  const int numComm = _eulerianProcessors.size();
  _pressures.resize(numComm);
  _pressureRequests.resize(numComm);
  int proc;
  for (int i = 0; i != numComm; ++i) {
    proc = _eulerianProcessors[i];
    _pressures[i].resize(numPoints);
#ifdef ELC_USE_CPP_INTERFACE
    _pressureRequests[i] = 
      _comm.Irecv(_pressures[i].data(), numPoints, _mpiNumber, proc, 
		  TagPressures);
#else
    MPI_Irecv(_pressures[i].data(), numPoints, _mpiNumber, 
	       proc, TagPressures, _comm, &_pressureRequests[i]);
#endif
  }
}



template<int N, typename T>
inline
void
LagrangianComm<N,T>::
waitForPressure() {
  MpiStatus status;
  const int numComm = _eulerianProcessors.size();
  int proc;
  for (int i = 0; i != numComm; ++i) {
#ifdef ELC_USE_CPP_INTERFACE
    _pressureRequests[i].Wait(status);
    assert(status.Get_count(_mpiNumber) == _numPoints);
#else
    MPI_Wait(&_pressureRequests[i], &status);
    int count;
    MPI_Get_count(&status, _mpiNumber, &count);
    assert(count == _numPoints);
#endif
    // Copy the pressures that were determined by this Lagrangian processor.
    for (int n = 0; n != _numPoints; ++n) {
      // CONTINUE: reference.
      if (_pressures[i][n] != std::numeric_limits<Number>::max()) {
	_compositePressures[n] = _pressures[i][n];
      }
    }
    _pressures[i].resize(0);
  }
}

END_NAMESPACE_ELC
