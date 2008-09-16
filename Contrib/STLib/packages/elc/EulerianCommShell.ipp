// -*- C++ -*-

#if !defined(__elc_EulerianCommShell_ipp__)
#error This file is an implementation detail of the class EulerianCommShell.
#endif

BEGIN_NAMESPACE_ELC


template<int N, typename T>
inline
void
EulerianCommShell<N,T>::
sendPressure() {
  //
  // Send the pressures.
  //
  const int numComm = _lagrangianData.size();
  _pressureRequests.resize(numComm);
  int faceIndex = 0;
  // For each Lagrangian processor with which we are communicating.
  for (int i = 0; i != numComm; ++i) {
    // We free these here to match the boundary case.
    _identifiers[i].resize(0);
    // The Lagrangian processor rank in the world communicator.
    const int lagrangianProcessor = _lagrangianData[i][0];
    // The pressures are defined at the faces.  Thus the number of 
    // pressures is the same as the number of faces.
    assert(_pressures[i].size() == _lagrangianData[i][2]);
    // Copy the pressures from the assembled mesh.
    for (int j = 0; j != _pressures[i].size(); ++j) {
      // Extract the pressure.
      _pressures[i][j] = _assembledPressures[faceIndex];
      ++faceIndex;
    }
#ifdef ELC_USE_CPP_INTERFACE
    _pressureRequests[i] = 
      _comm.Isend(_pressures[i].data(), _pressures[i].size(),
		  _mpiNumber, lagrangianProcessor, TagPressures);
#else
    MPI_Isend(_pressures[i].data(), _pressures[i].size(),
	      _mpiNumber, lagrangianProcessor, TagPressures, _comm, 
	      &_pressureRequests[i]);
#endif
  }
}



template<int N, typename T>
inline
void
EulerianCommShell<N,T>::
waitForPressure() {
  const int numComm = _lagrangianData.size();
  // For each Lagrangian processor with which we are communicating.
  for (int i = 0; i != numComm; ++i) {
#ifdef ELC_USE_CPP_INTERFACE
    _pressureRequests[i].Wait();
#else
    MPI_Wait(&_pressureRequests[i], MPI_STATUS_IGNORE);
#endif
    // Free the pressure memory for the i_th processor.
    _pressures[i].resize(0);
  }
}



template<int N, typename T>
inline
void
EulerianCommShell<N,T>::
initializePressure() {
  const int numComm = _lagrangianData.size();
  _pressures.resize(numComm);
  for (int i = 0; i != numComm; ++i) {
    const int numFacesToReceive = _lagrangianData[i][2];
    _pressures[i].resize(numFacesToReceive);
  }

  // Allocate memory.
  _assembledPressures.resize(_assembledConnectivities.size());

  // Fill the pressures with a flag value.
  _assembledPressures = std::numeric_limits<Number>::max();
}

END_NAMESPACE_ELC
