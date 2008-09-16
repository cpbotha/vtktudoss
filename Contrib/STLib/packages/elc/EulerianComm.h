// -*- C++ -*-

/*! 
  \file EulerianComm.h
  \brief Eulerian communicator for the Eulerian-Lagrangian coupling.
*/

#if !defined(__elc_EulerianComm_h__)
#define __elc_EulerianComm_h__

#include "ELComm.h"

#if defined(ELC_USE_CPP_INTERFACE) && !defined(PT2PT_BBOX_USE_CPP_INTERFACE)
#define PT2PT_BBOX_USE_CPP_INTERFACE
#endif
#include "../concurrent/pt2pt_bbox.h"

#include "../ads/iterator/TrivialOutputIterator.h"
#include "../ads/functor/select.h"

#include "../geom/mesh/iss/set.h"

#include "../third-party/loki/TypeManip.h"

#include <set>
#include <map>

// If we are debugging the whole elc package.
#if defined(DEBUG_elc) && !defined(DEBUG_EulerianComm)
#define DEBUG_EulerianComm
#endif

#ifdef DEBUG_EulerianComm
#ifndef DEBUG_ELComm
#define DEBUG_ELComm
#endif
#endif

BEGIN_NAMESPACE_ELC

//! Base class Eulerian communicator for the Eulerian-Lagrangian coupling.
/*!
  \param N is the space dimension.  2 and 3 are supported.
  \param T is the floating point number type.

  Implements the common functionality for boundaries and shells.
*/
template<int N, typename T>
class EulerianComm :
  public ELComm<N,T> {
  //
  // Private types.
  //

private:

  typedef ELComm<N,T> Base;

  //
  // Protected types.
  //

protected:

  //! The number type.
  typedef typename Base::Number Number;
  //! A Cartesian point.
  typedef typename Base::Point Point;
  //! A bounding box.
  typedef typename Base::BBox BBox;
  //! An MPI request.
  typedef typename Base::MpiRequest MpiRequest;
  //! Status for an MPI request.
  typedef typename Base::MpiStatus MpiStatus;
  //! An indexed face type.
  typedef ads::FixedArray<N,int> IndexedFace;

  //
  // Using base members.
  //

protected:

  //! The joint Eulerian/Lagrangian communicator.
  using Base::_comm;
  //! The MPI number type.
  using Base::_mpiNumber;
  //! The vertex identifier style.
  using Base::_vertexIdentifierStyle;

  //! The communication tag for pressures.
  using Base::TagPressures;

private:

  //! The communication tag for node identifiers.
  using Base::TagIdentifiers;
  //! The communication tag for node positions.
  using Base::TagPositions;
  //! The communication tag for node velocities.
  using Base::TagVelocities;
  //! The communication tag for face data.
  using Base::TagFaceData;

  //
  // Member data.
  //

protected:

#ifdef ELC_USE_CPP_INTERFACE
  //! The Eulerian communicator.
  MPI::Intracomm _eulerianCommunicator;
#else
  //! The Eulerian communicator.
  MPI_Comm _eulerianCommunicator;
#endif

  //! The Lagrangian root.
  int _lagrangianRoot;

  //! The node identifiers from Lagrangian processors.
  std::vector<ads::Array<1,int> > _identifiers;
  //! The node positions from Lagrangian processors.
  std::vector<ads::Array<1,Point> > _positions;
  //! The node velocities from Lagrangian processors.
  std::vector<ads::Array<1,Point> > _velocities;
  //! The node connectivities from Lagrangian processors.
  std::vector<ads::Array<1,IndexedFace> > _connectivities;
  //! The node pressures to be sent to Lagrangian processors.
  std::vector<ads::Array<1,Number> > _pressures;


  //! The mapping from node identifiers to node indices in the assembled boundary.
  std::map<int,int> _identifierToIndex;

  // CONTINUE: use something derived from an indexed simplex set instead.
  // That would facilitate interpolation.
  //! The assembled positions.
  ads::Array<1,Point> _assembledPositions;
  //! The assembled velocities.
  ads::Array<1,Point> _assembledVelocities;
  //! The assembled connectivities.
  ads::Array<1,IndexedFace> _assembledConnectivities;
  //! The assembled pressures.
  ads::Array<1,Number> _assembledPressures;

  //! Class for computing the point-to-point communication scheme.
  concurrent::PtToPt2Grp1Dom<N, T, int, ads::FixedArray<3,int> > 
  _pointToPoint;

  //! Data from the Lagrangian processors with which we communicate.
  std::vector<ads::FixedArray<3,int> > _lagrangianData;

  //! The face normals.
  ads::Array<1,Point> _faceNormals;
  //! The face centroids.
  ads::Array<1,Point> _faceCentroids;

private:

  std::vector<MpiRequest> _identifierRequests;
  std::vector<MpiRequest> _positionRequests;
  std::vector<MpiRequest> _velocityRequests;
  std::vector<MpiRequest> _connectivityRequests;

  //
  // Not implemented.
  //

private:

  // Default constructor not implemented.
  EulerianComm();

  // Copy constructor not implemented.
  EulerianComm(const EulerianComm&);

  // Assignment operator not implemented.
  EulerianComm&
  operator=(const EulerianComm&);

public:

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  // @{

#ifdef ELC_USE_CPP_INTERFACE
  //! Construct from the communicators and Lagrangian info.
  /*!
    \param comm is the communicator that contains the Eulerian and 
    Lagrangian processors.
    \param eulerian is the Eulerian communicator.  It is duplicated to avoid
    message conflicts.
    \param lagrangian_size is the number of Lagrangian processors.
    \param lagrangian_root is the rank of the Lagrangian root in \c comm.
    \param vertexIdentifierStyle is either LocalIndices or GlobalIdentifiers.
  */
  EulerianComm(const MPI::Comm& comm, const MPI::Intracomm& eulerian, 
	       const int lagrangianSize, const int lagrangianRoot,
	       VertexIdentifierStyle vertexIdentifierStyle) :
    Base(comm, vertexIdentifierStyle),
    _eulerianCommunicator(eulerian.Dup()),
    _lagrangianRoot(lagrangianRoot),
    _pointToPoint(comm, eulerian, lagrangianSize, lagrangianRoot)
  {}
#else
  //! Construct from the communicators and Lagrangian info.
  /*!
    \param comm is the communicator that contains the Eulerian and 
    Lagrangian processors.
    \param eulerian is the Eulerian communicator.  It is duplicated to avoid
    message conflicts.
    \param lagrangianSize is the number of Lagrangian processors.
    \param lagrangianRoot is the rank of the Lagrangian root in \c comm.
    \param vertexIdentifierStyle is either LocalIndices or GlobalIdentifiers.
  */
  EulerianComm(const MPI_Comm comm, const MPI_Comm eulerian, 
	       const int lagrangianSize, const int lagrangianRoot,
	       VertexIdentifierStyle vertexIdentifierStyle) :
    Base(comm, vertexIdentifierStyle),
    _eulerianCommunicator(),
    _lagrangianRoot(lagrangianRoot),
    _pointToPoint(comm, eulerian, lagrangianSize, lagrangianRoot) {
    MPI_Comm_dup(eulerian, &_eulerianCommunicator);
  }
#endif

  //! Destructor.  Free the duplicated communicator.
  virtual
  ~EulerianComm() {
#ifdef ELC_USE_CPP_INTERFACE
    _eulerianCommunicator.Free();
#else
    MPI_Comm_free(&_eulerianCommunicator);
#endif
  }

  // @}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  // @{

  //! Return the number of nodes.
  int
  getNumberOfNodes() const {
    return _assembledPositions.size();
  }

  //! Return the number of faces.
  int
  getNumberOfFaces() const {
    return _assembledConnectivities.size();
  }

  //! Return a const reference to the array of node positions.
  const ads::Array<1,Point>&
  getPositions() const {
    return _assembledPositions;
  }

  //! Return a const pointer to the node positions data.
  const Number*
  getPositionsData() const {
    return reinterpret_cast<const Number*>(_assembledPositions.data());
  }

  //! Return a const reference to the array of node velocities.
  const ads::Array<1,Point>&
  getVelocities() const {
    return _assembledVelocities;
  }

  //! Return a const pointer to the node velocities data.
  const Number*
  getVelocitiesData() const {
    return reinterpret_cast<const Number*>(_assembledVelocities.data());
  }

  //! Return a const reference to the array of node connectivities.
  const ads::Array<1,IndexedFace>&
  getConnectivities() const {
    return _assembledConnectivities;
  }

  //! Return a const pointer to the node connectivities data.
  const int*
  getConnectivitiesData() const {
    return reinterpret_cast<const int*>(_assembledConnectivities.data());
  }

  // @}
  //--------------------------------------------------------------------------
  //! \name Face Centroid/Normal Accessors.
  // @{

  //! Return a const reference to the array of face normals.
  const ads::Array<1,Point>&
  getFaceNormals() const {
    return _faceNormals;
  }

  //! Return a const pointer to the face normals data.
  const Number*
  getFaceNormalsData() const {
    return reinterpret_cast<const Number*>(_faceNormals.data());
  }

  //! Return a const reference to the n_th face normal.
  const Point&
  getFaceNormal(const int n) const {
    return _faceNormals[n];
  }

  //! Return a const pointer to the n_th face normal data.
  const Number*
  getFaceNormalData(const int n) const {
    return reinterpret_cast<const Number*>(_faceNormals[n].data());
  }

  //! Return a const reference to the array of face centroids.
  const ads::Array<1,Point>&
  getFaceCentroids() const {
    return _faceCentroids;
  }

  //! Return a const pointer to the face centroids data.
  const Number*
  getFaceCentroidsData() const {
    return reinterpret_cast<const Number*>(_faceCentroids.data());
  }

  //! Return a const reference to the n_th face centroid.
  const Point&
  getFaceCentroid(const int n) const {
    return _faceCentroids[n];
  }

  //! Return a const pointer to the n_th face centroid data.
  const Number*
  getFaceCentroidData(const int n) const {
    return reinterpret_cast<const Number*>(_faceCentroids[n].data());
  }

  // @}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  // @{

  //! Return a reference to the array of node pressures.
  ads::Array<1,Number>&
  getPressures() {
    return _assembledPressures;
  }

  //! Return a pointer to the node pressures data.
  Number*
  getPressuresData() {
    return reinterpret_cast<Number*>(_assembledPressures.data());
  }

  // @}
  //--------------------------------------------------------------------------
  //! \name Communication.
  // @{

  //! Post receives for the relevant portions of the mesh from the solid processors.
  /*!
    \param domain is the region of interest for this fluid processor.
  */
  void
  receiveMesh(const BBox& domain);

  //! Post receives for the relevant portions of the mesh from the solid processors.
  /*!
    \param domain is the region of interest for this fluid processor.

    This function just calls the above
    \c receive_mesh(const BBox& domain).
  */
  void
  receiveMesh(const Number* domain) {
    receive_mesh(BBox(domain));
  }

  //! Wait for the receives to complete.  Build the assembled mesh.
  /*!
    This function must be called before accessing the mesh.
  */
  void
  waitForMesh();


  //! Start sending the pressure to the relevant solid processors.
  /*!
    Call this function after the pressures have been set.
  */
  virtual
  void
  sendPressure() = 0;

  //! Wait for the pressure sends to be copied into communication buffers.
  /*!
    This function must be called after sendPressure().
  */
  virtual
  void
  waitForPressure() = 0;


  //! Compute the face normals.
  /*!
    Call this after waitForMesh() if you will use the face normals.
  */
  void
  computeFaceNormals();
  
  //! Compute the face centroids.
  /*!
    Call this after waitForMesh() if you will use the face centroids.
  */
  void
  computeFaceCentroids();
  
  //! Initialize the pressure at the nodes or the faces.
  /*!
    This function should be called after waitForMesh() and before accessing 
    the pressures.
    
    This function is implemented in the derived classes.
  */
  virtual
  void
  initializePressure() = 0;

  // @}

private:

  //! Compute the n_th face normal.
  void
  computeFaceNormal(const int n, Point* normal) const;

  //! Compute the n_th face centroid.
  void
  computeFaceCentroid(const int n, Point* centroid) const;

  //! Generate identifiers for the pieces of the boundary.
  void
  generateIdentifiers();

};


END_NAMESPACE_ELC

#define __elc_EulerianComm_ipp__
#include "EulerianComm.ipp"
#undef __elc_EulerianComm_ipp__

#endif
