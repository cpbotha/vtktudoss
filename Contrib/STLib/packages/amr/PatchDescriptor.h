// -*- C++ -*-

/*!
  \file amr/PatchDescriptor.h
  \brief Describe a patch.
*/

#if !defined(__amr_PatchDescriptor_h__)
#define __amr_PatchDescriptor_h__

#include "FieldDescriptor.h"

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_PatchDescriptor)
#define DEBUG_amr_PatchDescriptor
#endif

BEGIN_NAMESPACE_AMR

//! Describe a patch.
/*!
  \param _Traits Traits for the orthtree.
*/
template<typename _Traits>
class
PatchDescriptor {
  //
  // Public types.
  //
public:

  //! The spatial index.
  typedef typename _Traits::SpatialIndex SpatialIndex;
  //! A multi-index.
  typedef array::Array<int, _Traits::Dimension> Index;

  //
  // Member data.
  //
private:

  //! The index extents of a patch.
  Index _extents;
  //! The ghost cell width.
  int _ghostCellWidth;
  //! The field descriptor.
  FieldDescriptor _fieldDescriptor;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
private:

  //! The default constructor is not implemented.
  PatchDescriptor();

public:

  //! Construct from the array extents and the ghost cell width.
  PatchDescriptor(const Index& extents, const int ghostCellWidth,
		  const FieldDescriptor& fieldDescriptor) :
    _extents(extents),
    _ghostCellWidth(ghostCellWidth),
    _fieldDescriptor(fieldDescriptor) {
  }

  //! Copy constructor.
  PatchDescriptor(const PatchDescriptor& other) :
    _extents(other._extents),
    _ghostCellWidth(other._ghostCellWidth),
    _fieldDescriptor(other._fieldDescriptor) {
  }

  //! Assignment operator.
  PatchDescriptor&
  operator=(const PatchDescriptor& other) {
    if (this != &other) {
      _extents = other._extents;
      _ghostCellWidth = other._ghostCellWidth;
      _fieldDescriptor = other._fieldDescriptor;
    }
    return *this;
  }

  //! Destructor.
  ~PatchDescriptor() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Get the grid index extents.
  const Index&
  getExtents() const {
    return _extents;
  }

  //! Return the ghost cell width.
  int
  getGhostCellWidth() const {
    return _ghostCellWidth;
  }

  //! Get the field descriptor.
  const FieldDescriptor&
  getFieldDescriptor() const {
    return _fieldDescriptor;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
public:

  //! Return true if the data structures are equal.
  bool
  operator==(const PatchDescriptor& other) {
    return _extents == other._extents && 
      _ghostCellWidth == other._ghostCellWidth;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Mathematical.
  //@{
public:

  //! Compute the index lower bounds for the patch.
  Index
  computeLowerBounds(const SpatialIndex& spatialIndex) const {
    return _extents * spatialIndex.getCoordinates();
  }

  //@}
};

END_NAMESPACE_AMR

#define __amr_PatchDescriptor_ipp__
#include "PatchDescriptor.ipp"
#undef __amr_PatchDescriptor_ipp__

#endif
