// -*- C++ -*-

/*!
  \file amr/CellData.h
  \brief An multidimensional array of cell-centered data.
*/

#if !defined(__amr_CellData_h__)
#define __amr_CellData_h__

#include "MessageInputStream.h"
#include "MessageOutputStreamChecked.h"

#include "../array/MultiArray.h"

#include "../third-party/loki/TypeManip.h"

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_CellData)
#define DEBUG_amr_CellData
#endif

BEGIN_NAMESPACE_AMR

//! An multidimensional array of cell-centered data.
/*!
  \note This class stores a multi-array (array::MultiArray) of arrays 
  (array::Array which derives from std::tr1::array). The size of 
  std::tr1::array<T, N> is not necessarily N times the size of T because
  of vector alignment. For certain data types and depths this may waste 
  space. But the alignment rules may improve performance.
*/
template<class _Traits, std::size_t _Depth, std::size_t _GhostWidth, 
	 typename _Number = typename _Traits::Number>
class
CellData {
  //
  // Enumerations.
  //
public:

  // The field depth.
  enum {Depth = _Depth, GhostWidth = _GhostWidth};

  //
  // Public types.
  //
public:

  //! The number type.
  typedef _Number Number;
  //! The tuple of Depth numbers that form the fields.
  typedef array::Array<Number, Depth> FieldTuple;
  //! The array type.
  typedef array::MultiArray<FieldTuple, _Traits::Dimension> Array;
  //! The array view type.
  typedef typename Array::View ArrayView;
  //! The constant array view type.
  typedef typename Array::ConstView ArrayConstView;
  //! A list of sizes.
  typedef typename Array::SizeList SizeList;
  //! A spatial index.
  typedef typename _Traits::SpatialIndex SpatialIndex;


  //
  // Private types.
  //
private:

  //! A list of indices.
  typedef typename _Traits::IndexList IndexList;
  //! An index range.
  typedef typename Array::Range Range;

  //
  // Member data.
  //
private:

  Array _array;

  //
  // Not implemented.
  //
private:

  // Default constructor not implemented.
  CellData();

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Copy constructor.
  CellData(const CellData& other) :
    _array(other._array) {
  }

  //! Allocate memory and initialize the array.
  CellData(const SpatialIndex& spatialIndex, const SizeList& extents,
	   const FieldTuple& initialValues = FieldTuple(Number(0))) :
    // Build the array.
    _array(extents + 2 * GhostWidth,
	   IndexList(extents * spatialIndex.getCoordinates()) - GhostWidth) {
    // Set the initial values.
    _array.fill(initialValues);
  }

  //! Allocate memory and initialize the array.
  /*! The array extents are computed from \c cellData. */
  CellData(const SpatialIndex& spatialIndex, const CellData& cellData,
	   const FieldTuple& initialValues = FieldTuple(Number(0))) :
    // Build the array.
    _array(cellData.getArray().extents(),
	   IndexList(cellData.getInteriorExtents() * 
		     spatialIndex.getCoordinates()) - GhostWidth) {
    // Set the initial values.
    _array.fill(initialValues);
  }

  //! Assignment operator.
  CellData&
  operator=(const CellData& other);

  //! Destructor.
  ~CellData() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Get a const reference to the cell array.
  const Array&
  getArray() const {
    return _array;
  }

  //! Get the index extents for the interior portion of the array.
  /*! The interior portion excludes the ghost region. */
  SizeList
  getInteriorExtents() const {
    return _array.extents() - 2 * GhostWidth;
  }

  //! Get the index range for the interior portion of the array.
  /*! The interior portion excludes the ghost region. */
  Range
  getInteriorRange() const {
    return Range(getInteriorExtents(), _array.bases() + GhostWidth);
  }

  //! Get the index bases for the interior portion of the array.
  /*! The interior portion excludes the ghost region. */
  IndexList
  getInteriorBases() const {
    return _array.bases() + GhostWidth;
  }

  //! Get a constant array that references the interior portion of the array.
  /*! The interior portion excludes the ghost region. */
  ArrayConstView
  getInteriorArray() const {
    return _array.view(getInteriorRange());
  }

  //! Get an array that references the interior portion of the array.
  /*! The interior portion excludes the ghost region. */
  ArrayView
  getInteriorArray() {
    return _array.view(getInteriorRange());
  }

  //! Get the message stream size for this object.
  /*! \pre This must be initialized. */
  int
  getMessageStreamSize() const {
    return _array.size() * sizeof(FieldTuple);
  }

  //! Get the message stream size for this object.
  static
  int
  getMessageStreamSize(SizeList extents) {
    extents += 2 * GhostWidth;
    return product(extents) * sizeof(FieldTuple);
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Get a reference to the cell array.
  Array&
  getArray() {
    return _array;
  }

  //@}
  //--------------------------------------------------------------------------
  //! Prolongation, restriction, and synchronization.
  //@{
public:

  //! Prolongation from interior array data that is one level lower.
  void
  prolong(const CellData& source) {
    prolongConstant(source);
  }

  //! Restriction from interior array data that is one level higher.
  void
  restrict(const CellData& source) {
    restrictLinear(source);
  }

  //! Copy from the interior array data. The patches must be at the same level and adjacent.
  void
  copy(const CellData& source);

private:

  //! Prolongation from interior array data that is one level lower.
  /*! Perform constant extrapolation. */
  void
  prolongConstant(const CellData& source);

  //! Restriction from interior array data that is one level higher.
  /*! Perform averaging. */
  void
  restrictLinear(const CellData& source);

  //@}
  //--------------------------------------------------------------------------
  //! \name Message stream I/O.
  //@{
public:

  //! Write to the message stream.
  void
  write(MessageOutputStream& out) const {
    // Write the elements.
    out.write(_array.data(), _array.size());
  }

  //! Write to the checked message stream.
  void
  write(MessageOutputStreamChecked& out) const {
    // Write the elements.
    out.write(_array.data(), _array.size());
  }

  //! Read from the message stream.
  void
  read(MessageInputStream& in) {
    // Read the elements.
    in.read(_array.data(), _array.size());
  }

  //@}
};

//! Return true if the arrays are equal.
template<class _Traits, std::size_t _Depth, std::size_t _GhostWidth,
	 typename _Number>
inline
bool
operator==(const CellData<_Traits, _Depth, _GhostWidth, _Number>& x,
	   const CellData<_Traits, _Depth, _GhostWidth, _Number>& y) {
  return x.getArray() == y.getArray();
}

//! Write the cell data as an array.
/*!
  \relates CellData
*/
template<class _Traits, std::size_t _Depth, std::size_t _GhostWidth,
	 typename _Number>
inline
std::ostream&
operator<<(std::ostream& out,
	   const amr::CellData<_Traits, _Depth, _GhostWidth, _Number>& x) {
  // Write the array.
  out << x.getArray() << "\n";
  return out;
}

END_NAMESPACE_AMR

#define __amr_CellData_ipp__
#include "CellData.ipp"
#undef __amr_CellData_ipp__

#endif
