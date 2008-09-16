// -*- C++ -*-

/*!
  \file array/MultiArrayBase.h
  \brief Base class for multi-arrays.
*/

#if !defined(__array_MultiArrayBase_h__)
#define __array_MultiArrayBase_h__

#include "IndexTypes.h"
#include "IndexRange.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_MultiArrayBase)
#define DEBUG_array_MultiArrayBase
#endif

BEGIN_NAMESPACE_ARRAY

//! Base class for multi-arrays.
template<std::size_t _Dimension>
class
MultiArrayBase {
  //
  // Enumerations.
  //
public:

  //! The number of dimensions.
  enum {Dimension = _Dimension};

  //
  // Types.
  //
private:

  typedef IndexTypes<_Dimension> Types;

public:

  // Types for STL compliance.

  //! The size type.
  typedef typename Types::size_type size_type;
  //! Pointer difference type.
  typedef typename Types::difference_type difference_type;

  // Other types.

  //! An array index is a signed integer.
  typedef typename Types::Index Index;
  //! A list of indices.
  typedef typename Types::IndexList IndexList;
  //! A list of sizes.
  typedef typename Types::SizeList SizeList;
  //! The storage order.
  typedef typename Types::Storage Storage;
  //! A multi-index range.
  typedef IndexRange<Dimension> Range;

  //
  // Member data.
  //
protected:

  //! The array extents.
  SizeList _extents;
  //! The lower bound for each index.
  IndexList _bases;
  //! The storage order (from least to most significant).
  Storage _storage;
  //! The strides for indexing.
  IndexList _strides;
  //! The offset for indexing the bases.
  Index _offset;
  //! The number of elements.
  size_type _size;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  // The default copy constructor and assignment operator are fine.

  //! Construct from the array extents, the index bases, the storage order, and the strides.
  MultiArrayBase(const SizeList& extents, const IndexList& bases,
		 const Storage& storage, const IndexList& strides);

  //! Destructor does nothing.
  virtual
  ~MultiArrayBase() {
  }

protected:

  //! Rebuild the data structure.
  void
  rebuild(const SizeList& extents, const IndexList& bases,
	  const Storage& storage, const IndexList& strides);

private:

  //! Default constructor not implemented.
  MultiArrayBase();

  //@}
  //--------------------------------------------------------------------------
  //! \name Random Access Container.
  //@{
public:

  //! Return true if the range is empty.
  bool
  empty() const {
    return _size == 0;
  }

  //! Return the size (number of elements) of the range.
  size_type
  size() const {
    return _size;
  }

  //! Return the size of the range.
  /*! The the max_size and the size are the same. */
  size_type
  max_size() const {
    return size();
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Array indexing.
  //@{
public:

  //! The range extents.
  const SizeList&
  extents() const {
    return _extents;
  }

  //! The lower bound for each index.
  const IndexList&
  bases() const {
    return _bases;
  }

  //! Set the lower bounds for each index.
  void
  setBases(const IndexList& bases) {
    _bases = bases;
    _offset = dot(_strides, _bases);
  }

  //! The index ranges.
  Range
  range() const {
    return Range(extents(), bases());
  }

  //! The storage order.
  const Storage&
  storage() const {
    return _storage;
  }

  //! The strides for indexing.
  const IndexList&
  strides() const {
    return _strides;
  }

  //! The offset for indexing the bases.
  difference_type
  offset() const {
    return _offset;
  }

protected:

  //! Return the array index for the given index list.
  Index
  arrayIndex(const IndexList& indices) const {
    Index result = 0;
    for (size_type n = 0; n != Dimension; ++n) {
      result += _strides[n] * indices[n];
    }
    return result - _offset;
  }

  //@}
};
  
//----------------------------------------------------------------------------
//! \defgroup MultiArrayBaseEquality Equality Operators
//@{

//! Return true if the member data are equal.
/*! \relates MultiArrayBase */
template<std::size_t _Dimension>
inline
bool
operator==(const MultiArrayBase<_Dimension>& x, const MultiArrayBase<_Dimension>& y) {
  return x.extents() == y.extents() && x.bases() == y.bases() &&
    x.storage() == y.storage() && x.strides() == y.strides();
}

//! Return true if they are not equal.
/*! \relates MultiArrayBase */
template<std::size_t _Dimension>
inline
bool
operator!=(const MultiArrayBase<_Dimension>& x, const MultiArrayBase<_Dimension>& y) {
  return !(x == y);
}

//@}

END_NAMESPACE_ARRAY

#define __array_MultiArrayBase_ipp__
#include "MultiArrayBase.ipp"
#undef __array_MultiArrayBase_ipp__

#endif
