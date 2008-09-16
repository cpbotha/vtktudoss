// -*- C++ -*-

/*!
  \file array/IndexRange.h
  \brief An index range.
*/

#if !defined(__array_IndexRange_h__)
#define __array_IndexRange_h__

#include "IndexTypes.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_IndexRange)
#define DEBUG_array_IndexRange
#endif

BEGIN_NAMESPACE_ARRAY

//! An index range.
template<std::size_t _Dimension>
class
IndexRange {
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

  // Other types.

  //! An array index is a signed integer.
  typedef typename Types::Index Index;
  //! A list of indices.
  typedef typename Types::IndexList IndexList;
  //! A list of sizes.
  typedef typename Types::SizeList SizeList;

  //
  // Member data.
  //
protected:

  //! The array extents.
  SizeList _extents;
  //! The lower bound for each index.
  IndexList _bases;
  //! The steps.
  IndexList _steps;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  // The default copy constructor and assignment operator are fine.

  //! Default constructor. Unitilialized data.
  IndexRange() {
  }

  //! Construct from the extents, and optionally the bases and the steps.
  IndexRange(const SizeList& extents, 
	     const IndexList& bases = IndexList(Index(0)),
	     const IndexList& steps = IndexList(Index(1))) :
    _extents(extents),
    _bases(bases),
    _steps(steps) {
  }

#if 0
  //! Construct from the lower and upper bounds, and optionally the steps.
  IndexRange(const IndexList& lower, const IndexList& upper,
	     const IndexList& steps = IndexList(Index(1)));
#endif

  //! Destructor does nothing.
  ~IndexRange() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! The index range extents.
  const SizeList&
  extents() const {
    return _extents;
  }

  //! The bases.
  const IndexList&
  bases() const {
    return _bases;
  }

#if 0
  //! The (open) upper bound for each index.
  IndexList
  upper() const {
    return _bases + _steps * _extents;
  }
#endif

  //! The steps.
  const IndexList&
  steps() const {
    return _steps;
  }

  //@}
};
  
  
//---------------------------------------------------------------------------
// Free functions.

//! Return the intersection of the two index ranges.
/*! \relates IndexRange */
template<std::size_t _Dimension, typename _Size, typename _Index>
IndexRange<_Dimension>
overlap(const Array<_Size, _Dimension>& extents1,
	const Array<_Index, _Dimension>& bases1,
	const Array<_Size, _Dimension>& extents2,
	const Array<_Index, _Dimension>& bases2);

//! Return the intersection of the two ranges.
/*!
  \pre The ranges must have unit steps.
  \relates IndexRange
*/
template<std::size_t _Dimension>
IndexRange<_Dimension>
overlap(const IndexRange<_Dimension>& x, const IndexRange<_Dimension>& y);

//! Return true if the index is in the index range.
/*! \relates IndexRange */
template<std::size_t _Dimension>
bool
isIn(const IndexRange<_Dimension>& range, 
     const typename IndexRange<_Dimension>::IndexList& index);

// CONTINUE: REMOVE
#if 0
//! Return the intersection of the two index ranges.
/*! \relates IndexRange */
template<std::size_t _Dimension>
inline
IndexRange<_Dimension>
overlap(const MultiArrayBase<_Dimension>& x,
	const MultiArrayBase<_Dimension>& y) {
  return overlap(x.extents(), x.bases(), y.extents(), y.bases());
}
#endif

//---------------------------------------------------------------------------
// Equality.

//! Return true if the member data are equal.
/*! \relates IndexRange */
template<std::size_t _Dimension>
inline
bool
operator==(const IndexRange<_Dimension>& x, const IndexRange<_Dimension>& y) {
  return x.extents() == y.extents() && x.bases() == y.bases() &&
    x.steps() == y.steps();
}

//! Return true if they are not equal.
/*! \relates IndexRange */
template<std::size_t _Dimension>
inline
bool
operator!=(const IndexRange<_Dimension>& x, const IndexRange<_Dimension>& y) {
  return !(x == y);
}

//---------------------------------------------------------------------------
// File I/O.

//! Print the extents, bases, and steps.
/*! \relates IndexRange */
template<std::size_t _Dimension>
inline
std::ostream&
operator<<(std::ostream& out, const IndexRange<_Dimension>& x) {
  out << x.extents() << '\n'
      << x.bases() << '\n'
      << x.steps() << '\n';
  return out;
}

//! Read the extents, bases, and steps.
/*! \relates IndexRange */
template<std::size_t _Dimension>
inline
std::istream&
operator>>(std::istream& in, IndexRange<_Dimension>& x) {
  typedef IndexRange<_Dimension> IndexRange;
  typename IndexRange::SizeList extents;
  typename IndexRange::IndexList bases, steps;
  in >> extents >> bases >> steps;
  x = IndexRange(extents, bases, steps);
  return in;
}

END_NAMESPACE_ARRAY

#define __array_IndexRange_ipp__
#include "IndexRange.ipp"
#undef __array_IndexRange_ipp__

#endif
