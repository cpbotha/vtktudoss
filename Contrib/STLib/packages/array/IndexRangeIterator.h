// -*- C++ -*-

/*!
  \file array/IndexRangeIterator.h
  \brief An index range iterator.
*/

#if !defined(__array_IndexRangeIterator_h__)
#define __array_IndexRangeIterator_h__

#include "IndexRange.h"

#include <iterator>

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_IndexRangeIterator)
#define DEBUG_array_IndexRangeIterator
#endif

BEGIN_NAMESPACE_ARRAY

//! An index range iterator.
template<std::size_t _Dimension>
class
IndexRangeIterator {
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

  //! An index range.
  typedef IndexRange<Dimension> Range;
  //! The size type.
  typedef typename Types::size_type size_type;
  //! An array index is a signed integer.
  typedef typename Types::Index Index;
  //! A list of indices.
  typedef typename Types::IndexList IndexList;
  //! A list of sizes.
  typedef typename Types::SizeList SizeList;

  // Iterator types.

  //! Random access iterator category.
  typedef std::random_access_iterator_tag iterator_category;
  //! Value type.
  typedef IndexList value_type;
  //! Pointer difference type.
  typedef typename Types::difference_type difference_type;
  //! Const reference to the value type.
  typedef const value_type& reference;
  //! Const pointer to the value type.
  typedef const value_type* pointer;

  //
  // Member data.
  //
private:

  //! An index list.
  IndexList _indexList;
  //! The rank of the index list.
  Index _rank;
  //! The index range.
  Range _range;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Return an iterator to the beginning of the index range.
  static
  IndexRangeIterator
  begin(const Range& range);

  //! Return an iterator to the end of the index range.
  static
  IndexRangeIterator
  end(const Range& range);

  //! Copy constructor.
  IndexRangeIterator(const IndexRangeIterator& other);

  //! Assignment operator.
  IndexRangeIterator&
  operator=(const IndexRangeIterator& other);

  //! Destructor.
  ~IndexRangeIterator() {
  }

private:

  //! Default constructor. Uninitialized memory.
  IndexRangeIterator() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! The index range.
  const Range&
  range() const {
    return _range;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Validity.
  //@{
private:

  //! Return true if the iterator is valid.
  /*!
    It's valid if it is in the range [begin(), end()).
  */
  bool
  isValid() const;

  //! Return true if the iterator is at the beginning.
  bool
  isBegin() const;

  //! Return true if the iterator is at the end.
  bool
  isEnd() const;

  //@}
  //--------------------------------------------------------------------------
  //! \name Forward iterator requirements.
  //@{
public:

  reference
  operator*() const {
    return _indexList;
  }

  pointer
  operator->() const {
    return &_indexList;
  }

  //! Pre-increment.
  IndexRangeIterator&
  operator++();

  //! Post-increment.
  /*!
    \warning This function is inefficient. Use pre-increment instead.
  */
  IndexRangeIterator
  operator++(int);

  //@}
  //--------------------------------------------------------------------------
  //! \name Bidirectional iterator requirements.
  //@{
public:

  //! Pre-decrement.
  IndexRangeIterator&
  operator--();

  //! Post-decrement.
  /*!
    \warning This function is inefficient. Use pre-increment instead.
  */
  IndexRangeIterator
  operator--(int);

  //@}
  //--------------------------------------------------------------------------
  //! \name Random access iterator requirements.
  //@{
public:

  //! Iterator indexing.
  /*!
    \warning This function is inefficient.
  */
  value_type
  operator[](const difference_type n) const {
    return *(*this + n);
  }

  IndexRangeIterator&
  operator+=(const difference_type n) {
    _rank += n;
    calculateIndexList();
    return *this;
  }

  IndexRangeIterator
  operator+(const difference_type n) const {
    IndexRangeIterator tmp(*this);
    tmp += n;
    return tmp;
  }

  IndexRangeIterator&
  operator-=(const difference_type n) {
    _rank -= n;
    calculateIndexList();
    return *this;
  }

  IndexRangeIterator
  operator-(const difference_type n) const {
    IndexRangeIterator tmp(*this);
    tmp -= n;
    return tmp;
  }

  Index
  base() const { 
    return _rank;
  }

private:

  //! Calculate the index list from the rank.
  void
  calculateIndexList();

  //@}
};
  
//---------------------------------------------------------------------------
// Equality.

//! Return true if the iterators are equal.
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator==(const IndexRangeIterator<_Dimension>& x,
	   const IndexRangeIterator<_Dimension>& y) {
#ifdef DEBUG_array_IndexRangeIterator
  // The must be iterators over the same index range.
  assert(x.range() == y.range());
#endif
  return x.base() == y.base();
}

//! Return true if they are not equal.
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator!=(const IndexRangeIterator<_Dimension>& x,
	   const IndexRangeIterator<_Dimension>& y) {
  return !(x == y);
}


//! Return true if the first precedes the second.
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator<(const IndexRangeIterator<_Dimension>& x,
	  const IndexRangeIterator<_Dimension>& y) {
#ifdef DEBUG_array_IndexRangeIterator
  // The must be iterators over the same index range.
  assert(x.range() == y.range());
#endif
  return x.base() < y.base();
}

//! Return y < x.
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator>(const IndexRangeIterator<_Dimension>& x,
	  const IndexRangeIterator<_Dimension>& y) {
  return y < x;
}

//! Return !(y < x).
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator<=(const IndexRangeIterator<_Dimension>& x,
	   const IndexRangeIterator<_Dimension>& y) {
  return ! (y < x);
}

//! Return !(x < y).
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
bool
operator>=(const IndexRangeIterator<_Dimension>& x,
	   const IndexRangeIterator<_Dimension>& y) {
  return ! (x < y);
}

//! Return the difference between the two iterators.
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
typename IndexRangeIterator<_Dimension>::difference_type
operator-(const IndexRangeIterator<_Dimension>& x,
	  const IndexRangeIterator<_Dimension>& y) {
  return x.base() - y.base();
}

//! Advance the iterator.
/*! \relates IndexRangeIterator */
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>
operator+(const typename IndexRangeIterator<_Dimension>::difference_type& n,
	  const IndexRangeIterator<_Dimension>& x) {
  return x + n;
}

END_NAMESPACE_ARRAY

#define __array_IndexRangeIterator_ipp__
#include "IndexRangeIterator.ipp"
#undef __array_IndexRangeIterator_ipp__

#endif
