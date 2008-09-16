// -*- C++ -*-

/*!
  \file array/IndexListIterator.h
  \brief An index list iterator.
*/

#if !defined(__array_IndexListIterator_h__)
#define __array_IndexListIterator_h__

#include "MultiArrayBase.h"

#include <iterator>

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_IndexListIterator)
#define DEBUG_array_IndexListIterator
#endif

BEGIN_NAMESPACE_ARRAY

//! An index list iterator.
/*!
  CONTINUE: I think I should only use IndexRangeIterator instead of this class.
*/
template<std::size_t _Dimension>
class
IndexListIterator {
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
  typedef MultiArrayBase<_Dimension> Array;
  //! The size type.
  typedef typename Types::size_type size_type;
  //! An array index is a signed integer.
  typedef typename Types::Index Index;
  //! A list of indices.
  typedef typename Types::IndexList IndexList;
  //! A list of sizes.
  typedef typename Types::SizeList SizeList;
  //! The storage order.
  typedef typename Types::Storage Storage;

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
  //! Pointer to the multi-array base.
  const Array* _array;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Return an iterator to the beginning of the index range.
  static
  IndexListIterator
  begin(const Array& array);

  //! Return an iterator to the end of the index range.
  static
  IndexListIterator
  end(const Array& array);

  //! Copy constructor.
  IndexListIterator(const IndexListIterator& other);

  //! Assignment operator.
  IndexListIterator&
  operator=(const IndexListIterator& other);

  //! Destructor.
  ~IndexListIterator() {
  }

private:

  //! Default constructor. Uninitialized memory.
  IndexListIterator() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! The index range.
  const Array*
  array() const {
    return _array;
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
  IndexListIterator&
  operator++();

  //! Post-increment.
  /*!
    \warning This function is inefficient. Use pre-increment instead.
  */
  IndexListIterator
  operator++(int);

  //@}
  //--------------------------------------------------------------------------
  //! \name Bidirectional iterator requirements.
  //@{
public:

  //! Pre-decrement.
  IndexListIterator&
  operator--();

  //! Post-decrement.
  /*!
    \warning This function is inefficient. Use pre-increment instead.
  */
  IndexListIterator
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

  IndexListIterator&
  operator+=(const difference_type n) {
    _rank += n;
    calculateIndexList();
    return *this;
  }

  IndexListIterator
  operator+(const difference_type n) const {
    IndexListIterator tmp(*this);
    tmp += n;
    return tmp;
  }

  IndexListIterator&
  operator-=(const difference_type n) {
    _rank -= n;
    calculateIndexList();
    return *this;
  }

  IndexListIterator
  operator-(const difference_type n) const {
    IndexListIterator tmp(*this);
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
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
bool
operator==(const IndexListIterator<_Dimension>& x,
	   const IndexListIterator<_Dimension>& y) {
#ifdef DEBUG_array_IndexListIterator
  // The must be iterators over the same index range.
  assert(x.array() == y.array());
#endif
  return x.base() == y.base();
}

//! Return true if they are not equal.
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
bool
operator!=(const IndexListIterator<_Dimension>& x,
	   const IndexListIterator<_Dimension>& y) {
  return !(x == y);
}


//! Return true if the first precedes the second.
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
bool
operator<(const IndexListIterator<_Dimension>& x,
	  const IndexListIterator<_Dimension>& y) {
#ifdef DEBUG_array_IndexListIterator
  // The must be iterators over the same index range.
  assert(x.array() == y.array());
#endif
  return x.base() < y.base();
}

//! Return y < x.
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
bool
operator>(const IndexListIterator<_Dimension>& x,
	  const IndexListIterator<_Dimension>& y) {
  return y < x;
}

//! Return !(y < x).
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
bool
operator<=(const IndexListIterator<_Dimension>& x,
	   const IndexListIterator<_Dimension>& y) {
  return ! (y < x);
}

//! Return !(x < y).
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
bool
operator>=(const IndexListIterator<_Dimension>& x,
	   const IndexListIterator<_Dimension>& y) {
  return ! (x < y);
}

//! Return the difference between the two iterators.
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
typename IndexListIterator<_Dimension>::difference_type
operator-(const IndexListIterator<_Dimension>& x,
	  const IndexListIterator<_Dimension>& y) {
  return x.base() - y.base();
}

//! Advance the iterator.
/*! \relates IndexListIterator */
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>
operator+(const typename IndexListIterator<_Dimension>::difference_type& n,
	  const IndexListIterator<_Dimension>& x) {
  return x + n;
}

END_NAMESPACE_ARRAY

#define __array_IndexListIterator_ipp__
#include "IndexListIterator.ipp"
#undef __array_IndexListIterator_ipp__

#endif
