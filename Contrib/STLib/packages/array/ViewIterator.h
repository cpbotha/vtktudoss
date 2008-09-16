// -*- C++ -*-

/*!
  \file array/ViewIterator.h
  \brief An iterator for a view of an array.
*/

#if !defined(__array_ViewIterator_h__)
#define __array_ViewIterator_h__

#include "IndexTypes.h"

#include "../third-party/loki/TypeManip.h"

#include <iterator>

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_ViewIterator)
#define DEBUG_array_ViewIterator
#endif

BEGIN_NAMESPACE_ARRAY

//! An iterator for a view of an array.
template<typename _MultiArray, bool _IsConst>
class
ViewIterator {
  //
  // Enumerations.
  //
public:

  //! The number of dimensions.
  enum {Dimension = _MultiArray::Dimension};

  //
  // Types.
  //
private:

  typedef IndexTypes<Dimension> Types;

public:

  //! The multi-array type.
  typedef _MultiArray MultiArray;
  //! Reference to the multi-array.
  typedef typename Loki::Select<_IsConst, const MultiArray&, 
				MultiArray&>::Result MultiArrayReference;
  //! Pointer to the multi-array.
  typedef typename Loki::Select<_IsConst, const MultiArray*, 
				MultiArray*>::Result MultiArrayPointer;
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
  typedef typename MultiArray::value_type value_type;
  //! Pointer difference type.
  typedef typename Types::difference_type difference_type;
  //! Reference to the value type.
  typedef typename Loki::Select<_IsConst, const value_type&, value_type&>::
  Result reference;
  //! Pointer to the value type.
  typedef typename Loki::Select<_IsConst, const value_type*, value_type*>::
  Result pointer;

  //
  // Member data.
  //
private:

  //! An index list.
  IndexList _indexList;
  //! The rank of the index list.
  Index _rank;
  //! Pointer in the array.
  pointer _iterator;
  //! Pointer to the multi-array..
  MultiArrayPointer _array;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Return an iterator to the beginning of the index range.
  static
  ViewIterator
  begin(MultiArrayReference array);

  //! Return an iterator to the end of the index range.
  static
  ViewIterator
  end(MultiArrayReference array);

  // The default copy constructor, assignment operator and destructor are fine.

  //! Copy constructor from non-const.
  template<bool _IsConst2>
  ViewIterator(const ViewIterator<MultiArray, _IsConst2>& other);

  //! Assignment operator from non-const.
  template<bool _IsConst2>
  ViewIterator&
  operator=(const ViewIterator<MultiArray, _IsConst2>& other);

private:

  //! Default constructor. Uninitialized data.
  ViewIterator() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! The index list.
  IndexList
  indexList() const {
    return _indexList;
  }

  //! The multi-array.
  MultiArrayPointer
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
    return *_iterator;
  }

  pointer
  operator->() const {
    return _iterator;
  }

  //! Pre-increment.
  ViewIterator&
  operator++();

  //! Post-increment.
  /*!
    \warning This function is inefficient. Use pre-increment instead.
  */
  ViewIterator
  operator++(int);

  //@}
  //--------------------------------------------------------------------------
  //! \name Bidirectional iterator requirements.
  //@{
public:

  //! Pre-decrement.
  ViewIterator&
  operator--();

  //! Post-decrement.
  /*!
    \warning This function is inefficient. Use pre-increment instead.
  */
  ViewIterator
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

  ViewIterator&
  operator+=(const difference_type n) {
    _rank += n;
    update();
    return *this;
  }

  ViewIterator
  operator+(const difference_type n) const {
    ViewIterator tmp(*this);
    tmp += n;
    return tmp;
  }

  ViewIterator&
  operator-=(const difference_type n) {
    _rank -= n;
    update();
    return *this;
  }

  ViewIterator
  operator-(const difference_type n) const {
    ViewIterator tmp(*this);
    tmp -= n;
    return tmp;
  }

  Index
  rank() const { 
    return _rank;
  }

  pointer
  base() const { 
    return _iterator;
  }

private:

  //! Calculate the index list and the base iterator from the rank.
  void
  update();

  //@}
};
  
//---------------------------------------------------------------------------
// Equality.

//! Return true if the iterators are equal.
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
bool
operator==(const ViewIterator<_MultiArray, _IsConst1>& x,
	   const ViewIterator<_MultiArray, _IsConst2>& y) {
#ifdef DEBUG_array_ViewIterator
  // The must be iterators over the same index range.
  assert(x.array() == y.array());
#endif
  return x.rank() == y.rank();
}

//! Return true if they are not equal.
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
bool
operator!=(const ViewIterator<_MultiArray, _IsConst1>& x,
	   const ViewIterator<_MultiArray, _IsConst2>& y) {
  return !(x == y);
}


//! Return true if the first precedes the second.
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
bool
operator<(const ViewIterator<_MultiArray, _IsConst1>& x,
	  const ViewIterator<_MultiArray, _IsConst2>& y) {
#ifdef DEBUG_array_ViewIterator
  // The must be iterators over the same index range.
  assert(x.array() == y.array());
#endif
  return x.rank() < y.rank();
}

//! Return y < x.
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
bool
operator>(const ViewIterator<_MultiArray, _IsConst1>& x,
	  const ViewIterator<_MultiArray, _IsConst2>& y) {
  return y < x;
}

//! Return !(y < x).
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
bool
operator<=(const ViewIterator<_MultiArray, _IsConst1>& x,
	   const ViewIterator<_MultiArray, _IsConst2>& y) {
  return ! (y < x);
}

//! Return !(x < y).
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
bool
operator>=(const ViewIterator<_MultiArray, _IsConst1>& x,
	   const ViewIterator<_MultiArray, _IsConst2>& y) {
  return ! (x < y);
}

//! Return the difference between the two iterators.
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst1, bool _IsConst2>
inline
typename ViewIterator<_MultiArray, _IsConst1>::difference_type
operator-(const ViewIterator<_MultiArray, _IsConst1>& x,
	  const ViewIterator<_MultiArray, _IsConst2>& y) {
  return x.rank() - y.rank();
}

//! Advance the iterator.
/*! \relates ViewIterator */
template<typename _MultiArray, bool _IsConst>
inline
ViewIterator<_MultiArray, _IsConst>
operator+(const typename ViewIterator<_MultiArray, _IsConst>::
	  difference_type& n,
	  const ViewIterator<_MultiArray, _IsConst>& x) {
  return x + n;
}

END_NAMESPACE_ARRAY

#define __array_ViewIterator_ipp__
#include "ViewIterator.ipp"
#undef __array_ViewIterator_ipp__

#endif
