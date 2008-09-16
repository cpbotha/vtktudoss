// -*- C++ -*-

/*!
  \file array/MultiArrayConstRef.h
  \brief Multi-dimensional constant array that references memory and has contiguous storage.
*/

#if !defined(__array_MultiArrayConstRef_h__)
#define __array_MultiArrayConstRef_h__

#include "MultiArrayConstView.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_MultiArrayConstRef)
#define DEBUG_array_MultiArrayConstRef
#endif

BEGIN_NAMESPACE_ARRAY

//! Multi-dimensional constant %array that references memory and has contiguous storage.
/*!
  <b>Constructors, etc.</b>

  Since this %array references externally allocated memory, there is
  no default constructor.

  You can construct an %array from a const pointer to the data and its 
  index extents.
  Below we make a 2x4x8 %array with index range [0..1]x[0..3]x[0..7]
  \verbatim
  double data[2 * 4 * 8];
  ...
  array::MultiArrayConstRef<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArrayConstRef<double, 3> a(data, extents); \endverbatim
  
  You can also specify the index bases. Below we make a
  2x4x8 %array with index range [-1..0]x[2..5]x[-3..4]
  \verbatim
  double data[2 * 4 * 8];
  ...
  array::MultiArrayConstRef<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArrayConstRef<double, 3>::IndexList bases(-1, 2, -3)
  array::MultiArrayConstRef<double, 3> a(data, extents, bases); \endverbatim
  
  The copy constructors create shallow copies of the argument, i.e. the 
  array data is referenced.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArrayConstRef<int, 3> b(a); \endverbatim
  The argument may be a MultiArray, MultiArrayRef, or a MultiArrayConstRef.
  The dimension and value type must be the same.

  Since this is a constant %array class, there are no assignment operators.

  You can use rebuild() to make a constant reference to another %array.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArrayConstRef<int, 3> b(a);
  array::MultiArray<int, 3> c(extents);
  b.rebuild(c); \endverbatim

  <b>Container Member Functions</b>
  
  MultiArrayConstRef inherits the following functionality for treating the 
  %array as a constant random access container.
  
  - MultiArrayBase::empty()
  - MultiArrayBase::size()
  - MultiArrayBase::max_size()

  It defines the following functions.

  - begin()
  - end()
  - rbegin()
  - rend()
  - operator[]()

  <b>%Array Indexing Member Functions</b>

  MultiArrayConstRef inherits the following %array indexing functionality.

  - MultiArrayBase::extents()
  - MultiArrayBase::bases()
  - MultiArrayBase::setBases()
  - MultiArrayBase::range()
  - MultiArrayBase::storage()
  - MultiArrayBase::strides()
  - MultiArrayBase::offset()
  - MultiArrayView::operator()()
  - MultiArrayView::view()

  <b>Free Functions</b>

  - \ref MultiArrayConstRefEquality
  - \ref MultiArrayConstRefFile
*/
template<typename _T, std::size_t _Dimension, typename _DataPointer = const _T*>
class
MultiArrayConstRef : 
  virtual public MultiArrayConstView<_T, _Dimension, _DataPointer> {
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

  typedef MultiArrayConstView<_T, _Dimension, _DataPointer> Base;
  typedef MultiArrayTypes<_T, _Dimension> Types;

public:

  // Types for STL compliance.

  //! The element type of the array.  
  typedef typename Types::value_type value_type;
  //! A pointer to a constant array element.
  typedef typename Types::const_pointer const_pointer;
  //! A iterator on constant elements in the array.
  typedef typename Types::const_iterator const_iterator;
  //! A reverse iterator on constant elements in the array.
  typedef typename Types::const_reverse_iterator const_reverse_iterator;
  //! A reference to a constant array element.
  typedef typename Types::const_reference const_reference;
  //! The size type.
  typedef typename Types::size_type size_type;
  //! Pointer difference type.
  typedef typename Types::difference_type difference_type;

  // Other types.

  //! The parameter type.
  /*! This is used for passing the value type as an argument. */
  typedef typename Types::Parameter Parameter;
  //! An array index is a signed integer.
  typedef typename Types::Index Index;
  //! A list of indices.
  typedef typename Types::IndexList IndexList;
  //! A list of sizes.
  typedef typename Types::SizeList SizeList;
  //! The storage order.
  typedef typename Types::Storage Storage;
  //! An index range.
  typedef typename Base::Range Range;
  //! A constant view of this array.
  typedef typename Base::ConstView ConstView;
  //! The data pointer type used to construct this array.
  typedef _DataPointer DataPointer;

  //
  // Using member data.
  //
protected:

  using Base::_data;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  // The default copy constructor is fine.

  //! Copy constructor for other data pointer types.
  template<typename _DP>
  MultiArrayConstRef(const MultiArrayConstRef<value_type, Dimension, _DP>& 
		     other);

  //! Construct from a pointer to the memory, the array extents, and optionally the storage order.
  MultiArrayConstRef(DataPointer data, const SizeList& extents, 
		     const Storage& storage = Storage(RowMajor()));

  //! Construct from a pointer to the memory, the array extents, the index bases, and optionally the storage order.
  MultiArrayConstRef(DataPointer data, const SizeList& extents,
		     const IndexList& bases,
		     const Storage& storage = Storage(RowMajor()));

  //! Destructor does not deallocate memory.
  virtual
  ~MultiArrayConstRef() {
  }

  //! Rebuild the data structure.
  /*! \note The size (number of elements) cannot change. */
  void
  rebuild(const SizeList& extents) {
    rebuild(extents, bases(), storage());
  }

  //! Rebuild the data structure.
  /*! \note The size (number of elements) cannot change. */
  void
  rebuild(const SizeList& extents, const IndexList& bases) {
    rebuild(extents, bases, storage());
  }

  //! Rebuild the data structure.
  /*! \note The size (number of elements) cannot change. */
  void
  rebuild(const SizeList& extents, const IndexList& bases,
	  const Storage& storage) {
    assert(product(extents) == size());
    rebuild(_data, extents, bases, storage);
  }

  //! Copy the data structure. Shallow copy of the elements.
  template<typename _DataPointer2>
  void
  rebuild(const MultiArrayConstRef<value_type, Dimension, _DataPointer2>& x) {
    Base::rebuild(x.data(), x.extents(), x.bases(), x.storage(), x.strides());
  }

protected:

  using Base::rebuild;

  //! Rebuild the data structure.
  void
  rebuild(DataPointer data, const SizeList& extents, const IndexList& bases,
	  const Storage& storage) {
    Base::rebuild(data, extents, bases, storage, 
		  computeStrides(extents, storage));
  }

  //! Compute the strides.
  /*!
    This is static so it can be called in the initializer list.
  */
  static
  IndexList
  computeStrides(const SizeList& extents, const Storage& storage);

private:

  //! Default constructor not implemented.
  MultiArrayConstRef() {
  }

  //! Assignment operator not implemented. You cannot assign to const data.
  MultiArrayConstRef&
  operator=(const MultiArrayConstRef& other);

  //@}
  //--------------------------------------------------------------------------
  //! \name Random Access Container.
  //@{
public:

  using Base::empty;
  using Base::size;
  using Base::max_size;

  //! Return a const iterator to the first value.
  const_iterator
  begin() const {
    return data();
  }

  //! Return a const iterator to one past the last value.
  const_iterator
  end() const {
    return data() + size();
  }

  //! Return a const reverse iterator to the end of the sequence.
  const_reverse_iterator
  rbegin() const {
    return const_reverse_iterator(end());
  }

  //! Return a const reverse iterator to the beginning of the sequence.
  const_reverse_iterator
  rend() const {
    return const_reverse_iterator(begin());
  }

  //! Container indexing.
  const_reference
  operator[](const size_type n) const {
    return data()[n];
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Array indexing.
  //@{
public:

  using Base::extents;
  using Base::bases;
  using Base::setBases;
  using Base::range;
  using Base::storage;
  using Base::strides;
  using Base::offset;
  using Base::data;
  using Base::view;

  //@}
};
  
//----------------------------------------------------------------------------
//! \defgroup MultiArrayConstRefEquality Equality and Comparison Operators
//@{

//! Return true if the arrays have the same extents and elements.
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer1,
	 typename _DataPointer2>
inline
bool
operator==(const MultiArrayConstRef<_T, _Dimension, _DataPointer1>& x,
	   const MultiArrayConstRef<_T, _Dimension, _DataPointer2>& y) {
   return x.extents() == y.extents() &&
     std::equal(x.begin(), x.end(), y.begin());
}

//! Return true if they are not equal.
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer1,
	 typename _DataPointer2>
inline
bool
operator!=(const MultiArrayConstRef<_T, _Dimension, _DataPointer1>& x,
	   const MultiArrayConstRef<_T, _Dimension, _DataPointer2>& y) {
  return !(x == y);
}


//! Lexicographical comparison of the elements.
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer1,
	 typename _DataPointer2>
inline
bool
operator<(const MultiArrayConstRef<_T, _Dimension, _DataPointer1>& x,
	  const MultiArrayConstRef<_T, _Dimension, _DataPointer2>& y) {
  return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
}

//! Return y < x.
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer1,
	 typename _DataPointer2>
inline
bool
operator>(const MultiArrayConstRef<_T, _Dimension, _DataPointer1>& x,
	  const MultiArrayConstRef<_T, _Dimension, _DataPointer2>& y) {
  return y < x;
}

//! Return !(y < x).
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer1,
	 typename _DataPointer2>
inline
bool
operator<=(const MultiArrayConstRef<_T, _Dimension, _DataPointer1>& x,
	   const MultiArrayConstRef<_T, _Dimension, _DataPointer2>& y) {
  return ! (y < x);
}

//! Return !(x < y).
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer1,
	 typename _DataPointer2>
inline
bool
operator>=(const MultiArrayConstRef<_T, _Dimension, _DataPointer1>& x,
	   const MultiArrayConstRef<_T, _Dimension, _DataPointer2>& y) {
  return ! (x < y);
}

//@}
//----------------------------------------------------------------------------
//! \defgroup MultiArrayConstRefFile File I/O
//@{

//! Print the index bases, array extents, and elements.
/*! \relates MultiArrayConstRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
std::ostream&
operator<<(std::ostream& out, 
	   const MultiArrayConstRef<_T, _Dimension, _DataPointer>& x) {
  out << x.extents() << '\n'
      << x.bases() << '\n'
      << x.storage() << '\n';
  std::copy(x.begin(), x.end(), std::ostream_iterator<_T>(out, "\n"));
  return out;
}

//@}

END_NAMESPACE_ARRAY

#define __array_MultiArrayConstRef_ipp__
#include "MultiArrayConstRef.ipp"
#undef __array_MultiArrayConstRef_ipp__

#endif
