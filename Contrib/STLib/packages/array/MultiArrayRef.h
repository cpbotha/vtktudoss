// -*- C++ -*-

/*!
  \file array/MultiArrayRef.h
  \brief Multi-dimensional array that references memory and has contiguous storage.
*/

#if !defined(__array_MultiArrayRef_h__)
#define __array_MultiArrayRef_h__

#include "MultiArrayConstRef.h"
#include "MultiArrayView.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_MultiArrayRef)
#define DEBUG_array_MultiArrayRef
#endif

BEGIN_NAMESPACE_ARRAY

//! Multi-dimensional %array that references memory and has contiguous storage.
/*!
  <b>Constructors, etc.</b>

  Since this %array references externally allocated memory, there is
  no default constructor.

  You can construct an %array from a pointer to the data and its index extents.
  Below we make a 2x4x8 %array with index range [0..1]x[0..3]x[0..7]
  \verbatim
  double data[2 * 4 * 8];
  array::MultiArrayRef<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArrayRef<double, 3> a(data, extents); \endverbatim
  
  You can also specify the index bases. Below we make a
  2x4x8 %array with index range [-1..0]x[2..5]x[-3..4]
  \verbatim
  double data[2 * 4 * 8];
  array::MultiArrayRef<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArrayRef<double, 3>::IndexList bases(-1, 2, -3)
  array::MultiArrayRef<double, 3> a(data, extents, bases); \endverbatim
  
  The copy constructors create shallow copies of the argument, i.e. the 
  array data is referenced.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArrayRef<int, 3> b(a); \endverbatim
  The argument may be a MultiArray, or a MultiArrayRef.
  The dimension and value type must be the same.

  The assignment operators copy the element values. The argument must have
  the same index ranges as the %array, though they can differ in the value
  type.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  {
    int* data = new int[product(extents)];
    array::MultiArrayRef<int, 3> b(data, extents);
    b = a;
  }
  {
    double* data = new double[product(extents)];
    array::MultiArray<double, 3> c(data, extents);
    c = a;
  } \endverbatim
  The argument may be any of the multidimensional %array types.

  You can use rebuild() to make a reference to another %array.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArrayRef<int, 3> b(a);
  array::MultiArray<int, 3> c(extents);
  b.rebuild(c); \endverbatim

  <b>Container Member Functions</b>
  
  MultiArrayRef inherits the following functionality for treating the %array as
  a constant random access container.
  
  - MultiArrayBase::empty()
  - MultiArrayBase::size()
  - MultiArrayBase::max_size()
  - MultiArrayConstRef::begin()
  - MultiArrayConstRef::end()
  - MultiArrayConstRef::rbegin()
  - MultiArrayConstRef::rend()

  It defines the following functions.

  - begin()
  - end()
  - rbegin()
  - rend()
  - operator[]()
  - fill()

  <b>%Array Indexing Member Functions</b>

  MultiArrayRef inherits the following %array indexing functionality.

  - MultiArrayBase::extents()
  - MultiArrayBase::bases()
  - MultiArrayBase::setBases()
  - MultiArrayBase::range()
  - MultiArrayBase::storage()
  - MultiArrayBase::strides()
  - MultiArrayBase::offset()
  - MultiArrayView::operator()()
  - MultiArrayView::view()

  It defines the following functions.

  - data()

  <b>Free Functions</b>

  - \ref MultiArrayRefAssignmentOperatorsScalar
  - \ref MultiArrayConstRefEquality
  - \ref MultiArrayConstRefFile
*/
template<typename _T, std::size_t _Dimension, typename _DataPointer = _T*>
class
MultiArrayRef : public MultiArrayConstRef<_T, _Dimension, _DataPointer>,
		public MultiArrayView<_T, _Dimension, _DataPointer> {
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

  typedef MultiArrayTypes<_T, _Dimension> Types;
  typedef MultiArrayConstRef<_T, _Dimension, _DataPointer> Base;
  typedef MultiArrayView<_T, _Dimension, _DataPointer> ViewBase;
  typedef MultiArrayConstView<_T, _Dimension, _DataPointer> VirtualBase;

public:

  // Types for STL compliance.

  //! The element type of the array.  
  typedef typename Types::value_type value_type;
  //! A pointer to an array element.
  typedef typename Types::pointer pointer;
  //! A pointer to a constant array element.
  typedef typename Types::const_pointer const_pointer;
  //! A iterator on elements in the array.
  typedef typename Types::iterator iterator;
  //! A iterator on constant elements in the array.
  typedef typename Types::const_iterator const_iterator;
  //! A reverse iterator on elements in the array.
  typedef typename Types::reverse_iterator reverse_iterator;
  //! A reverse iterator on constant elements in the array.
  typedef typename Types::const_reverse_iterator const_reverse_iterator;
  //! A reference to an array element.
  typedef typename Types::reference reference;
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
  typedef typename ViewBase::ConstView ConstView;
  //! A view of this array.
  typedef typename ViewBase::View View;

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

  //! Copy constructor for other data pointer types.
  template<typename _DP>
  MultiArrayRef(const MultiArrayRef<value_type, Dimension, _DP>& other);

  //! Copy constructor.
  MultiArrayRef(const MultiArrayRef& other);

  //! Construct from a pointer to the memory, the array extents, and optionally the storage order.
  MultiArrayRef(DataPointer data, const SizeList& extents, 
		const Storage& storage = Storage(RowMajor()));

  //! Construct from a pointer to the memory, the array extents, the index bases, and optionally the storage order.
  MultiArrayRef(DataPointer data, const SizeList& extents,
		const IndexList& bases,
		const Storage& storage = Storage(RowMajor()));

  //! Assignment operator for other array views.
  /*! \pre The arrays must have the same index range. */
  template<typename _T2, typename _DataPointer2>
  MultiArrayRef&
  operator=(const MultiArrayConstView<_T2, Dimension, _DataPointer2>& other);

  //! Assignment operator for arrays with contiguous memory.
  /*! 
    \pre The arrays must have the same index range.
    \note This version is faster than the assignment operator that takes a
    MultiArrayConstView as an argument because arrays with contiguous memory
    have faster iterators.
  */
  template<typename _T2, typename _DataPointer2>
  MultiArrayRef&
  operator=(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& other);

  //! Assignment operator.
  /*! \pre The arrays must have the same index range. */
  MultiArrayRef&
  operator=(const MultiArrayRef& other);

  //! Destructor does not deallocate memory.
  virtual
  ~MultiArrayRef() {
  }

  using Base::rebuild;

  //! Copy the data structure. Shallow copy of the elements.
  void
  rebuild(const MultiArrayRef& x) {
    rebuild(x.data(), x.extents(), x.bases(), x.storage(), x.strides());
  }

protected:

  using Base::computeStrides;

private:

  //! Default constructor not implemented.
  MultiArrayRef();

  //@}
  //--------------------------------------------------------------------------
  //! \name Random Access Container.
  //@{
public:

  using Base::empty;
  using Base::size;
  using Base::max_size;
  using Base::begin;
  using Base::end;
  using Base::rbegin;
  using Base::rend;
  using Base::operator[];

  //! Return an iterator to the first value.
  iterator
  begin() {
    return _data;
  }

  //! Return an iterator to one past the last value.
  iterator
  end() {
    return _data + size();
  }

  //! Return a reverse iterator to the end of the sequence.
  reverse_iterator
  rbegin() {
    return reverse_iterator(end());
  }

  //! Return a reverse iterator to the beginning of the sequence.
  reverse_iterator
  rend() {
    return reverse_iterator(begin());
  }

  //! Container indexing.
  reference
  operator[](const size_type n) {
    return _data[n];
  }

  //! Fill the array with the specified value.
  template<typename _T2>
  void
  fill(const _T2& value) {
    std::fill(begin(), end(), value);
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Array indexing accessors.
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
  using ViewBase::view;

  //! Return a pointer to the beginning of the data.
  pointer
  data() {
    return _data;
  }

  //@}
};
  
//----------------------------------------------------------------------------
//! \defgroup MultiArrayRefAssignmentOperatorsScalar Assignment Operators with Scalar Operand
//@{

//! Array-scalar addition.
/*! \relates MultiArrayRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator+=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value);

//! Array-scalar subtraction.
/*! \relates MultiArrayRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator-=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value);

//! Array-scalar multiplication.
/*! \relates MultiArrayRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator*=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value);

//! Array-scalar division.
/*!
  \relates MultiArrayRef 
  \note This does not check for division by zero as the value type may not be
  as number type.
*/
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator/=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value);

//! Array-scalar modulus.
/*!
  \relates MultiArrayRef 
  \note This does not check for division by zero as the value type may not be
  as number type.
*/
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator%=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value);

//! Left shift.
/*! \relates MultiArrayRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator<<=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	    const int offset);

//! Right shift.
/*! \relates MultiArrayRef */
template<typename _T, std::size_t _Dimension, typename _DataPointer>
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator>>=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	    const int offset);

//@}

//---------------------------------------------------------------------------
// File I/O.

// CONTINUE: Add input.

END_NAMESPACE_ARRAY

#define __array_MultiArrayRef_ipp__
#include "MultiArrayRef.ipp"
#undef __array_MultiArrayRef_ipp__

#endif
