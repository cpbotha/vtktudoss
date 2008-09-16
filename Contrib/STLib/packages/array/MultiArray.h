// -*- C++ -*-

/*!
  \file array/MultiArray.h
  \brief Multi-dimensional array that allocates its memory and has contiguous storage.
*/

#if !defined(__array_MultiArray_h__)
#define __array_MultiArray_h__

#include "MultiArrayRef.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_MultiArray)
#define DEBUG_array_MultiArray
#endif

BEGIN_NAMESPACE_ARRAY

//! Multi-dimensional %array that allocates it memory and has contiguous storage.
/*!
  <b>Constructors, etc.</b>
 
  The default constructor creates an empty %array. Below we make a 3-D %array
  of double precision floating point numbers.
  \verbatim
  array::MultiArray<double, 3> a; \endverbatim

  You can construct an %array from its index extents. Below we make a
  2x4x8 %array with index range [0..1]x[0..3]x[0..7]
  \verbatim
  array::MultiArray<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArray<double, 3> a(extents); \endverbatim
  
  You can also specify the index bases. Below we make a
  2x4x8 %array with index range [-1..0]x[2..5]x[-3..4]
  \verbatim
  array::MultiArray<double, 3>::SizeList extents(2, 4, 8)
  array::MultiArray<double, 3>::IndexList bases(-1, 2, -3)
  array::MultiArray<double, 3> a(extents, bases); \endverbatim
  
  The copy constructors create (deep) copies of the argument.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArray<int, 3> b(a);
  array::MultiArray<double, 3> c = a; \endverbatim
  The argument may be a MultiArray, a MultiArrayRef, or a MultiArrayConstRef .
  The dimension must be the same, but the value type may differ.

  The assignment operators copy the element values. The argument must have
  the same index ranges as the %array, though they can differ in the value
  type.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArray<int, 3> b(extents);
  b = a;
  array::MultiArray<double, 3> c(extents);
  c = a; \endverbatim
  The argument may be any of the multidimensional %array types.

  You can change the shape of an %array with rebuild(). You can specify 
  the extents or both the extents and the bases. If the new shape has different
  number of elements than the old, new memory will be allocated.
  \verbatim
  array::MultiArray<int, 3> a;
  a.rebuild(extents);
  a.rebuild(extents, bases); \endverbatim
  
  You can also use rebuild() to make a copy of another %array. Again, this
  will allocate new memory if necessary.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArray<int, 3> b;
  b.rebuild(a);
  array::MultiArray<double, 3> c;
  c.rebuild(a); \endverbatim

  <b>Container Member Functions</b>
  
  MultiArray inherits the following functionality for treating the %array as
  as a random access container.
  
  - MultiArrayBase::empty()
  - MultiArrayBase::size()
  - MultiArrayBase::max_size()
  - MultiArrayRef::begin()
  - MultiArrayRef::end()
  - MultiArrayRef::rbegin()
  - MultiArrayRef::rend()
  - MultiArrayRef::operator[]()
  - MultiArrayRef::fill()

  <b>%Array Indexing Member Functions</b>

  MultiArray inherits the following %array indexing functionality.

  - MultiArrayBase::extents()
  - MultiArrayBase::bases()
  - MultiArrayBase::setBases()
  - MultiArrayBase::range()
  - MultiArrayBase::storage()
  - MultiArrayBase::strides()
  - MultiArrayBase::offset()
  - MultiArrayView::operator()()
  - MultiArrayView::view()
  - MultiArrayRef::data()

  <b>Free Functions</b>

  - \ref MultiArrayRefAssignmentOperatorsScalar
  - \ref MultiArrayConstRefEquality
  - \ref MultiArrayConstRefFile
*/
template<typename _T, std::size_t _Dimension>
class
MultiArray : public MultiArrayRef<_T, _Dimension> {
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
  typedef MultiArrayRef<_T, _Dimension> Base;
  typedef MultiArrayConstView<_T, _Dimension, _T*> VirtualBase;

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
  typedef typename Base::ConstView ConstView;
  //! A view of this array.
  typedef typename Base::View View;

  //
  // Using member data.
  //
protected:

  using Base::_data;
  
  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Default constructor. Empty array.
  MultiArray();

  //! Copy constructor for different types.
  template<typename _T2, typename _DataPointer2>
  MultiArray(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& other);

  //! Copy constructor.
  MultiArray(const MultiArray& other);

  //! Construct from the array extents, and optionally the storage order.
  MultiArray(const SizeList& extents,
	     const Storage& storage = Storage(RowMajor()));

  //! Construct from the array extents, the index bases, and optionally the storage order.
  MultiArray(const SizeList& extents, const IndexList& bases,
	     const Storage& storage = Storage(RowMajor()));

  //! Assignment operator for other array views.
  /*! \pre The arrays must have the same index range. */
  template<typename _T2, typename _DataPointer2>
  MultiArray&
  operator=(const MultiArrayConstView<_T2, Dimension, _DataPointer2>& other);

  //! Assignment operator for arrays with contiguous memory.
  /*! 
    \pre The arrays must have the same index range.
    \note This version is faster than the assignment operator that takes a
    MultiArrayConstView as an argument because arrays with contiguous memory
    have faster iterators.
  */
  template<typename _T2, typename _DataPointer2>
  MultiArray&
  operator=(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& other);

  //! Assignment operator.
  /*! \pre The arrays must have the same index range. */
  MultiArray&
  operator=(const MultiArray& other);

  //! Destructor. Deallocate the memory.
  virtual
  ~MultiArray() {
    delete[] _data;
    _data = 0;
  }

  //! Rebuild the data structure. Re-allocate memory if the size changes.
  void
  rebuild(const SizeList& extents) {
    rebuild(extents, bases(), storage());
  }

  //! Rebuild the data structure. Re-allocate memory if the size changes.
  void
  rebuild(const SizeList& extents, const IndexList& bases) {
    rebuild(extents, bases, storage());
  }

  //! Rebuild the data structure. Re-allocate memory if the size changes.
  void
  rebuild(const SizeList& extents, const IndexList& bases,
	  const Storage& storage);

  //! Copy the data structure. Deepy copy of the elements.
  template<typename _T2, typename _DataPointer2>
  void
  rebuild(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& x) {
    // Set the array shape and allocate memory if necessary.
    rebuild(x.extents(), x.bases(), x.storage());
    // Copy the elements.
    std::copy(x.begin(), x.end(), begin());
  }

protected:

  using Base::computeStrides;

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
  using Base::fill;

  //@}
  //--------------------------------------------------------------------------
  //! \name Array Indexing.
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
  
//---------------------------------------------------------------------------
// File I/O.

// CONTINUE: Add input.

END_NAMESPACE_ARRAY

#define __array_MultiArray_ipp__
#include "MultiArray.ipp"
#undef __array_MultiArray_ipp__

#endif
