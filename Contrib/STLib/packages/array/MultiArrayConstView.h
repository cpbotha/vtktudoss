// -*- C++ -*-

/*!
  \file array/MultiArrayConstView.h
  \brief Multi-dimensional constant view of an array.
*/

#if !defined(__array_MultiArrayConstView_h__)
#define __array_MultiArrayConstView_h__

#include "MultiArrayTypes.h"
#include "MultiArrayBase.h"
#include "ViewIterator.h"

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_MultiArrayConstView)
#define DEBUG_array_MultiArrayConstView
#endif

BEGIN_NAMESPACE_ARRAY

//! Multi-dimensional constant view of an %array.
/*!
  <b>Constructors, etc.</b>

  Since this %array references externally allocated memory, there is
  no default constructor. This class uses the automatically-generated 
  copy constructor; the array data is referenced. You can create
  an instance of this class with the view() member function.

  The copy constructors create shallow copies of the argument, i.e. the 
  array data is referenced.
  \verbatim
  array::MultiArray<int, 3> a(extents);
  array::MultiArrayConstView<int, 3> b(a); \endverbatim
  The argument may be any multidimensional %array type, however
  the dimension and value type must be the same.

  <b>Container Member Functions</b>
  
  MultiArrayConstView inherits the following functionality for treating 
  the %array as a random access container.
  
  - MultiArrayBase::empty()
  - MultiArrayBase::size()
  - MultiArrayBase::max_size()

  It defines the following functions.

  - begin()
  - end()
  - rbegin()
  - rend()

  <b>%Array Indexing Member Functions</b>

  MultiArrayConstView inherits the following %array indexing functionality.

  - MultiArrayBase::extents()
  - MultiArrayBase::bases()
  - MultiArrayBase::setBases()
  - MultiArrayBase::range()
  - MultiArrayBase::storage()
  - MultiArrayBase::strides()
  - MultiArrayBase::offset()

  It defines the following functions.

  - operator()()
  - data()
  - view()

  <b>Free Functions</b>

*/
template<typename _T, std::size_t _Dimension, typename _DataPointer = const _T*>
class
MultiArrayConstView : public MultiArrayBase<_Dimension> {
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

  typedef MultiArrayBase<_Dimension> Base;
  typedef MultiArrayTypes<_T, _Dimension> Types;

public:

  // Types for STL compliance.

  //! The element type of the array.  
  typedef typename Types::value_type value_type;
  //! A pointer to a constant array element.
  typedef typename Types::const_pointer const_pointer;
  //! An iterator on constant elements in the array.
  typedef ViewIterator<MultiArrayConstView, true> const_iterator;
  //! A reverse iterator on constant elements in the array.
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
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
  typedef MultiArrayConstView<_T, _Dimension, const _T*> ConstView;
  //! The data pointer type used to construct this array.
  typedef _DataPointer DataPointer;

  //
  // Member data.
  //
protected:

  //! Pointer to the beginning of a contiguous block of data.
  DataPointer _data;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  // The default copy constructor is fine.

  //! Copy constructor for other data pointer types.
  template<typename _DP>
  MultiArrayConstView(const MultiArrayConstView<value_type, Dimension, _DP>& other);

  //! Construct from a pointer to the memory, the array extents, the index bases, the storage order, and the strides.
  MultiArrayConstView(DataPointer data, const SizeList& extents, 
		      const IndexList& bases, const Storage& storage,
		      const IndexList& strides);

  //! Destructor does not deallocate memory.
  virtual
  ~MultiArrayConstView() {
  }

protected:

  //! Rebuild the data structure.
  void
  rebuild(DataPointer data, const SizeList& extents, const IndexList& bases,
	  const Storage& storage, const IndexList& strides) {
    _data = data;
    Base::rebuild(extents, bases, storage, strides);
  }

  //! Copy the data structure. Shallow copy of the elements.
  template<typename _DataPointer2>
  void
  rebuild(const MultiArrayConstView<value_type, Dimension, _DataPointer2>& x) {
    rebuild(x.data(), x.extents(), x.bases(), x.storage(), x.strides());
  }

private:

  //! Default constructor not implemented.
  /*!
    This class is a virtual base for other classes. Making the default 
    constructor private makes sure this class is appropriately constructed.
  */
  MultiArrayConstView();

  //! Assignment operator not implemented. You cannot assign to const data.
  MultiArrayConstView&
  operator=(const MultiArrayConstView& other);

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
    return const_iterator::begin(*this);
  }

  //! Return a const iterator to one past the last value.
  const_iterator
  end() const {
    return const_iterator::end(*this);
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

  //! Array indexing.
  const_reference
  operator()(const IndexList& indices) const {
    return _data[arrayIndex(indices)];
  }

  //! Return a const pointer to the beginning of the data.
  const_pointer
  data() const {
    return _data;
  }

  //! Make a sub-array view with the specified index range.
  /*! The bases for the view are the same as that for the index range. */
  ConstView
  view(const Range& range) const {
    return ConstView(&(*this)(range.bases()),
		     range.extents(), range.bases(), storage(),
		     strides() * range.steps());
  }

protected:

  using Base::arrayIndex;

  //@}
};

// CONTINUE: Add equality and file output.

END_NAMESPACE_ARRAY

#define __array_MultiArrayConstView_ipp__
#include "MultiArrayConstView.ipp"
#undef __array_MultiArrayConstView_ipp__

#endif
