// -*- C++ -*-

#if !defined(__array_MultiArrayRef_ipp__)
#error This file is an implementation detail of the class MultiArrayRef.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Copy constructor.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _DP>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>::
MultiArrayRef(const MultiArrayRef<value_type, Dimension, _DP>& other) :
  VirtualBase(other),
  // These do nothing.
  Base(other),
  ViewBase(other) {
}

// Copy constructor.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>::
MultiArrayRef(const MultiArrayRef& other) :
  VirtualBase(other),
  // These do nothing.
  Base(other),
  ViewBase(other) {
}

// Construct from a pointer to the memory, the array extents, and optionally
// the storage order.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>::
MultiArrayRef(DataPointer data, const SizeList& extents,
	      const Storage& storage) :
  VirtualBase(data, extents, IndexList(Index(0)), storage,
	      computeStrides(extents, storage)),
  // These do nothing.
  Base(data, extents, storage),
  ViewBase(data, extents, bases(), storage, strides()) {
}

// Construct from a pointer to the memory, the array extents, the index bases,
// and optionally the storage order.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>::
MultiArrayRef(DataPointer data, const SizeList& extents,
	      const IndexList& bases, const Storage& storage) :
  VirtualBase(data, extents, bases, storage, computeStrides(extents, storage)),
  // These do nothing.
  Base(data, extents, bases, storage),
  ViewBase(data, extents, bases, storage, strides()) {
}

// Assignment operator for other array views.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _T2, typename _DataPointer2>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
MultiArrayRef<_T, _Dimension, _DataPointer>::
operator=(const MultiArrayConstView<_T2, Dimension, _DataPointer2>& other) {
  typedef IndexRangeIterator<Dimension> Iterator;

#ifdef DEBUG_array_MultiArrayView
  // The arrays must have the same index range.
  assert(extents() == other.extents() && bases() == other.bases());
#endif

  // Copy the elements.
  if (storage() == other.storage()) {
    // If the storage orders are the same, we can use the built-in iterators.
    std::copy(other.begin(), other.end(), begin());
  }
  else {
    // If the storage orders differ, iterate over the index range and do 
    // array indexing.
    const Range range(extents(), bases());
    const Iterator end = Iterator::end(range);
    for (Iterator i = Iterator::begin(range); i != end; ++i) {
      (*this)(*i) = other(*i);
    }
  }
  return *this;
}

// Assignment operator for arrays with contiguous memory.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _T2, typename _DataPointer2>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
MultiArrayRef<_T, _Dimension, _DataPointer>::
operator=(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& other) {
  typedef IndexRangeIterator<Dimension> Iterator;

#ifdef DEBUG_array_MultiArrayView
  // The arrays must have the same index range.
  assert(extents() == other.extents() && bases() == other.bases());
#endif

  // Copy the elements.
  if (storage() == other.storage()) {
    // If the storage orders are the same, we can use the built-in iterators.
    std::copy(other.begin(), other.end(), begin());
  }
  else {
    // If the storage orders differ, iterate over the index range and do 
    // array indexing.
    const Range range(extents(), bases());
    const Iterator end = Iterator::end(range);
    for (Iterator i = Iterator::begin(range); i != end; ++i) {
      (*this)(*i) = other(*i);
    }
  }
  return *this;
}

// Assignment operator.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
MultiArrayRef<_T, _Dimension, _DataPointer>::
operator=(const MultiArrayRef& other) {
  typedef IndexRangeIterator<Dimension> Iterator;

  if (this != &other) {
#ifdef DEBUG_array_MultiArrayView
    // The arrays must have the same index range.
    assert(extents() == other.extents() && bases() == other.bases());
#endif

    // Copy the elements.
    if (storage() == other.storage()) {
      // If the storage orders are the same, we can use the built-in iterators.
      std::copy(other.begin(), other.end(), begin());
    }
    else {
      // If the storage orders differ, iterate over the index range and do 
      // array indexing.
      const Range range(extents(), bases());
      const Iterator end = Iterator::end(range);
      for (Iterator i = Iterator::begin(range); i != end; ++i) {
	(*this)(*i) = other(*i);
      }
    }
  }
  return *this;
}

//----------------------------------------------------------------------------
// Assignment operators with scalar operand.
//----------------------------------------------------------------------------

// Array-scalar addition.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator+=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i += value;
  }
  return x;
}

// Array-scalar subtraction.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator-=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i -= value;
  }
  return x;
}

// Array-scalar multiplication.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator*=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i *= value;
  }
  return x;
}

// Array-scalar division.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator/=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i /= value;
  }
  return x;
}

// Array-scalar modulus.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator%=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i %= value;
  }
  return x;
}

// Left shift.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator<<=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	    const int offset) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i <<= offset;
  }
  return x;
}

// Right shift.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayRef<_T, _Dimension, _DataPointer>&
operator>>=(MultiArrayRef<_T, _Dimension, _DataPointer>& x,
	    const int offset) {
  typedef typename MultiArrayRef<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i >>= offset;
  }
  return x;
}

END_NAMESPACE_ARRAY
