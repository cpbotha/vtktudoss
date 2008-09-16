// -*- C++ -*-

#if !defined(__array_MultiArrayView_ipp__)
#error This file is an implementation detail of the class MultiArrayView.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Copy constructor.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _DP>
inline
MultiArrayView<_T, _Dimension, _DataPointer>::
MultiArrayView(const MultiArrayView<value_type, Dimension, _DP>& other) :
  Base(other) {
}

// Construct from a pointer to the memory, the array extents, the index bases,
// the storage order, and the strides.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayView<_T, _Dimension, _DataPointer>::
MultiArrayView(DataPointer data, const SizeList& extents, 
	       const IndexList& bases, const Storage& storage,
	       const IndexList& strides) :
  Base(data, extents, bases, storage, strides) {
}

// Assignment operator for other array views.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _T2, typename _DataPointer2>
inline
MultiArrayView<_T, _Dimension, _DataPointer>&
MultiArrayView<_T, _Dimension, _DataPointer>::
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
MultiArrayView<_T, _Dimension, _DataPointer>&
MultiArrayView<_T, _Dimension, _DataPointer>::
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
MultiArrayView<_T, _Dimension, _DataPointer>&
MultiArrayView<_T, _Dimension, _DataPointer>::
operator=(const MultiArrayView& other) {
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

//--------------------------------------------------------------------------
// Manipulators.

#if 0
// For the overlapping elements, set the first array to the unary function of 
// the second.
template<typename _T1, std::size_t _Dimension, typename _DataPointer1,
	 typename _T2, typename _DataPointer2,
	 typename _UnaryFunction>
inline
void
applyUnaryToOverlap
(MultiArrayConstView<_T1, _Dimension, _DataPointer1>* a,
 const MultiArrayConstView<_T2, _Dimension, _DataPointer2>& b,
 _UnaryFunction f) {
  typedef IndexRange<_Dimension> Range;
  typedef IndexRangeIterator<_Dimension> Iterator;

  const Range range = intersect(Range(a->extents(), a->bases()),
				Range(b.extents(), b.bases()));
  const Iterator end = Iterator::end(range);
  for (Iterator i = Iterator::begin(range); i != end; ++i) {
    (*a)(*i) = f(b(*i));
  }
}
#endif

//----------------------------------------------------------------------------
// Assignment operators with scalar operand.

// Array-scalar addition.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayView<_T, _Dimension, _DataPointer>&
operator+=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
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
MultiArrayView<_T, _Dimension, _DataPointer>&
operator-=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
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
MultiArrayView<_T, _Dimension, _DataPointer>&
operator*=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
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
MultiArrayView<_T, _Dimension, _DataPointer>&
operator/=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
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
MultiArrayView<_T, _Dimension, _DataPointer>&
operator%=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	   typename boost::call_traits<_T>::param_type value) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
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
MultiArrayView<_T, _Dimension, _DataPointer>&
operator<<=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	    const int offset) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
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
MultiArrayView<_T, _Dimension, _DataPointer>&
operator>>=(MultiArrayView<_T, _Dimension, _DataPointer>& x,
	    const int offset) {
  typedef typename MultiArrayView<_T, _Dimension, _DataPointer>::iterator
    iterator;
  const iterator end = x.end();
  for (iterator i = x.begin(); i != end; ++i) {
    *i >>= offset;
  }
  return x;
}

END_NAMESPACE_ARRAY
