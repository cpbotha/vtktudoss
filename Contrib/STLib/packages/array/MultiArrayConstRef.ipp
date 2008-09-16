// -*- C++ -*-

#if !defined(__array_MultiArrayConstRef_ipp__)
#error This file is an implementation detail of the class MultiArrayConstRef.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Copy constructor.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _DP>
inline
MultiArrayConstRef<_T, _Dimension, _DataPointer>::
MultiArrayConstRef(const MultiArrayConstRef<value_type, Dimension, _DP>& other) :
  Base(other) {
}

// Construct from a pointer to the memory, the array extents, and optionally
// the storage order.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayConstRef<_T, _Dimension, _DataPointer>::
MultiArrayConstRef(DataPointer data, const SizeList& extents,
		   const Storage& storage) :
  Base(data, extents, IndexList(Index(0)), storage, 
       computeStrides(extents, storage)) {
}

// Construct from a pointer to the memory, the array extents, the index bases,
// and optionally the storage order.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayConstRef<_T, _Dimension, _DataPointer>::
MultiArrayConstRef(DataPointer data, const SizeList& extents,
		   const IndexList& bases, const Storage& storage) :
  Base(data, extents, bases, storage, computeStrides(extents, storage)) {
}

// Compute the strides.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
typename MultiArrayConstRef<_T, _Dimension, _DataPointer>::IndexList
MultiArrayConstRef<_T, _Dimension, _DataPointer>::
computeStrides(const SizeList& extents, const Storage& storage) {
  IndexList strides;
  Index s = 1;
  for (size_type i = 0; i != Dimension; ++i) {
    strides[storage[i]] = s;
    s *= extents[storage[i]];
  }
  return strides;
}


END_NAMESPACE_ARRAY
