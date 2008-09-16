// -*- C++ -*-

#if !defined(__array_MultiArrayConstView_ipp__)
#error This file is an implementation detail of the class MultiArrayConstView.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Copy constructor.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
template<typename _DP>
inline
MultiArrayConstView<_T, _Dimension, _DataPointer>::
MultiArrayConstView(const MultiArrayConstView<value_type, Dimension, _DP>& other) :
  Base(other),
  _data(other.data()) {
}

// Construct from a pointer to the memory, the array extents, the index bases,
// the storage order, and the strides.
template<typename _T, std::size_t _Dimension, typename _DataPointer>
inline
MultiArrayConstView<_T, _Dimension, _DataPointer>::
MultiArrayConstView(DataPointer data, const SizeList& extents, 
		    const IndexList& bases, const Storage& storage,
		    const IndexList& strides) :
  Base(extents, bases, storage, strides),
  _data(data) {
}

END_NAMESPACE_ARRAY
