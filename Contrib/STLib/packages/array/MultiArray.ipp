// -*- C++ -*-

#if !defined(__array_MultiArray_ipp__)
#error This file is an implementation detail of the class MultiArray.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Default constructor.
template<typename _T, std::size_t _Dimension>
inline
MultiArray<_T, _Dimension>::
MultiArray() :
  // Null pointer to memory and zero extents.
  VirtualBase(0, SizeList(size_type(0)), IndexList(Index(0)), RowMajor(),
	      IndexList(Index(1))),
  // This does nothing.
  Base(0, extents()) {
}

// Copy constructor for different types.
template<typename _T, std::size_t _Dimension>
template<typename _T2, typename _DataPointer2>
inline
MultiArray<_T, _Dimension>::
MultiArray(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& other) :
  VirtualBase(new value_type[other.size()], other.extents(), other.bases(),
	      other.storage(), other.strides()),
  // This does nothing.
  Base(_data, extents(), bases(), storage()) {
  // Copy the elements.
  std::copy(other.begin(), other.end(), begin());
}

// Copy constructor.
template<typename _T, std::size_t _Dimension>
inline
MultiArray<_T, _Dimension>::
MultiArray(const MultiArray& other) :
  VirtualBase(new value_type[other.size()], other.extents(), other.bases(),
	      other.storage(), other.strides()),
  // This does nothing.
  Base(_data, extents(), bases(), storage()) {
  // Copy the elements.
  std::copy(other.begin(), other.end(), begin());
}

// Construct from the array extents, and optionally the storage order.
template<typename _T, std::size_t _Dimension>
inline
MultiArray<_T, _Dimension>::
MultiArray(const SizeList& extents, const Storage& storage) :
  VirtualBase(new value_type[product(extents)], extents, IndexList(Index(0)),
	      storage, computeStrides(extents, storage)),
  // This does nothing.
  Base(_data, extents, storage) {
}


// Construct from the array extents, the index bases, and optionally the storage order.
template<typename _T, std::size_t _Dimension>
inline
MultiArray<_T, _Dimension>::
MultiArray(const SizeList& extents, const IndexList& bases,
	   const Storage& storage) :
  VirtualBase(new value_type[product(extents)], extents, bases, storage,
	      computeStrides(extents, storage)),
  // This does nothing.
  Base(_data, extents, bases, storage) {
}

// Assignment operator for other array views.
template<typename _T, std::size_t _Dimension>
template<typename _T2, typename _DataPointer2>
inline
MultiArray<_T, _Dimension>&
MultiArray<_T, _Dimension>::
operator=(const MultiArrayConstView<_T2, Dimension, _DataPointer2>& other) {
  Base::operator=(other);
  return *this;
}

// Assignment operator for arrays with contiguous memory.
template<typename _T, std::size_t _Dimension>
template<typename _T2, typename _DataPointer2>
inline
MultiArray<_T, _Dimension>&
MultiArray<_T, _Dimension>::
operator=(const MultiArrayConstRef<_T2, Dimension, _DataPointer2>& other) {
  Base::operator=(other);
  return *this;
}

// Assignment operator.
template<typename _T, std::size_t _Dimension>
inline
MultiArray<_T, _Dimension>&
MultiArray<_T, _Dimension>::
operator=(const MultiArray& other) {
  if (this != &other) {
    Base::operator=(other);
  }
  return *this;
}

// Rebuild the data structure. Re-allocate memory if the size changes.
template<typename _T, std::size_t _Dimension>
inline
void
MultiArray<_T, _Dimension>::
rebuild(const SizeList& extents, const IndexList& bases,
	const Storage& storage) {
  const size_type newSize = product(extents);
  if (newSize == size()) {
    Base::rebuild(data(), extents, bases, storage);
  }
  else {
    Base::rebuild(new value_type[newSize], extents, bases, storage);
  }
}

END_NAMESPACE_ARRAY
