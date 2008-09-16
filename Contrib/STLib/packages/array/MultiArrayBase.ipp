// -*- C++ -*-

#if !defined(__array_MultiArrayBase_ipp__)
#error This file is an implementation detail of the class MultiArrayBase.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Construct from the array extents, the index bases, the storage order, and 
// the strides.
template<std::size_t _Dimension>
inline
MultiArrayBase<_Dimension>::
MultiArrayBase(const SizeList& extents, const IndexList& bases,
	       const Storage& storage, const IndexList& strides) :
  _extents(extents),
  _bases(bases),
  _storage(storage),
  _strides(strides),
  _offset(dot(_strides, _bases)),
  _size(product(_extents)) {
}

// Rebuild the data structure.
template<std::size_t _Dimension>
inline
void
MultiArrayBase<_Dimension>::
rebuild(const SizeList& extents, const IndexList& bases,
	const Storage& storage, const IndexList& strides) {
  _extents = extents;
  _bases = bases;
  _storage = storage;
  _strides = strides;
  _offset = dot(_strides, _bases);
  _size = product(_extents);
}

END_NAMESPACE_ARRAY
