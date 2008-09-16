// -*- C++ -*-

#if !defined(__array_IndexListIterator_ipp__)
#error This file is an implementation detail of the class IndexListIterator.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Return an iterator to the beginning of the index range.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>
IndexListIterator<_Dimension>::
begin(const Array& array) {
  IndexListIterator x;
  x._indexList = array.bases();
  x._rank = 0;
  x._array = &array;
  return x;
}

// Return an iterator to the end of the index range.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>
IndexListIterator<_Dimension>::
end(const Array& array) {
  IndexListIterator x;
  x._indexList = array.bases();
  const size_type n = array.storage()[Dimension - 1];
  x._indexList[n] = array.bases()[n] + array.extents()[n];
  x._rank = array.size();
  x._array = &array;
  return x;
}

// Copy constructor.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>::
IndexListIterator(const IndexListIterator& other) :
  _indexList(other._indexList),
  _rank(other._rank),
  _array(other._array) {
}

// Assignment operator.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>&
IndexListIterator<_Dimension>::
operator=(const IndexListIterator& other) {
  if (this != &other) {
    _indexList = other._indexList;
    _rank = other._rank;
    _array = other._array;
  }
  return *this;
}

//--------------------------------------------------------------------------
// Validity.

// Return true if the iterator is valid.
template<std::size_t _Dimension>
inline
bool
IndexListIterator<_Dimension>::
isValid() const {
  if(! (0 <= _rank && _rank < Index(_array->size()))) {
    return false;
  }
  for (size_type i = 0; i != Dimension; ++i) {
    if(! (_array->bases()[i] <= _indexList[i] && 
	  _indexList[i] < _array->bases()[i] + 
	  Index(_array->extents()[i]))) {
      return false;
    }
  }
  return true;
}

// Return true if the iterator is at the beginning.
template<std::size_t _Dimension>
inline
bool
IndexListIterator<_Dimension>::
isBegin() const {
  return _indexList == _array->bases();
}

// Return true if the iterator is at the end.
template<std::size_t _Dimension>
inline
bool
IndexListIterator<_Dimension>::
isEnd() const {
  for (size_type i = 0; i < Dimension - 1; ++i) {
    const size_type n = _array->storage()[i];
    if (_indexList[n] != _array->bases()[n]) {
      return false;
    }
  }
  const size_type n = _array->storage()[Dimension - 1];
  if (_indexList[n] != _array->bases()[n] + 
      Index(_array->extents()[n])) {
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------
// Forward iterator requirements.

// Pre-increment.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>&
IndexListIterator<_Dimension>::
operator++() {
#ifdef DEBUG_array_IndexListIterator
  assert(isValid());
#endif

  // Increment the index list.
  ++_indexList[_array->storage()[0]];
  for (size_type i = 0; i < Dimension - 1; ++i) {
    const size_type n = _array->storage()[i];
    if (_indexList[n] == _array->bases()[n] + 
	Index(_array->extents()[n])) {
      _indexList[n] = _array->bases()[n];
      ++_indexList[_array->storage()[i + 1]];
    }
    else {
      break;
    }
  }
  // Increment the rank.
  ++_rank;

#ifdef DEBUG_array_IndexListIterator
  assert(isValid() || isEnd());
#endif

  return *this;
}

// Post-increment.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>
IndexListIterator<_Dimension>::
operator++(int) {
  IndexListIterator tmp(*this);
  ++*this;
  return tmp;
}

//--------------------------------------------------------------------------
// Bidirectional iterator requirements.

// Pre-decrement.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>&
IndexListIterator<_Dimension>::
operator--() {
#ifdef DEBUG_array_IndexListIterator
  assert((isValid() || isEnd()) && ! isBegin());
#endif

  // Decrement the index list.
  --_indexList[_array->storage()[0]];
  for (size_type i = 0; i < Dimension - 1; ++i) {
    const size_type n = _array->storage()[i];
    if (_indexList[n] == _array->bases()[n] - 1) {
      _indexList[n] = _array->bases()[n] + _array->extents()[n] - 1;
      --_indexList[_array->storage()[i + 1]];
    }
    else {
      break;
    }
  }
  // Decrement the rank.
  --_rank;

#ifdef DEBUG_array_IndexListIterator
  assert(isValid());
#endif

  return *this;
}


// Post-decrement.
template<std::size_t _Dimension>
inline
IndexListIterator<_Dimension>
IndexListIterator<_Dimension>::
operator--(int) {
  IndexListIterator tmp(*this);
  --*this;
  return tmp;
}

// Calculate the index list from the rank.
template<std::size_t _Dimension>
inline
void
IndexListIterator<_Dimension>::
calculateIndexList() {
#ifdef DEBUG_array_IndexListIterator
  assert(0 <= _rank && _rank <= Index(_array->size()));
#endif
  Index r = _rank;
  for (std::size_t i = 0; i != Dimension; ++i) {
    // Traverse from most significant to least.
    const std::size_t n = _array->storage()[Dimension - i - 1];
    _indexList[n] = r / _array->strides()[n];
    r -= _indexList[n] * _array->strides()[n];
    _indexList[n] += _array->bases()[n];
  }
#ifdef DEBUG_array_IndexListIterator
  assert(r == 0);
#endif
}

END_NAMESPACE_ARRAY
