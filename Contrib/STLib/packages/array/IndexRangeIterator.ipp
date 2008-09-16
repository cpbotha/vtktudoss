// -*- C++ -*-

#if !defined(__array_IndexRangeIterator_ipp__)
#error This file is an implementation detail of the class IndexRangeIterator.
#endif

BEGIN_NAMESPACE_ARRAY

//--------------------------------------------------------------------------
// Constructors etc.

// Return an iterator to the beginning of the index range.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>
IndexRangeIterator<_Dimension>::
begin(const Range& range) {
  IndexRangeIterator x;
  x._indexList = range.bases();
  x._rank = 0;
  x._range = range;
  return x;
}

// Return an iterator to the end of the index range.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>
IndexRangeIterator<_Dimension>::
end(const Range& range) {
  IndexRangeIterator x;
  x._indexList = range.bases();
  x._indexList[0] = range.bases()[0] + range.extents()[0] * range.steps()[0];
  x._rank = product(range.extents());
  x._range = range;
  return x;
}

// Copy constructor.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>::
IndexRangeIterator(const IndexRangeIterator& other) :
  _indexList(other._indexList),
  _rank(other._rank),
  _range(other._range) {
}

// Assignment operator.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>&
IndexRangeIterator<_Dimension>::
operator=(const IndexRangeIterator& other) {
  if (this != &other) {
    _indexList = other._indexList;
    _rank = other._rank;
    _range = other._range;
  }
  return *this;
}

//--------------------------------------------------------------------------
// Validity.

// Return true if the iterator is valid.
template<std::size_t _Dimension>
inline
bool
IndexRangeIterator<_Dimension>::
isValid() const {
  if(! (0 <= _rank && _rank < Index(product(_range.extents())))) {
    return false;
  }
  for (size_type i = 0; i != Dimension; ++i) {
    if(! (_range.bases()[i] <= _indexList[i] && 
	  _indexList[i] < _range.bases()[i] + 
	  Index(_range.extents()[i]) * _range.steps()[i])) {
      return false;
    }
  }
  return true;
}

// Return true if the iterator is at the beginning.
template<std::size_t _Dimension>
inline
bool
IndexRangeIterator<_Dimension>::
isBegin() const {
  if (_rank == 0) {
    assert(_indexList == _range.bases());
    return true;
  }
  return false;
}

// Return true if the iterator is at the end.
template<std::size_t _Dimension>
inline
bool
IndexRangeIterator<_Dimension>::
isEnd() const {
  if (_indexList[0] != _range.bases()[0] + 
      Index(_range.extents()[0]) * _range.steps()[0]) {
    return false;
  }
  for (size_type i = 1; i != Dimension; ++i) {
    if (_indexList[i] != _range.bases()[i]) {
      return false;
    }
  }
  assert(_rank == Index(product(_range.extents())));
  return true;
}

//--------------------------------------------------------------------------
// Forward iterator requirements.

// Pre-increment.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>&
IndexRangeIterator<_Dimension>::
operator++() {
#ifdef DEBUG_array_IndexRangeIterator
  assert(isValid());
#endif

  // Increment the index list.
  _indexList[Dimension - 1] += _range.steps()[Dimension - 1];
  for (size_type i = Dimension - 1; i != 0; --i) {
    if (_indexList[i] == _range.bases()[i] + 
	Index(_range.extents()[i]) * _range.steps()[i]) {
      _indexList[i] = _range.bases()[i];
      _indexList[i - 1] += _range.steps()[i - 1];
    }
    else {
      break;
    }
  }
  // Increment the rank.
  ++_rank;

#ifdef DEBUG_array_IndexRangeIterator
  assert(isValid() || isEnd());
#endif

  return *this;
}

// Post-increment.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>
IndexRangeIterator<_Dimension>::
operator++(int) {
  IndexRangeIterator tmp(*this);
  ++*this;
  return tmp;
}

//--------------------------------------------------------------------------
// Bidirectional iterator requirements.

// Pre-decrement.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>&
IndexRangeIterator<_Dimension>::
operator--() {
#ifdef DEBUG_array_IndexRangeIterator
  assert((isValid() || isEnd()) && ! isBegin());
#endif

  // Decrement the index list.
  _indexList[Dimension - 1] -= _range.steps()[Dimension - 1];
  for (size_type i = Dimension - 1; i != 0; --i) {
    if (_indexList[i] == _range.bases()[i] - _range.steps()[i]) {
      _indexList[i] = _range.bases()[i] + 
	(_range.extents()[i] - 1) * _range.steps()[i];
      _indexList[i - 1] -= _range.steps()[i - 1];
    }
    else {
      break;
    }
  }
  // Decrement the rank.
  --_rank;

#ifdef DEBUG_array_IndexRangeIterator
  assert(isValid());
#endif

  return *this;
}


// Post-decrement.
template<std::size_t _Dimension>
inline
IndexRangeIterator<_Dimension>
IndexRangeIterator<_Dimension>::
operator--(int) {
  IndexRangeIterator tmp(*this);
  --*this;
  return tmp;
}

// Calculate the index list from the rank.
template<std::size_t _Dimension>
inline
void
IndexRangeIterator<_Dimension>::
calculateIndexList() {
#ifdef DEBUG_array_IndexRangeIterator
  assert(0 <= _rank && _rank <= Index(product(_range.extents())));
#endif
  // The strides.
  IndexList strides;
  strides[Dimension - 1] = 1;
  for (size_type i = Dimension - 1; i != 0; --i) {
    strides[i - 1] = strides[i] * _range.extents()[i];
  }
  Index r = _rank;
  // Traverse from most significant to least.
  for (std::size_t i = 0; i != Dimension; ++i) {
    _indexList[i] = r / strides[i];
    r -= _indexList[i] * strides[i];
    _indexList[i] *= _range.steps()[i];
    _indexList[i] += _range.bases()[i];
  }
#ifdef DEBUG_array_IndexRangeIterator
  assert(r == 0);
#endif
}

END_NAMESPACE_ARRAY
