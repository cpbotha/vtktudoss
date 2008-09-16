// -*- C++ -*-

#if !defined(__array_Array_ipp__)
#error This file is an implementation detail of array.
#endif

BEGIN_NAMESPACE_ARRAY

//----------------------------------------------------------------------------
// I/O
//----------------------------------------------------------------------------

template<typename _T, std::size_t _N>
inline
std::ostream&
operator<<(std::ostream& out, const Array<_T, _N>& x) {
  typename Array<_T, _N>::const_iterator i = x.begin(), i_end = x.end();
  if (i != i_end) {
    out << *i;
    ++i;
  }
  for (; i != i_end; ++i) {
    out << " " << *i;
  }
  return out;
}

template<typename _T, std::size_t _N>
inline
std::istream&
operator>>(std::istream& in, Array<_T, _N>& x) {
  typename Array<_T, _N>::iterator i = x.begin(), i_end = x.end();
  for (; i != i_end; ++i) {
    in >> *i;
  }
  return in;
}

//----------------------------------------------------------------------------
// Math operators.
//----------------------------------------------------------------------------

// Return the sum of the components.
template<typename _T, std::size_t _N>
inline
_T
sum(const Array<_T, _N>& x) {
  return std::accumulate(x.begin(), x.end(), _T(0));
}

// Return the product of the components.
template<typename _T, std::size_t _N>
inline
_T
product(const Array<_T, _N>& x) {
  return std::accumulate(x.begin(), x.end(), _T(1), std::multiplies<_T>());
}

// Return the minimum component.  Use < for comparison.
template<typename _T, std::size_t _N>
inline
_T
min(const Array<_T, _N>& x) {
  BOOST_STATIC_ASSERT(_N > 0);
  return *std::min_element(x.begin(), x.end());
}
    
// Return the maximum component.  Use > for comparison.
template<typename _T, std::size_t _N>
inline
_T
max(const Array<_T, _N>& x) {
  BOOST_STATIC_ASSERT(_N > 0);
  return *std::max_element(x.begin(), x.end());
}
    
// Return the dot product of the two vectors.
template<typename _T1, typename _T2, std::size_t _N>
inline
_T1
dot(const Array<_T1, _N>& x, const Array<_T2, _N>& y) {
  return std::inner_product(x.begin(), x.end(), y.begin(), _T1(0));
}

//----------------------------------------------------------------------------
// I/O
//----------------------------------------------------------------------------
    
// Write the array elements in binary format.
template<typename _T, std::size_t _N>
inline
void
writeElementsBinary(std::ostream& out, const Array<_T, _N>& x) {
  out.write(reinterpret_cast<const char*>(x.data()), _N * sizeof(_T));
}
    
// Read the array elements in binary format.
template<typename _T, std::size_t _N>
inline
void
readElementsBinary(std::istream& in, Array<_T, _N>& x) {
  in.read(reinterpret_cast<char*>(x.data()), _N * sizeof(_T));
}

END_NAMESPACE_ARRAY
