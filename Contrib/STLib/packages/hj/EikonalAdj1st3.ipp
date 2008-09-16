// -*- C++ -*-

#if !defined(__hj_EikonalAdj1st3_ipp__)
#error This file is an implementation detail of the class EikonalAdj1st.
#endif

BEGIN_NAMESPACE_HJ

template<typename T>
inline
typename EikonalAdj1st<3,T>::Number
EikonalAdj1st<3,T>::
diff_adj(const Index& i, const Index& d) const {
#ifdef DEBUG_EikonalAdj1st
  assert(debug::is_adjacent(d));
#endif

  return diff_a1(_solution(i[0] + d[0], i[1] + d[1], i[2] + d[2]), 
		 _inverseSpeed(i));
}


template<typename T>
inline
typename EikonalAdj1st<3,T>::Number
EikonalAdj1st<3,T>::
diff_adj_adj(const Index& i,	const Index& a, const Index& b) const {
#ifdef DEBUG_EikonalAdj1st
  assert(debug::is_adjacent(a));
  assert(debug::is_adjacent(b));
#endif

  if (_status(i[0] + b[0], i[1] + b[1], i[2] + b[2]) == KNOWN) {
    return diff_a1_a1(_solution(i[0] + a[0], i[1] + a[1], i[2] + a[2]), 
		      _solution(i[0] + b[0], i[1] + b[1], i[2] + b[2]), 
		      _inverseSpeed(i));
  }
  return std::numeric_limits<Number>::max();
}


template<typename T>
inline
typename EikonalAdj1st<3,T>::Number
EikonalAdj1st<3,T>::
diff_adj_adj_adj(const Index& i, const Index& a, const Index& b,
		 const Index& c) const {
#ifdef DEBUG_EikonalAdj1st
  assert(debug::is_adjacent(a));
  assert(debug::is_adjacent(b));
  assert(debug::is_adjacent(c));
#endif

  if (_status(i[0] + b[0], i[1] + b[1], i[2] + b[2]) == KNOWN &&
      _status(i[0] + c[0], i[1] + c[1], i[2] + c[2]) == KNOWN) {
    return diff_a1_a1_a1(_solution(i[0] + a[0], i[1] + a[1], i[2] + a[2]), 
			 _solution(i[0] + b[0], i[1] + b[1], i[2] + b[2]), 
			 _solution(i[0] + c[0], i[1] + c[1], i[2] + c[2]), 
			 _inverseSpeed(i));
  }
  return std::numeric_limits<Number>::max();
}

END_NAMESPACE_HJ

// End of file.
