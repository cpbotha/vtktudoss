// -*- C++ -*-

#if !defined(__hj_EikonalAdjDiag1st2_ipp__)
#error This file is an implementation detail of the class EikonalAdjDiag1st.
#endif

BEGIN_NAMESPACE_HJ

template<typename T>
inline
typename EikonalAdjDiag1st<2,T>::Number
EikonalAdjDiag1st<2,T>::
diff_adj(const Index& i, const Index& d) const {
#ifdef DEBUG_EikonalAdjDiag1st
  assert(debug::is_adjacent(d));
#endif

  return diff_a1(_solution(i, d), _inverseSpeed(i));
}


template<typename T>
inline
typename EikonalAdjDiag1st<2,T>::Number
EikonalAdjDiag1st<2,T>::
diff_adj_diag(const Index& i, const Index& a, const Index& b) const {
#ifdef DEBUG_EikonalAdjDiag1st
  assert(debug::is_adjacent(a) && debug::is_diagonal(b));
#endif

  if (_status(i, b) == KNOWN) {
    return diff_a1_d1(_solution(i, a), 
		      _solution(i, b), 
		      _inverseSpeed(i));
  }
  return std::numeric_limits<Number>::max();
}


template<typename T>
inline
typename EikonalAdjDiag1st<2,T>::Number
EikonalAdjDiag1st<2,T>::
diff_diag(const Index& i, const Index& d) const {
#ifdef DEBUG_EikonalAdjDiag1st
  assert(debug::is_diagonal(d));
#endif

  return diff_d1(_solution(i, d), _inverseSpeed(i));
}


template<typename T>
inline
typename EikonalAdjDiag1st<2,T>::Number
EikonalAdjDiag1st<2,T>::
diff_diag_adj(const Index& i, const Index& a, const Index& b) const {
#ifdef DEBUG_EikonalAdjDiag1st
  assert(debug::is_diagonal(a) && debug::is_adjacent(b));
#endif

  if (_status(i, b) == KNOWN) {
    return diff_a1_d1(_solution(i, b), 
		      _solution(i, a), 
		      _inverseSpeed(i));
  }
  return std::numeric_limits<Number>::max();
}

END_NAMESPACE_HJ

// End of file.
