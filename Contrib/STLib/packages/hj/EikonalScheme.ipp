// -*- C++ -*-

#if !defined(__hj_EikonalScheme_ipp__)
#error This file is an implementation detail of the class EikonalScheme.
#endif

BEGIN_NAMESPACE_HJ

//
// File I/O
//

template<typename T>
inline
void
print_inverse_speed_array(std::ostream& out, 
			  const ads::Array<2,T>& inverseSpeed) {
  out << "Inverse speed:" << '\n';
  for (int j = inverseSpeed.ubound(1) - 1; j >= inverseSpeed.lbound(1); --j){
    for (int i = inverseSpeed.lbound(0); i < inverseSpeed.ubound(0); ++i){
      out << inverseSpeed(i, j) << " ";
    }
    out << '\n';
  }
}


template<typename T>
inline
void
print_inverse_speed_array(std::ostream& out, 
			  const ads::Array<3,T>& inverseSpeed) {
  out << "Inverse speed:" << '\n';
  for (int k = inverseSpeed.ubound(2) - 1; k >= inverseSpeed.lbound(2); --k){
    for (int j = inverseSpeed.ubound(1) - 1; j >= inverseSpeed.lbound(1); --j){
      for (int i = inverseSpeed.lbound(0); i < inverseSpeed.ubound(0); ++i){
	out << inverseSpeed(i, j, k) << " ";
      }
      out << '\n';
    }
    out << '\n';
  }
}


template<int N, typename T>
inline
void
EikonalScheme<N,T>::
put(std::ostream& out) const {
  print_inverse_speed_array(out, _inverseSpeed);
}

END_NAMESPACE_HJ

// End of file.
