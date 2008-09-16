// -*- C++ -*-

#if !defined(__stochastic_StateVariables_ipp__)
#error This file is an implementation detail of StateVariables.
#endif

BEGIN_NAMESPACE_STOCHASTIC

//----------------------------------------------------------------------------
// Free functions.
//----------------------------------------------------------------------------

// Return true if the states are equal.
template<typename _PType, typename _RCType>
inline
bool
operator==(const StateVariables<_PType, _RCType>& x,
	   const StateVariables<_PType, _RCType>& y) {
  // Check the time.
  if (x.getTime() != y.getTime()) {
    return false;
  }

  // Check the populations.
  if (x.getPopulations() != y.getPopulations()) {
    return false;
  }

  return true;
}


// Write the populations in ascii format.
template<typename _PType, typename _RCType>
inline
void
writePopulationsAscii(std::ostream& out,
		      const StateVariables<_PType, _RCType>& x) {
  out << x.getPopulations();
}
 

// Read the populations in ascii format.
template<typename _PType, typename _RCType>
inline
void
readPopulationsAscii(std::istream& in, StateVariables<_PType, _RCType>* state) {
  typename StateVariables<_PType, _RCType>::PopulationsContainer populations;
  in >> populations;
  state->setPopulations(populations);
}


// Write the state in ascii format.
template<typename _PType, typename _RCType>
inline
void
writeAscii(std::ostream& out, const StateVariables<_PType, _RCType>& state) {
  // The time.
  out << state.getTime() << "\n";
  out << state.getPopulations();
}
 
END_NAMESPACE_STOCHASTIC
