// -*- C++ -*-

#if !defined(__stochastic_ReactionSet_ipp__)
#error This file is an implementation detail of ReactionSet.
#endif

BEGIN_NAMESPACE_STOCHASTIC

//----------------------------------------------------------------------------
// Manipulators.
//----------------------------------------------------------------------------


template<typename T>
template<typename _InputIterator>
inline
void
ReactionSet<T>::
rebuild(_InputIterator reactionsBeginning, _InputIterator reactionsEnd) {
  _reactions.clear();
  _reactions.insert(_reactions.end(), reactionsBeginning, reactionsEnd);
}


//----------------------------------------------------------------------------
// Free functions
//----------------------------------------------------------------------------

// Return true if the states are equal.
template<typename T>
inline
bool
operator==(const ReactionSet<T>& x, const ReactionSet<T>& y) {
  if (x.getNumberOfReactions() != y.getNumberOfReactions()) {
    return false;
  }
  for (int n = 0; n != x.getNumberOfReactions(); ++n) {
    if (x.getReaction(n) != y.getReaction(n)) {
      return false;
    }
  }

  return true;
}


// Write the reactions in ascii format.
template<typename T>
inline
void
writeAscii(std::ostream& out, const ReactionSet<T>& x) {
  out << x.getSize() << "\n";
  for (int n = 0; n != x.getSize(); ++n) {
    writeAscii(out, x.getReaction(n));
  }
}
 

// Read the reactions in ascii format.
template<typename T>
inline
void
readAscii(std::istream& in, ReactionSet<T>* x) {
  typedef typename ReactionSet<T>::Reaction Reaction;

  int numberOfReactions;
  in >> numberOfReactions;
  assert(numberOfReactions >= 0);

  std::vector<Reaction> reactions(numberOfReactions);
  for (int n = 0; n != numberOfReactions; ++n) {
    readAscii(in, &reactions[n]);
  }
  
  x->rebuild(reactions.begin(), reactions.end());
}


// Read the reactants and products in ascii format.
template<typename T>
inline
void
readReactantsAndProductsAscii(std::istream& in, ReactionSet<T>* x) {
  std::size_t numberOfReactions;
  in >> numberOfReactions;
  readReactantsAndProductsAscii(in, numberOfReactions, x);
}

// Read the reactants and products in ascii format.
template<typename T>
inline
void
readReactantsAndProductsAscii(std::istream& in, 
			      const std::size_t numberOfReactions, 
			      ReactionSet<T>* x) {
  typedef typename ReactionSet<T>::Reaction Reaction;
  std::vector<Reaction> reactions(numberOfReactions);
  for (std::size_t n = 0; n != numberOfReactions; ++n) {
    readReactantsAscii(in, &reactions[n]);
    readProductsAscii(in, &reactions[n]);
  }
  
  x->rebuild(reactions.begin(), reactions.end());
}

// Read the rate constants in ascii format.
template<typename T>
inline
void
readRateConstantsAscii(std::istream& in, ReactionSet<T>* x) {
  std::vector<T> rateConstants;
  in >> rateConstants;
  assert(rateConstants.size() == x->getSize());
  x->setRateConstants(rateConstants.begin(), rateConstants.end());
}

END_NAMESPACE_STOCHASTIC
