// -*- C++ -*-

#if !defined(__stochastic_TauLeapingDynamic_ipp__)
#error This file is an implementation detail of TauLeapingDynamic.
#endif

BEGIN_NAMESPACE_STOCHASTIC

template<class _State>
inline
void
TauLeapingDynamic<_State>::
initialize(const ReactionSet& reactionSet) {
  //
  // Compute the orders for the species.
  //
  typedef typename ReactionSet::Reaction Reaction;
  // Initialize the arrays.
  std::fill(_highestOrder.begin(), _highestOrder.end(), 0);
  std::fill(_highestIndividualOrder.begin(), _highestIndividualOrder.end(), 
	    0);

  int sum, index, order;
  // Loop over the reactions.
  for (std::size_t n = 0; n != reactionSet.getSize(); ++n) {
    const Reaction& reaction = reactionSet.getReaction(n);
    // The sum of the reactant coefficients.
    sum = ads::computeSum(reaction.getReactants());
    // Loop over the reactant species.
    for (int i = 0; i != reaction.getReactants().size(); ++i) {
      index = reaction.getReactants().getIndex(i);
      order = reaction.getReactants()[i];
      if (sum > _highestOrder[index]) {
	_highestOrder[index] = sum;
      }
      if (order > _highestIndividualOrder[index]) {
	_highestIndividualOrder[index] = order;
      }
    }
  }
  //
  // Store the reactants for each reaction.
  //
  std::vector<int> sizes, indices;
  // Loop over the reactions.
  for (std::size_t n = 0; n != reactionSet.getSize(); ++n) {
    const Reaction& reaction = reactionSet.getReaction(n);
    sizes.push_back(reaction.getReactants().size());
    for (int i = 0; i != reaction.getReactants().size(); ++i) {
      indices.push_back(reaction.getReactants().getIndex(i));
    }
  }
  _reactants.rebuild(sizes.begin(), sizes.end(), indices.begin(), 
		     indices.end());
}

template<class _State>
inline
typename TauLeapingDynamic<_State>::Number
TauLeapingDynamic<_State>::
computeStep(const std::vector<ads::SparseArray<1, PopulationType> >& 
	    listOfStateChangeVectors,
	    const std::vector<Number>& propensities,
	    const std::vector<PopulationType>& populations) {
#ifdef DEBUG_stlib
  assert(_epsilon > 0);
#endif

  // Compute mu and sigmaSquared.
  computeMuAndSigmaSquared(listOfStateChangeVectors, propensities);

  // Initialize tau to infinity.
  Number tau = std::numeric_limits<Number>::max();

  // Take the minimum over the active species.
  Number numerator, temp, a, b;
  for (std::size_t i = 0; i != _activeSpecies.size(); ++i) {
    const int n = _activeSpecies[i];
#ifdef DEBUG_stlib
    // CONTINUE
    assert(_highestOrder[n] != 0);
    assert(populations[n] > 0);
#endif
    // CONTINUE: Tau-leaping should not be used for slow reactions.
    //numerator = std::max(_epsilon * populations[n] / computeG(n), 1.0);
    numerator = _epsilon * populations[n] / computeG(n, populations[n]);

    if (_mu[n] != 0) {
      a = numerator / std::abs(_mu[n]);
    }
    else {
      a = std::numeric_limits<Number>::max();
    }

    if (_sigmaSquared[n] != 0) {
      b = numerator * numerator / _sigmaSquared[n];
    }
    else {
      b = std::numeric_limits<Number>::max();
    }

    temp = std::min(a, b);
    if (temp < tau) {
      tau = temp;
    }
  }

  return tau;
}


// Compute the g described in "Efficient step size selection for the 
// tau-leaping simulation method".
template<class _State>
inline
typename TauLeapingDynamic<_State>::Number
TauLeapingDynamic<_State>::
computeG(const int speciesIndex, const PopulationType population) const {
  const int order = _highestOrder[speciesIndex];

  if (order == 1) {
    return 1.0;
  }
  else if (order == 2) {
    if (_highestIndividualOrder[speciesIndex] == 1) {
      return 2.0;
    }
    else if (_highestIndividualOrder[speciesIndex] == 2) {
      return 2.0 + 1.0 / (population - 1);
    }
  }
  else if (order == 3) {
    if (_highestIndividualOrder[speciesIndex] == 1) {
      return 3.0;
    }
    else if (_highestIndividualOrder[speciesIndex] == 2) {
      return 1.5 * (2.0 + 1.0 / (population - 1.0));
    }
    else if (_highestIndividualOrder[speciesIndex] == 3) {
      return 3.0 + 1.0 / (population - 1) + 2.0 / (population - 2);
    }
  }

  // Catch any other cases with an assertion failure.
  assert(false);
  return 0;
}

template<class _State>
inline
void
TauLeapingDynamic<_State>::
computeMuAndSigmaSquared
(const std::vector<ads::SparseArray<1, PopulationType> >& 
 listOfStateChangeVectors,
 const std::vector<Number>& propensities) {
  // Initialize.
  std::fill(_mu.begin(), _mu.end(), 0.0);
  std::fill(_sigmaSquared.begin(), _sigmaSquared.end(), 0.0);

  Number propensity, value;
  int index, reaction;
  // Loop over the active reactions.
  for (std::size_t i = 0; i != _activeReactions.size(); ++i) {
    reaction = _activeReactions[i];
    propensity = propensities[reaction];
#ifdef DEBUG_stlib
    assert(propensity >= 0);
#endif
    const ads::SparseArray<1, PopulationType>& stateChange = 
      listOfStateChangeVectors[reaction];
    // For each species that is modified by this reaction.
    for (int j = 0; j != stateChange.size(); ++j) {
      index = stateChange.getIndex(j);
      value = stateChange[j];
      _mu[index] += value * propensity;
      _sigmaSquared[index] += value * value * propensity;
    }
  }
}

END_NAMESPACE_STOCHASTIC
