// -*- C++ -*-

#if !defined(__stochastic_reactionPropensityInfluence_ipp__)
#error This file is an implementation detail of reactionPropensityInfluence.
#endif

BEGIN_NAMESPACE_STOCHASTIC

// Compute the reaction-propensity influence.
template<typename ForwardIterator>
inline
void
computeReactionPropensityInfluence(const int numberOfSpecies,
				   ForwardIterator reactionsBeginning,
				   ForwardIterator reactionsEnd,
				   ads::StaticArrayOfArrays<int>* influence,
				   const bool includeSelf) {
  // The index and coefficient of a species in a reaction.
  typedef std::map<int, int>::value_type IndexAndCoefficient;

  //
  // First build a mapping from reactants to reactions.
  //
  ads::StaticArrayOfArrays<int> reactantsToReactions;
  if (numberOfSpecies != 0) {
    // Determine the sizes for the array of arrays.
    std::vector<int> sizes(numberOfSpecies, 0);
    // For each reaction.
    for (ForwardIterator reaction = reactionsBeginning; 
	 reaction != reactionsEnd; ++reaction) {
      // For each reactant.
      for (int i = 0; i != reaction->getReactants().size(); ++i) {
	// Count the reaction.
	++sizes[reaction->getReactants().getIndex(i)];
      }
    }
   
    // Allocate memory.
    reactantsToReactions.rebuild(std::accumulate(sizes.begin(), sizes.end(), 0),
				 sizes.begin(), sizes.end());

    // Record the mapping from reactants to reactions.
    std::fill(sizes.begin(), sizes.end(), 0);
    int reactionIndex = 0;
    for (ForwardIterator reaction = reactionsBeginning; 
	 reaction != reactionsEnd; ++reaction, ++reactionIndex) {
      // For each reactant.
      for (int i = 0; i != reaction->getReactants().size(); ++i) {
	const int reactant = reaction->getReactants().getIndex(i);
	reactantsToReactions(reactant, sizes[reactant]) = reactionIndex;
	++sizes[reactant];
      }
    }
  }

  // The number of influences for each reaction.
  std::vector<int> sizes;
  // The reaction indices.
  std::vector<int> reactionIndices;

  // We store the reactants and products in a std::map for easy searching.
  std::map<int, int> reactants, products;
  // The influenced species and reactions for a reaction.
  std::set<int> influencedSpecies, influencedReactions;

  // For each reaction.
  int reactionIndex = 0;
  for (ForwardIterator reaction = reactionsBeginning; reaction != reactionsEnd;
       ++reaction, ++reactionIndex) {
    //
    // First determine the influence on the species.  If a reaction changes the
    // population of a species, then it influences the species.  I check 
    // the case that a species is both a reactant and product, but does 
    // not change in population.  For example, the reaction
    // x_0 + x_1 -> x_0 + x_2
    // influences species 1 and 2, but not species 0.
    //

    // Get the reactants.
    reactants.clear();
    for (int i = 0; i != reaction->getReactants().size(); ++i) {
      reactants.insert(IndexAndCoefficient(reaction->getReactants().getIndex(i),
					   reaction->getReactants()[i]));
    }
    // Get the products.
    products.clear();
    for (int i = 0; i != reaction->getProducts().size(); ++i) {
      products.insert(IndexAndCoefficient(reaction->getProducts().getIndex(i),
					  reaction->getProducts()[i]));
    }
    influencedSpecies.clear();
    influencedReactions.clear();

    // First check the reactants.
    std::map<int, int>::const_iterator j;
    // For each reactant.
    for (std::map<int, int>::const_iterator i = reactants.begin();
	 i != reactants.end(); ++i) {
      // See if the same species is also a product.
      j = products.find(i->first);
      // If the species is not a product or if the product coefficient differs
      // from the reactant coefficient.
      if (j == products.end() || j->second != i->second) {
	// The species is influenced.
	influencedSpecies.insert(i->first);
      }
    }

    // Then check the products.
    // For each product.
    for (std::map<int, int>::const_iterator i = products.begin();
	 i != products.end(); ++i) {
      // See if the same species is also a reactant.
      j = reactants.find(i->first);
      // If the species is not a reactant or if the reactant coefficient 
      // differs from the product coefficient.
      if (j == reactants.end() || j->second != i->second) {
	// The species is influenced.
	influencedSpecies.insert(i->first);
      }
    }

    // Now we have the species that are influenced by this reaction.  We then 
    // determine the reaction propensities that are influenced.

    for (std::set<int>::const_iterator species = influencedSpecies.begin();
	 species != influencedSpecies.end(); ++species) {
      for (int n = 0; n != reactantsToReactions.size(*species); ++n) {
	const int r = reactantsToReactions(*species, n);
	if (! includeSelf && r == reactionIndex) {
	  continue;
	}
	influencedReactions.insert(r);
      }
    }

    // Now we have the reactions whose propensities are influenced by 
    // this reaction.  Record the number of influenced reactions and their
    // indices.
    sizes.push_back(influencedReactions.size());
    for (std::set<int>::const_iterator i = influencedReactions.begin();
	 i != influencedReactions.end(); ++i) {
      reactionIndices.push_back(*i);
    }
  }

  // Build the static array of arrays.
  influence->rebuild(sizes.begin(), sizes.end(), reactionIndices.begin(),
		     reactionIndices.end());
}

END_NAMESPACE_STOCHASTIC

// End of file.
