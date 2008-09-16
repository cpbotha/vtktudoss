// -*- C++ -*-

#if !defined(__stochastic_api_ipp__)
#error This file is an implementation detail.
#endif

BEGIN_NAMESPACE_STOCHASTIC

// Build a new solver that uses the direct method.
template<typename _Solver>
inline
_Solver*
newSolverDirect(const std::size_t numberOfSpecies,
		const std::size_t numberOfReactions,
		const int packedReactions[], 
		const double propensityFactors[]) {
  typedef typename _Solver::State State;
  typedef typename _Solver::PropensitiesFunctor PropensitiesFunctor;
  typedef typename _Solver::ExponentialGenerator ExponentialGenerator;
  typedef typename _Solver::DiscreteFiniteGenerator DiscreteFiniteGenerator;
  typedef typename _Solver::Number Number;
  typedef ReactionSet<Number> ReactionSet;
  typedef typename ReactionSet::Reaction Reaction;

  // Construct the reaction set.
  std::vector<Reaction> reactionsVector(numberOfReactions);
  {
    std::vector<int> reactantIndices, reactantStoichiometries,
      productIndices, productStoichiometries;
    const int* data = packedReactions;
    for (std::size_t i = 0; i != numberOfReactions; ++i) {
      reactantIndices.clear();
      reactantStoichiometries.clear();
      productIndices.clear();
      productStoichiometries.clear();

      int numberOfReactants = *data++;
      for (int j = 0; j != numberOfReactants; ++j) {
	reactantIndices.push_back(*data++);
	reactantStoichiometries.push_back(*data++);
      }
      int numberOfProducts = *data++;
      for (int j = 0; j != numberOfProducts; ++j) {
	productIndices.push_back(*data++);
	productStoichiometries.push_back(*data++);
      }
      reactionsVector[i] = 
	Reaction(reactantIndices.begin(), reactantIndices.end(),
		 reactantStoichiometries.begin(), reactantStoichiometries.end(),
		 productIndices.begin(), productIndices.end(),
		 productStoichiometries.begin(), productStoichiometries.end(),
		 propensityFactors[i]);
    }
  }
  ReactionSet reactions(reactionsVector.begin(), reactionsVector.end());
  
  //
  // Build the state change vectors.
  //

  typename State::ScvContainer stateChangeVectors;
  buildStateChangeVectors(numberOfSpecies, reactions.getBeginning(),
			  reactions.getEnd(), &stateChangeVectors);

  //
  // Build the array of reaction influences.
  //
  
  ads::StaticArrayOfArrays<int> reactionInfluence;
  computeReactionPropensityInfluence
    (numberOfSpecies, reactions.getBeginning(), reactions.getEnd(),
     &reactionInfluence, true);

  // The propensities functor.
  PropensitiesFunctor propensitiesFunctor(reactions);

  // Construct the state class.
  State state(numberOfSpecies, stateChangeVectors);

  // Check the validity of the initial state.
  assert(state.isValid());

  // Construct the simulation class.
  return new _Solver(state, propensitiesFunctor, reactionInfluence);
}


// Delete the solver.
template<typename _Solver>
inline
void
deleteSolver(_Solver* solver) {
  assert(solver);
  delete solver;
  solver = 0;
}

// Generate the state vector for the Mersenne Twister from the seed.
// Return a new seed.
inline
unsigned
generateMt19937State(const unsigned seed, unsigned state[]) {
  return numerical::DiscreteUniformGeneratorMt19937::generateState(seed, state);
}

// Get the state of the Mersenne twister.
template<typename _Solver>
inline
void
getMt19937State(const _Solver* const solver, unsigned state[]) {
  solver->getDiscreteUniformGenerator().getState(state);
}

// Set the state of the Mersenne twister.
template<typename _Solver>
inline
void
setMt19937State(_Solver* const solver, const unsigned state[]) {
  solver->getDiscreteUniformGenerator().setState(state);
}

// CONTINUE: Think about the return value.
// Generate a trajectory.
template<typename _Solver>
inline
int
generateTrajectory(_Solver* const solver, const int initialPopulationsArray[],
		   const double startTime, std::size_t maximumAllowedSteps,
		   const std::size_t numberOfFrames, const double frameTimes[],
		   int framePopulations[], std::size_t frameReactionCounts[]) {
  typedef typename _Solver::Number Number;

  // CONTINUE
  //Py_BEGIN_ALLOW_THREADS

  const std::size_t numberOfSpecies = solver->getState().getNumberOfSpecies();
  const std::size_t numberOfReactions = 
    solver->getState().getNumberOfReactions();
  
  // If they did not specify a maximum allowed number of steps.
  if (maximumAllowedSteps == 0) {
    maximumAllowedSteps = std::numeric_limits<std::size_t>::max();
  }

  // Copy the initial populations.
  typename _Solver::State::PopulationsContainer
    initialPopulations(numberOfSpecies);
  std::copy(initialPopulationsArray, initialPopulationsArray + numberOfSpecies, 
	    initialPopulations.begin());
  // Set the initial population and the starting time.
  solver->initialize(initialPopulations, startTime);

  // Run the simulation.
  int* populationIterator = framePopulations;
  std::size_t* reactionCountIterator = frameReactionCounts;
  for (std::size_t i = 0; i != numberOfFrames; ++i) {
    // Make a termination condition.
    EssTerminationConditionEndTimeReactionCount<Number>
      terminationCondition(frameTimes[i], maximumAllowedSteps);
    // Simulate up to the termination condition.
    solver->simulate(terminationCondition);
    // Record the populations and reaction counts.
    for (std::size_t species = 0; species != numberOfSpecies; ++species) {
      *populationIterator++ = solver->getState().getPopulation(species);
    }
    for (std::size_t reaction = 0; reaction != numberOfReactions; ++reaction) {
      *reactionCountIterator++ = solver->getState().getReactionCount(reaction);
    }
  }

  // CONTINUE
  //Py_END_ALLOW_THREADS

  return true;
}


inline
int
simulate(int numberOfSpecies, int initialPopulationsArray[],
	 int numberOfReactions, int packedReactions[], 
	 double propensityFactors[],
	 double startTime, std::size_t maximumAllowedSteps,
	 int numberOfFrames, double frameTimes[],
	 int framePopulations[], std::size_t frameReactionCounts[],
	 unsigned mt19937state[]) {
  typedef double Number;
  typedef numerical::ExponentialGeneratorZiggurat<> ExponentialGenerator;
  typedef numerical::DiscreteFiniteGeneratorRejectionBinsSplitting<true, true>
    DiscreteFiniteGenerator;

  typedef State<> State;

  typedef ReactionSet<Number> ReactionSet;
  typedef ReactionSet::Reaction Reaction;

  typedef PropensitiesSingle<Number> PropensitiesFunctor;

  const int N = numerical::mt19937mn::N;
  assert(maximumAllowedSteps <= double(std::numeric_limits<long>::max()));
  const std::size_t maxSteps = (maximumAllowedSteps > 0 ? 
				std::size_t(maximumAllowedSteps) :
				std::numeric_limits<std::size_t>::max());

  // Record the initial populations.
  State::PopulationsContainer initialPopulations(numberOfSpecies);
  std::copy(initialPopulationsArray, initialPopulationsArray + numberOfSpecies,
	    initialPopulations.begin());

  // Construct the reaction set.
  std::vector<Reaction> reactionsVector(numberOfReactions);
  {
    std::vector<int> reactantIndices, reactantStoichiometries,
      productIndices, productStoichiometries;
    const int* data = packedReactions;
    for (int i = 0; i != numberOfReactions; ++i) {
      reactantIndices.clear();
      reactantStoichiometries.clear();
      productIndices.clear();
      productStoichiometries.clear();

      int numberOfReactants = *data++;
      for (int j = 0; j != numberOfReactants; ++j) {
	reactantIndices.push_back(*data++);
	reactantStoichiometries.push_back(*data++);
      }
      int numberOfProducts = *data++;
      for (int j = 0; j != numberOfProducts; ++j) {
	productIndices.push_back(*data++);
	productStoichiometries.push_back(*data++);
      }
      reactionsVector[i] = 
	Reaction(reactantIndices.begin(), reactantIndices.end(),
		 reactantStoichiometries.begin(), reactantStoichiometries.end(),
		 productIndices.begin(), productIndices.end(),
		 productStoichiometries.begin(), productStoichiometries.end(),
		 propensityFactors[i]);
    }
  }
  ReactionSet reactions(reactionsVector.begin(), reactionsVector.end());
  
  //
  // Build the state change vectors.
  //

  State::ScvContainer stateChangeVectors;
  buildStateChangeVectors(initialPopulations.size(), reactions.getBeginning(),
			  reactions.getEnd(), &stateChangeVectors);

  //
  // Build the array of reaction influences.
  //
  
  ads::StaticArrayOfArrays<int> reactionInfluence;
  computeReactionPropensityInfluence
    (initialPopulations.size(), reactions.getBeginning(), reactions.getEnd(),
     &reactionInfluence, true);

  //
  // The propensities functor.
  //

  PropensitiesFunctor propensitiesFunctor(reactions);


  // Construct the state class.
  State state(initialPopulations, stateChangeVectors);

  // Check the validity of the initial state.
  if (! state.isValid()) {
    std::cerr << "Error: The initial state of the simulation is not valid.\n";
    return false;
  }

  // Construct the simulation class.
  Direct<DiscreteFiniteGenerator, ExponentialGenerator, PropensitiesFunctor,
    State> direct(state, propensitiesFunctor, reactionInfluence);

  // If a state vector was specified, instead of a null pointer.
  if (mt19937state) {
    // Set the state of the Mersenne twister.
    direct.getDiscreteUniformGenerator().setState(mt19937state);
  }

  // Run the simulation.
  direct.initialize(initialPopulations, startTime);
  int* populationIterator = framePopulations;
  std::size_t* reactionCountIterator = frameReactionCounts;
  for (int i = 0; i != numberOfFrames; ++i) {
    // Make a termination condition.
    EssTerminationConditionEndTimeReactionCount<Number>
      terminationCondition(frameTimes[i], maxSteps);
    // Simulate up to the termination condition.
    direct.simulate(terminationCondition);
    // Record the populations and reaction counts.
    for (int species = 0; species != numberOfSpecies; ++species) {
      *populationIterator++ = direct.getState().getPopulation(species);
    }
    for (int reaction = 0; reaction != numberOfReactions; ++reaction) {
      *reactionCountIterator++ = direct.getState().getReactionCount(reaction);
    }
  }

  // If a state vector was specified, instead of a null pointer.
  if (mt19937state) {
    // Get the new state of the Mersenne twister.
    for (int i = 0; i != N; ++i) {
      mt19937state[i] = direct.getDiscreteUniformGenerator().getState(i);
    }
  }

  return true;
}

inline
int
simulate(int numberOfSpecies, int initialPopulationsArray[],
	 int numberOfReactions, int packedReactions[], 
	 double propensityFactors[],
	 double startTime, std::size_t maximumAllowedSteps,
	 int numberOfFrames, double frameTimes[],
	 int framePopulations[], std::size_t frameReactionCounts[]) {
  return simulate(numberOfSpecies, initialPopulationsArray,
		  numberOfReactions, packedReactions, 
		  propensityFactors, startTime, maximumAllowedSteps,
		  numberOfFrames, frameTimes,
		  framePopulations, frameReactionCounts, 0);
}

END_NAMESPACE_STOCHASTIC

