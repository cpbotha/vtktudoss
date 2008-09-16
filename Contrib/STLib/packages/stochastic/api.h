// -*- C++ -*-

#if !defined(__stochastic_api_h__)
#define __stochastic_api_h__

// CONTINUE: I used this for Py_BEGIN_ALLOW_THREADS
//#include "Python.h"

#include "Direct.h"
#include "Propensities.h"
#include "reactionPropensityInfluence.h"

// CONTINUE REMOVE
#include "../numerical/random/exponential/ExponentialGeneratorZiggurat.h"
#include "../numerical/random/discreteFinite/DiscreteFiniteGeneratorRejectionBinsSplitting.h"

#include "../ads/array/SparseArray.h"

#include <iostream>

BEGIN_NAMESPACE_STOCHASTIC

//! Build a new solver that uses the direct method.
template<typename _Solver>
_Solver*
newSolverDirect(std::size_t numberOfSpecies,
		std::size_t numberOfReactions,
		const int packedReactions[], 
		const double propensityFactors[]);

//! Delete the solver.
template<typename _Solver>
void
deleteSolver(_Solver* solver);

//! Generate the state vector for the Mersenne Twister from the seed.
/*! \return A new seed. */
unsigned
generateMt19937State(unsigned seed, unsigned state[]);

//! Get the state of the Mersenne twister.
template<typename _Solver>
void
getMt19937State(const _Solver* solver, unsigned state[]);

//! Set the state of the Mersenne twister.
template<typename _Solver>
void
setMt19937State(_Solver* solver, const unsigned state[]);

// Generate a trajectory.
template<typename _Solver>
int
generateTrajectory(_Solver* solver, const int initialPopulationsArray[],
		   double startTime, std::size_t maximumAllowedSteps,
		   std::size_t numberOfFrames, const double frameTimes[],
		   int framePopulations[], std::size_t frameReactionCounts[]);

int
simulate(int numberOfSpecies, int initialPopulationsArray[],
	 int numberOfReactions, int packedReactions[], 
	 double propensityFactors[],
	 double startTime, std::size_t maximumAllowedSteps,
	 int numberOfFrames, double frameTimes[],
	 int framePopulations[], std::size_t frameReactionCounts[],
	 unsigned mt19937state[]);

int
simulate(int numberOfSpecies, int initialPopulationsArray[],
	 int numberOfReactions, int packedReactions[], 
	 double propensityFactors[],
	 double startTime, std::size_t maximumAllowedSteps,
	 int numberOfFrames, double frameTimes[],
	 int framePopulations[], std::size_t frameReactionCounts[]);

END_NAMESPACE_STOCHASTIC

#define __stochastic_api_ipp__
#include "api.ipp"
#undef __stochastic_api_ipp__

#endif
