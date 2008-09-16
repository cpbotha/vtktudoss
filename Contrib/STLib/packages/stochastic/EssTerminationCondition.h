// -*- C++ -*-

/*! 
  \file stochastic/EssTerminationCondition.h
  \brief Termination conditions form exact stochastic simulations.
*/

#if !defined(__stochastic_EssTerminationCondition_h__)
#define __stochastic_EssTerminationCondition_h__

#include "State.h"

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_EssTerminationCondition)
#define DEBUG_stochastic_EssTerminationCondition
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Terminate when no more reactions can fire.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename _T = double>
class EssTerminationConditionExhaust {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _T Number;

  //
  // Default constructors will do.
  //

  //--------------------------------------------------------------------------
  //! \name Termination condition.
  //@{
public:

  //! Return false.
  /*!
    \param state The current state.

    \return false.
  */
  template<class _State>
  bool
  operator()(const _State& state) const {
    return false;
  }

  //! Return false.
  /*!
    \param state The current state.
    \param timeToNextReaction The time increment to the next reaction.

    \return false.
  */
  template<class _State>
  bool
  operator()(const _State& state, const Number timeToNextReaction) const {
    return false;
  }

  //! Return infinity.
  Number
  getEndTime() const {
    return std::numeric_limits<Number>::max();
  }

  //@}
};



  
//! Terminate when the end time is reached.
/*!
  \param _T The number type.  The default is \c double .
*/
template<typename _T = double>
class EssTerminationConditionEndTime {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _T Number;

  //
  // Member data.
  //
private:

  Number _endTime;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssTerminationConditionEndTime();
  //! Assignment operator not implemented.
  EssTerminationConditionEndTime&
  operator=(const EssTerminationConditionEndTime&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the end time.
  EssTerminationConditionEndTime(const Number endTime) :
    _endTime(endTime)
  {}

  //! Copy constructor.
  EssTerminationConditionEndTime(const EssTerminationConditionEndTime& other) :
    _endTime(other._endTime)
  {}

  //! Destructor.
  ~EssTerminationConditionEndTime() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Termination condition.
  //@{
public:

  //! Return true if the simulation should terminate before the next step.
  /*!
    \param state The current state.

    \return True if taking the next step would exceed the end time for 
    the simulation.
  */
  template<class _State>
  bool
  operator()(const _State& state) const {
    return state.getTime() >= _endTime;
  }

  //! Return true if the simulation should terminate before the next reaction.
  /*!
    \param state The current state.
    \param timeToNextReaction The time increment to the next reaction.

    \return True if firing the next reaction would exceed the end time for 
    the simulation.
  */
  template<class _State>
  bool
  operator()(const _State& state, const Number timeToNextReaction) const {
    return state.getTime() + timeToNextReaction > _endTime;
  }

  //! Return the end time.
  Number
  getEndTime() const {
    return _endTime;
  }

  //@}
};



  
//! Terminate when a specified number of reactions have fired.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename _T = double>
class EssTerminationConditionReactionCount {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _T Number;

  //
  // Member data.
  //
private:

  std::size_t _reactionLimit;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssTerminationConditionReactionCount();
  //! Assignment operator not implemented.
  EssTerminationConditionReactionCount&
  operator=(const EssTerminationConditionReactionCount&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the reaction limit.
  EssTerminationConditionReactionCount(const std::size_t reactionLimit) :
    _reactionLimit(reactionLimit)
  {}

  //! Copy constructor.
  EssTerminationConditionReactionCount
  (const EssTerminationConditionReactionCount& other) :
    _reactionLimit(other._reactionLimit) {
  }

  //! Destructor.
  ~EssTerminationConditionReactionCount() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Termination condition.
  //@{
public:

  //! Return true if the simulation should terminate before the next step.
  /*!
    \param state The current state.

    \return True if taking another step would exceed the reaction count
    limit.
  */
  template<class _State>
  bool
  operator()(const _State& state) const {
    return state.getReactionCount() >= _reactionLimit;
  }

  //! Return true if the simulation should terminate before the next reaction.
  /*!
    \param state The current state.
    \param timeToNextReaction The time increment to the next reaction.

    \note The time to the next reaction is not used.  It is a function
    parameter for constistency with the other termination conditions.
    
    \return True if firing the next reaction would exceed the reaction count
    limit.
  */
  template<class _State>
  bool
  operator()(const _State& state, const Number timeToNextReaction) const {
    return state.getReactionCount() >= _reactionLimit;
  }

  //! Return infinity.
  Number
  getEndTime() const {
    return std::numeric_limits<Number>::max();
  }

  //@}
};

  
//! Terminate when the end time is reached or when a specified number of reactions have fired.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename _T = double>
class EssTerminationConditionEndTimeReactionCount {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef _T Number;

  //
  // Member data.
  //
private:

  Number _endTime;
  std::size_t _reactionLimit;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssTerminationConditionEndTimeReactionCount();
  //! Assignment operator not implemented.
  EssTerminationConditionEndTimeReactionCount&
  operator=(const EssTerminationConditionEndTimeReactionCount&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the end time and the reaction limit.
  EssTerminationConditionEndTimeReactionCount
  (const Number endTime, const std::size_t reactionLimit) :
    _endTime(endTime),
    _reactionLimit(reactionLimit)
  {}

  //! Copy constructor.
  EssTerminationConditionEndTimeReactionCount
  (const EssTerminationConditionEndTimeReactionCount& other) :
    _endTime(other._endTime),
    _reactionLimit(other._reactionLimit) {
  }

  //! Destructor.
  ~EssTerminationConditionEndTimeReactionCount() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Termination condition.
  //@{
public:

  //! Return true if the simulation should terminate before the next step.
  /*!
    \param state The current state.

    \return True if taking another step would exceed the end time or the 
    reaction count limit.
  */
  template<class _State>
  bool
  operator()(const _State& state) const {
    return state.getTime() >= _endTime ||
      state.getReactionCount() >= _reactionLimit;
  }

  //! Return true if the simulation should terminate before the next reaction.
  /*!
    \param state The current state.
    \param timeToNextReaction The time increment to the next reaction.

    \return True if firing the next reaction would exceed the end time or the 
    reaction count limit.
  */
  template<class _State>
  bool
  operator()(const _State& state, const Number timeToNextReaction) const {
    return state.getTime() + timeToNextReaction > _endTime ||
      state.getReactionCount() >= _reactionLimit;
  }

  //! Return the end time.
  Number
  getEndTime() const {
    return _endTime;
  }

  //@}
};

  
END_NAMESPACE_STOCHASTIC

#endif
