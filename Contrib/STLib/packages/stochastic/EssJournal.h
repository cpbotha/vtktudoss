// -*- C++ -*-

/*! 
  \file stochastic/EssJournal.h
  \brief Journals for exact stochastic simulations.
*/

#if !defined(__stochastic_EssJournalNull_h__)
#define __stochastic_EssJournalNull_h__

#include "defs.h"

#include <iostream>

// If we are debugging the whole stochastic namespace.
#if defined(DEBUG_stochastic) && !defined(DEBUG_stochastic_EssJournal)
#define DEBUG_stochastic_EssJournal
#endif

BEGIN_NAMESPACE_STOCHASTIC

//! Write the array size and the elements.
template<class _Container>
inline
void
writeArrayAscii(std::ostream& out, const _Container& array) {
  out << array.size() << "\n";
  for (std::size_t i = 0; i != array.size(); ++i) {
    out << array[i] << "\n";
  }
}


//============================================================================
//! Journal for exact stochastic simulations.  Record nothing.
/*!
  \param T The number type.  The default is \c double .

  This class has empty journalling functions.
*/
template<typename T = double>
class EssJournalNull {
  //
  // Public types.
  //
public:

  //! The number type.
  typedef T Number;

  //
  // Not implemented.
  //
private:  

  //! Copy constructor not implemented.
  EssJournalNull(const EssJournalNull&);
  //! Assignment operator not implemented.
  EssJournalNull&
  operator=(const EssJournalNull&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Default constructor.
  EssJournalNull()
  {}

  //! Construct from an output stream.
  EssJournalNull(std::ostream* out)
  {}

  //! Destructor.
  ~EssJournalNull() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Journalling functions.
  //@{
public:

  //! Do nothing.
  /*!
    \param time The current time.
    \param index The reaction index.
  */
  void
  recordReaction(const Number time, const int index) const
  {}

  //! Do nothing.
  /*!
    \param time The current time.
    \param timeToNextReaction The time increment to the next reaction.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordPopulations(const Number time, const Number timeToNextReaction,
		    const _Container& populations) const
  {}

  //! Do nothing.
  /*!
    \param time The current time.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordInitialPopulations(const Number time,
			   const _Container& populations) const
  {}

  //! Do nothing.
  /*!
    \param time The current time.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordFinalPopulations(const Number time,
			 const _Container& populations) const
  {}

  //@}
};



  
//============================================================================
//! Journal for exact stochastic simulations.  Record the reactions.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename T = double>
class EssJournalReactions :
  public EssJournalNull<T> {
  //
  // Private types.
  //
private:

  //! The base class.
  typedef EssJournalNull<T> Base;

  //
  // Public types.
  //
public:

  //! The number type.
  typedef typename Base::Number Number;

  //
  // Member data.
  //
private:

  //! The output stream.
  std::ostream* _out;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssJournalReactions();
  //! Copy constructor not implemented.
  EssJournalReactions(const EssJournalReactions&);
  //! Assignment operator not implemented.
  EssJournalReactions&
  operator=(const EssJournalReactions&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from an output stream.
  EssJournalReactions(std::ostream* out = 0) :
    Base(),
    _out(out)
  {}

  //! Destructor.
  ~EssJournalReactions() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Journalling functions.
  //@{
public:

  //! Record the reaction.
  /*!
    \param time The current time.
    \param index The reaction index.
  */
  void
  recordReaction(const Number time, const int index) {
    if (_out) {
      *_out << time << " " << index << "\n";
    }
  }

  //! Do nothing.
  using Base::recordPopulations;

  //! Do nothing.
  using Base::recordInitialPopulations;

  //! Do nothing.
  using Base::recordFinalPopulations;

  //@}
};

  
  

//============================================================================
//! Journal for exact stochastic simulations.  Record the final populations.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename T = double>
class EssJournalPopulationsFinal :
  public EssJournalNull<T> {
  //
  // Private types.
  //
private:

  //! The base class.
  typedef EssJournalNull<T> Base;

  //
  // Public types.
  //
public:

  //! The number type.
  typedef typename Base::Number Number;

  //
  // Member data.
  //
private:

  //! The output stream.
  std::ostream* _out;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssJournalPopulationsFinal();
  //! Copy constructor not implemented.
  EssJournalPopulationsFinal(const EssJournalPopulationsFinal&);
  //! Assignment operator not implemented.
  EssJournalPopulationsFinal&
  operator=(const EssJournalPopulationsFinal&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from an output stream.
  EssJournalPopulationsFinal(std::ostream* out = 0) :
    Base(),
    _out(out)
  {}

  //! Destructor.
  ~EssJournalPopulationsFinal() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Journalling functions.
  //@{
public:

  //! Do nothing.
  using Base::recordReaction;

  //! Do nothing.
  using Base::recordPopulations;

  //! Do nothing.
  using Base::recordInitialPopulations;

  //! Record the final populations.
  /*!
    \param time The current time.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordFinalPopulations(const Number time,
			 const _Container& populations) {
    if (_out) {
      *_out << time << "\n";
      writeArrayAscii(*_out, populations);
    }
  }

  //@}
};

  
  


//============================================================================
//! Journal for exact stochastic simulations.  Record the initial and final populations.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename T = double>
class EssJournalPopulationsInitialFinal :
  public EssJournalNull<T> {
  //
  // Private types.
  //
private:

  //! The base class.
  typedef EssJournalNull<T> Base;

  //
  // Public types.
  //
public:

  //! The number type.
  typedef typename Base::Number Number;

  //
  // Member data.
  //
private:

  //! The output stream.
  std::ostream* _out;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssJournalPopulationsInitialFinal();
  //! Copy constructor not implemented.
  EssJournalPopulationsInitialFinal(const EssJournalPopulationsInitialFinal&);
  //! Assignment operator not implemented.
  EssJournalPopulationsInitialFinal&
  operator=(const EssJournalPopulationsInitialFinal&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from an output stream.
  EssJournalPopulationsInitialFinal(std::ostream* out = 0) :
    Base(),
    _out(out)
  {}

  //! Destructor.
  ~EssJournalPopulationsInitialFinal() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Journalling functions.
  //@{
public:

  //! Do nothing.
  using Base::recordReaction;

  //! Do nothing.
  using Base::recordPopulations;

  //! Record the initial populations.
  /*!
    \param time The current time.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordInitialPopulations(const Number time,
			   const _Container& populations) {
    if (_out) {
      *_out << time << "\n";
      writeArrayAscii(*_out, populations);
    }
  }

  //! Record the final populations.
  /*!
    \param time The current time.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordFinalPopulations(const Number time,
			 const _Container& populations) {
    if (_out) {
      *_out << time << "\n";
      writeArrayAscii(*_out, populations);
    }
  }

  //@}
};

  
  


//============================================================================
//! Journal for exact stochastic simulations.  Record the populations at regular time intervals.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename T = double>
class EssJournalPopulationsTime :
  public EssJournalNull<T> {
  //
  // Private types.
  //
private:

  //! The base class.
  typedef EssJournalNull<T> Base;

  //
  // Public types.
  //
public:

  //! The number type.
  typedef typename Base::Number Number;

  //
  // Member data.
  //
private:

  //! The output stream.
  std::ostream* _out;
  //! The time.
  Number _time;
  //! The time interval between recording the populations.
  Number _timeInterval;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssJournalPopulationsTime();
  //! Copy constructor not implemented.
  EssJournalPopulationsTime(const EssJournalPopulationsTime&);
  //! Assignment operator not implemented.
  EssJournalPopulationsTime&
  operator=(const EssJournalPopulationsTime&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from an output stream.
  EssJournalPopulationsTime(std::ostream* out, const Number initialTime,
			    const Number timeInterval) :
    Base(),
    _out(out),
    _time(initialTime),
    _timeInterval(timeInterval)
  {}

  //! Destructor.
  ~EssJournalPopulationsTime() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Journalling functions.
  //@{
public:

  //! Do nothing.
  using Base::recordReaction;

  //! Record the populations if it is time.
  /*!
    \param time The current time.
    \param timeToNextReaction The time increment to the next reaction.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordPopulations(const Number time, const Number timeToNextReaction,
		    const _Container& populations) {
    if (_out) {
      while (time + timeToNextReaction > _time) {
	*_out << time << "\n";
	writeArrayAscii(*_out, populations);
	_time += _timeInterval;
      }
    }
  }

  //! Do nothing.
  using Base::recordInitialPopulations;

  //! Do nothing.
  using Base::recordFinalPopulations;

  //@}
};

  
  


//============================================================================
//! Journal for exact stochastic simulations.  Record the populations at regular reaction count intervals.
/*!
  \param T The number type.  The default is \c double .
*/
template<typename T = double>
class EssJournalPopulationsCount :
  public EssJournalNull<T> {
  //
  // Private types.
  //
private:

  //! The base class.
  typedef EssJournalNull<T> Base;

  //
  // Public types.
  //
public:

  //! The number type.
  typedef typename Base::Number Number;

  //
  // Member data.
  //
private:

  //! The output stream.
  std::ostream* _out;
  //! The count.
  int _count;
  //! The interval between recording the populations.
  int _interval;

  //
  // Not implemented.
  //
private:  

  //! Default constructor not implemented.
  EssJournalPopulationsCount();
  //! Copy constructor not implemented.
  EssJournalPopulationsCount(const EssJournalPopulationsCount&);
  //! Assignment operator not implemented.
  EssJournalPopulationsCount&
  operator=(const EssJournalPopulationsCount&);

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from an output stream.
  EssJournalPopulationsCount(std::ostream* out, const int interval) :
    Base(),
    _out(out),
    _count(interval),
    _interval(interval) {
    assert(interval > 0);
  }

  //! Destructor.
  ~EssJournalPopulationsCount() 
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Journalling functions.
  //@{
public:

  //! Do nothing.
  using Base::recordReaction;

  //! Record the populations if it is time.
  /*!
    \param time The current time.
    \param timeToNextReaction The time increment to the next reaction.
    \param populations The species populations.
  */
  template<class _Container>
  void
  recordPopulations(const Number time, const Number timeToNextReaction,
		    const _Container& populations) {
    if (_out) {
      if (--_count == 0) {
	*_out << time << "\n";
	writeArrayAscii(*_out, populations);
	_count = _interval;
      }
    }
  }

  //! Do nothing.
  using Base::recordInitialPopulations;

  //! Do nothing.
  using Base::recordFinalPopulations;

  //@}
};

  
  


END_NAMESPACE_STOCHASTIC

#define __stochastic_EssJournal_ipp__
#include "EssJournal.ipp"
#undef __stochastic_EssJournal_ipp__

#endif
