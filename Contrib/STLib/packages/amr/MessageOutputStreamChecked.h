// -*- C++ -*-

/*!
  \file amr/MessageOutputStreamChecked.h
  \brief Message output stream.
*/

#if !defined(__amr_MessageOutputStreamChecked_h__)
#define __amr_MessageOutputStreamChecked_h__

#include "MessageOutputStream.h"

#include <cassert>
#include <cstring>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_MessageOutputStreamChecked)
#define DEBUG_amr_MessageOutputStreamChecked
#endif

BEGIN_NAMESPACE_AMR

//! Message output stream.
/*!
*/
class
MessageOutputStreamChecked : public MessageOutputStream {
  //
  // Private types.
  //
private:

  typedef MessageOutputStream Base;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct and reserve memory.
  MessageOutputStreamChecked(const SizeType capacity = 1) :
    Base(capacity) {
  }

  //! Copy constructor.
  MessageOutputStreamChecked(const MessageOutputStreamChecked& other) :
    Base(other) {
  }

  //! Assignment operator.
  MessageOutputStreamChecked&
  operator=(const MessageOutputStreamChecked& other) {
    if (this != &other) {
      Base::operator=(other);
    }
    return *this;
  }

  //! Destructor.
  ~MessageOutputStreamChecked() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Output.
  //@{
public:

  //! Output a boolean.
  MessageOutputStreamChecked&
  operator<<(const bool x) {
    write(x);
    return *this;
  }

  //! Output a character.
  MessageOutputStreamChecked&
  operator<<(const char x) {
    write(x);
    return *this;
  }

  //! Output an unsigned character.
  MessageOutputStreamChecked&
  operator<<(const unsigned char x) {
    write(x);
    return *this;
  }

  //! Output a short integer.
  MessageOutputStreamChecked&
  operator<<(const short x) {
    write(x);
    return *this;
  }

  //! Output an unsigned short integer.
  MessageOutputStreamChecked&
  operator<<(const unsigned short x) {
    write(x);
    return *this;
  }

  //! Output an integer.
  MessageOutputStreamChecked&
  operator<<(const int x) {
    write(x);
    return *this;
  }

  //! Output an unsigned integer. 
  MessageOutputStreamChecked&
  operator<<(const unsigned x) {
    write(x);
    return *this;
  }

  //! Output a long integer.
  MessageOutputStreamChecked&
  operator<<(const long x) {
    write(x);
    return *this;
  }

  //! Output an unsigned long integer.
  MessageOutputStreamChecked&
  operator<<(const unsigned long x) {
    write(x);
    return *this;
  }

  //! Output a float.
  MessageOutputStreamChecked&
  operator<<(const float x) {
    write(x);
    return *this;
  }

  //! Output a double.
  MessageOutputStreamChecked&
  operator<<(const double x) {
    write(x);
    return *this;
  }

  //! Output a built-in type or plain-old data.
  template<typename _T>
  void
  write(const _T x) {
    reserveAdditional(sizeof(_T));
    Base::write(x);
  }

  //! Output an array of variables.
  template<typename _T>
  void
  write(const _T* x, const int n) {
    reserveAdditional(n * sizeof(_T));
    Base::write(x, n);
  }

  //@}
};

//! Write output to the stream.
/*!
  \relates MessageOutputStreamChecked
*/
template<typename _T>
inline
MessageOutputStreamChecked&
operator<<(MessageOutputStreamChecked& out, const _T& x) {
  x.write(out);
  return out;
}

END_NAMESPACE_AMR

#define __amr_MessageOutputStreamChecked_ipp__
#include "MessageOutputStreamChecked.ipp"
#undef __amr_MessageOutputStreamChecked_ipp__

#endif
