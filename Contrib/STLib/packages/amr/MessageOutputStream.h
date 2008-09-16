// -*- C++ -*-

/*!
  \file amr/MessageOutputStream.h
  \brief Message output stream.
*/

#if !defined(__amr_MessageOutputStream_h__)
#define __amr_MessageOutputStream_h__

#include "MessageStreamBuffer.h"

#include <cassert>
#include <cstring>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_MessageOutputStream)
#define DEBUG_amr_MessageOutputStream
#endif

BEGIN_NAMESPACE_AMR

//! Message output stream.
/*!
*/
class
MessageOutputStream : public MessageStreamBuffer {
  //
  // Private types.
  //
private:

  typedef MessageStreamBuffer Base;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct and reserve memory.
  MessageOutputStream(const SizeType capacity = 1) :
    Base(capacity) {
  }

  //! Copy constructor.
  MessageOutputStream(const MessageOutputStream& other) :
    Base(other) {
  }

  //! Assignment operator.
  MessageOutputStream&
  operator=(const MessageOutputStream& other) {
    if (this != &other) {
      Base::operator=(other);
    }
    return *this;
  }

  //! Destructor.
  ~MessageOutputStream() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Return the remaining capacity.
  SizeType
  getRemainingCapacity() const {
    return getCapacity() - getSize();
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
public:

  //! Return true if buffer contents are equal.
  bool
  operator==(const MessageOutputStream& other) const {
    return Base::operator==(other);
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Reserve at least the indicated additional capacity.
  void
  reserveAdditional(const SizeType size) {
    reserve(getSize() + size);
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Output.
  //@{
public:

  //! Output a boolean.
  MessageOutputStream&
  operator<<(const bool x) {
    write(x);
    return *this;
  }

  //! Output a character.
  MessageOutputStream&
  operator<<(const char x) {
    write(x);
    return *this;
  }

  //! Output an unsigned character.
  MessageOutputStream&
  operator<<(const unsigned char x) {
    write(x);
    return *this;
  }

  //! Output a short integer.
  MessageOutputStream&
  operator<<(const short x) {
    write(x);
    return *this;
  }

  //! Output an unsigned short integer.
  MessageOutputStream&
  operator<<(const unsigned short x) {
    write(x);
    return *this;
  }

  //! Output an integer.
  MessageOutputStream&
  operator<<(const int x) {
    write(x);
    return *this;
  }

  //! Output an unsigned integer. 
  MessageOutputStream&
  operator<<(const unsigned x) {
    write(x);
    return *this;
  }

  //! Output a long integer.
  MessageOutputStream&
  operator<<(const long x) {
    write(x);
    return *this;
  }

  //! Output an unsigned long integer.
  MessageOutputStream&
  operator<<(const unsigned long x) {
    write(x);
    return *this;
  }

  //! Output a float.
  MessageOutputStream&
  operator<<(const float x) {
    write(x);
    return *this;
  }

  //! Output a double.
  MessageOutputStream&
  operator<<(const double x) {
    write(x);
    return *this;
  }

  //! Output a built-in type or plain-old data.
  template<typename _T>
  void
  write(const _T x) {
#ifdef DEBUG_amr_MessageOutputStream
    assert(int(_size + sizeof(_T)) <= _capacity);
#endif
    memcpy(_data + _size, &x, sizeof(_T));
    _size += sizeof(_T);
  }

  //! Output an array of variables.
  template<typename _T>
  void
  write(const _T* x, const int n) {
#ifdef DEBUG_amr_MessageOutputStream
    assert(int(_size + n * sizeof(_T)) <= _capacity);
#endif
    memcpy(_data + _size, x, n * sizeof(_T));
    _size += n * sizeof(_T);
  }

  //@}
};

//! Write output to the stream.
/*!
  \relates MessageOutputStream
*/
template<typename _T>
inline
MessageOutputStream&
operator<<(MessageOutputStream& out, const _T& x) {
  x.write(out);
  return out;
}

END_NAMESPACE_AMR

#define __amr_MessageOutputStream_ipp__
#include "MessageOutputStream.ipp"
#undef __amr_MessageOutputStream_ipp__

#endif
