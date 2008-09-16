// -*- C++ -*-

/*!
  \file amr/MessageStreamBuffer.h
  \brief Message stream buffer.
*/

#if !defined(__amr_MessageStreamBuffer_h__)
#define __amr_MessageStreamBuffer_h__

#include "defs.h"

#include <cassert>
#include <cstring>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_MessageStreamBuffer)
#define DEBUG_amr_MessageStreamBuffer
#endif

BEGIN_NAMESPACE_AMR

//! Message stream buffer.
/*!
*/
class
MessageStreamBuffer {
  //
  // Public types.
  //
public:

  //! The size type.
  typedef int SizeType;

  //
  // Member data.
  //
protected:

  SizeType _size;
  SizeType _capacity;
  char* _data;

  //
  // Not implemented.
  //
private:

  MessageStreamBuffer();

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
protected:

  //! Construct and reserve memory.
  MessageStreamBuffer(const SizeType capacity) :
    _size(0),
    _capacity(capacity),
    _data() {
    assert(capacity > 0);
    _data = new char[_capacity];
  }

  //! Copy constructor.
  MessageStreamBuffer(const MessageStreamBuffer& other) :
    _size(other._size),
    _capacity(other._capacity),
    _data() {
    _data = new char[_capacity];
    memcpy(_data, other._data, _size);
  }

  //! Assignment operator.
  MessageStreamBuffer&
  operator=(const MessageStreamBuffer& other) {
    if (this != &other) {
      clear();
      if (_capacity < other._size) {
	reserve(other._size);
      }
      _size = other._size;
      memcpy(_data, other._data, _size);
    }
    return *this;
  }

  //! Destructor.
  ~MessageStreamBuffer() {
    delete[] _data;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Return the size.
  SizeType
  getSize() const {
    return _size;
  }

  //! Return the capacity.
  SizeType
  getCapacity() const {
    return _capacity;
  }

  //! Return a const pointer to the data.
  const char*
  getData() const {
    return _data;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
protected:

  //! Return true if buffer contents are equal.
  bool
  operator==(const MessageStreamBuffer& other) const {
    return _size == other._size && (memcmp(_data, other._data, _size) == 0);
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Return a pointer to the data.
  char*
  getData() {
    return _data;
  }

  //! Resize the buffer.
  void
  resize(const int size) {
    clear();
    reserve(size);
    _size = size;
  }

  //! Clear the buffer.
  void
  clear() {
    _size = 0;
  }

  //! Reserve at least the indicated capacity.
  /*!
    The contents of the buffer are preserved.
  */
  void
  reserve(const SizeType newCapacity) {
#ifdef DEBUG_amr_MessageStreamBuffer
    assert(newCapacity > 0);
#endif
    if (_capacity < newCapacity) {
      // Increace the capacity by factors of 2.
      while (_capacity < newCapacity) {
	_capacity *= 2;
      }
      // Allocate new memory.
      char* _tmp = new char[_capacity];
      // Copy into the new memory.
      memcpy(_tmp, _data, _size);
      // Free the old memory.
      delete[] _data;
      // Point to the new memory.
      _data = _tmp;
    }
  }

  //@}
};

END_NAMESPACE_AMR

#define __amr_MessageStreamBuffer_ipp__
#include "MessageStreamBuffer.ipp"
#undef __amr_MessageStreamBuffer_ipp__

#endif
