// -*- C++ -*-

/*!
  \file amr/FieldDescriptor.h
  \brief Describe a data field.
*/

#if !defined(__amr_FieldDescriptor_h__)
#define __amr_FieldDescriptor_h__

#include "defs.h"

#include <string>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_FieldDescriptor)
#define DEBUG_amr_FieldDescriptor
#endif

BEGIN_NAMESPACE_AMR

//! Describe a data field.
class
FieldDescriptor {

  //
  // Member data.
  //
private:

  int _numberOfComponents;
  std::string _name;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the number of components and the name.
  FieldDescriptor(const int numberOfComponents, std::string name) :
    _numberOfComponents(numberOfComponents),
    _name(name) {
  }

  //! Copy constructor.
  FieldDescriptor(const FieldDescriptor& other) :
    _numberOfComponents(other._numberOfComponents),
    _name(other._name) {
  }

  //! Assignment operator.
  FieldDescriptor&
  operator=(const FieldDescriptor& other) {
    if (this != &other) {
      _numberOfComponents = other._numberOfComponents;
      _name = other._name;
    }
    return *this;
  }

  //! Destructor.
  ~FieldDescriptor() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Get the number of components.
  int
  getNumberOfComponents() const {
    return _numberOfComponents;
  }

  //! Get the name.
  const std::string&
  getName() const {
    return _name;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
public:

  //! Return true if the data structures are equal.
  bool
  operator==(const FieldDescriptor& other) {
    return _numberOfComponents == other._numberOfComponents &&
      _name == other._name;
  }

  //@}
};

END_NAMESPACE_AMR

#endif
