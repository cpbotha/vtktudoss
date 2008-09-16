// -*- C++ -*-

/*!
  \file amr/Element.h
  \brief Element that stores data.
*/

#if !defined(__amr_Element_h__)
#define __amr_Element_h__

#include "defs.h"

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_Element)
#define DEBUG_amr_Element
#endif

BEGIN_NAMESPACE_AMR

//! Element that stores data.
/*!
  \param _Data The data held in each node.
*/
template<typename _Data>
class
Element {
  //
  // Public types.
  //
public:

  //! The data type.
  typedef _Data Data;

  //
  // More public types.
  //
public:

  //
  // Member data.
  //
private:

  Data _data;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Make a default node.
  Element() :
    _data() {
  }

  //! Copy constructor.
  Element(const Element& other) :
    _data(other._data) {
  }

  //! Assignment operator.
  Element&
  operator=(const Element& other) {
    if (this != &other) {
      _data = other._data;
    }
    return *this;
  }

  //! Destructor.
  ~Element() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{
public:

  //! Return a const reference to the data.
  const Data&
  getData() const {
    return _data;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
public:

  //! Return true if the data structures are equal.
  bool
  operator==(const Element& other) {
    return _data == other._data;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{
public:

  //! Return a reference to the data.
  Data&
  getData() {
    return _data;
  }

  //! Link to the other nodes in the tree.
  /*!
    \param spatialIndex The spatial index for this element.

    This node has no adjacent links, so nothing is done.

    The spatial index is a template parameter so that this class does not 
    need to know about them.
  */
  template<typename _SpatialIndex>
  void
  link(const _SpatialIndex& spatialIndex) {
  }

  //! Unlink from the other nodes in the tree.
  /*!
    This node has no adjacent links, so nothing is done.
  */
  void
  unlink() {
  }

  //@}
};

END_NAMESPACE_AMR

#define __amr_Element_ipp__
#include "Element.ipp"
#undef __amr_Element_ipp__

#endif
