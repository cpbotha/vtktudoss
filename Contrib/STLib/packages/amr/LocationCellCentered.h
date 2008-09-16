// -*- C++ -*-

/*!
  \file amr/LocationCellCentered.h
  \brief Computes the location of cell centers in arrays.
*/

#if !defined(__amr_LocationCellCentered_h__)
#define __amr_LocationCellCentered_h__

#include "defs.h"

#include "../array/Array.h"

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_LocationCellCentered)
#define DEBUG_amr_LocationCellCentered
#endif

BEGIN_NAMESPACE_AMR

//! Computes the location of cell centers in arrays.
/*!
  \param _Traits Traits for the orthtree.
*/
template<typename _Traits>
class
LocationCellCentered {
  //
  // Public types.
  //
public:

  //! A Cartesian position.
  typedef array::Array<typename _Traits::Number, _Traits::Dimension> Point;

  //
  // Member data.
  //
private:

  Point _offsets, _origin;

  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{
public:

  //! Construct from the Cartesian lower corner, the Cartesian extents and the array extents.
  template<typename _ExtentList>
  LocationCellCentered(const Point& lowerCorner, const Point& cartesianExtents,
		       const _ExtentList& arrayExtents) :
    _offsets(cartesianExtents / arrayExtents),
    _origin(lowerCorner + 0.5 * _offsets) {
  }

  //! Copy constructor.
  LocationCellCentered(const LocationCellCentered& other) :
    _offsets(other._offsets),
    _origin(other._origin) {
  }

  //! Assignment operator.
  LocationCellCentered&
  operator=(const LocationCellCentered& other) {
    if (this != &other) {
      _offsets = other._offsets;
      _origin = other._origin;
    }
    return *this;
  }

  //! Destructor.
  ~LocationCellCentered() {
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Functor.
  //@{
public:

  //! Compute the location for the specified cell center.
  template<typename _IndexList>
  void
  operator()(const _IndexList& index, Point* location) const {
    for (std::size_t n = 0; n != _Traits::Dimension; ++n) {
      (*location)[n] = _origin[n] + index[n] * _offsets[n];
    }
  }

  //! Compute the location for the specified cell center.
  template<typename _IndexList>
  Point
  operator()(const _IndexList& index) const {
    Point location;
    for (std::size_t n = 0; n != _Traits::Dimension; ++n) {
      location[n] = _origin[n] + index[n] * _offsets[n];
    }
    return location;
  }

  //@}
  //--------------------------------------------------------------------------
  //! \name Equality.
  //@{
public:

  //! Return true if the data structures are equal.
  bool
  operator==(const LocationCellCentered& other) {
    return _offsets == other._offsets && _origin == other._origin;
  }

  //@}
};

END_NAMESPACE_AMR

#endif
