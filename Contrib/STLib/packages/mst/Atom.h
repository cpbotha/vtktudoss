// -*- C++ -*-

/*! 
  \file Atom.h
  \brief An atom is represented by a position and a radius.
*/

#if !defined(__mst_Atom_h__)
#define __mst_Atom_h__

#include "defs.h"

#include "../geom/kernel/Ball.h"
#include "../geom/kernel/content.h"

BEGIN_NAMESPACE_MST


//! An atom.
template<typename T>
class Atom : 
  public geom::Ball<3,T> {

  //
  // Private types.
  //

private:

  typedef geom::Ball<3,T> Base;

  //
  // Public types.
  //

public:

  //! The number type.
  typedef typename Base::Number Number;
  //! A Cartesian point.
  typedef typename Base::Point Point;

public:
    
  //--------------------------------------------------------------------------
  //! \name Constructors etc.
  //@{

  //! Default constructor.  Uninitialized memory.
  Atom() :
    Base()
  {}

  //! Copy constructor.
  Atom(const Atom& other) :
    Base(other)
  {}

  //! Assignment operator.
  Atom&
  operator=(const Atom& other) {
    // Avoid assignment to self.
    if (&other != this) {
      Base::operator=(other);
    }
    // Return *this so assignments can chain.
    return *this;
  }

  //! Trivial destructor.
  ~Atom()
  {}

  //! Construct from the center and radius.
  Atom(const Point& center, const Number radius) :
    Base(center, radius)
  {}

  //@}
  //--------------------------------------------------------------------------
  //! \name Accessors.
  //@{

  //! Return the center.
  using Base::getCenter;

  //! Return the radius.
  using Base::getRadius;

  //@}
  //--------------------------------------------------------------------------
  //! \name Manipulators.
  //@{

  //! Set the center.
  using Base::setCenter;

  //! Set the radius.
  using Base::setRadius;

  //@}
  //--------------------------------------------------------------------------
  //! \name Mathematical functions.
  //@{
  
  //! Return true if the point is inside this atom.
  using Base::isInside;
  
  //@}
};


//! Return true if b clips a.  c is the distance between the centers.
template<typename T>
inline
bool
doesClip(const Atom<T>& a, const Atom<T>& b, const T c) {
  // If the two balls do not intersect.
  if (a.getRadius() + b.getRadius() <= c) {
    return false;
  }
  // If the first contains the second.
  if (a.getRadius() - b.getRadius() >= c) {
    return false;
  }
  return true;
}


//! Return true if b clips a.
template<typename T>
inline
bool
doesClip(const Atom<T>& a, const Atom<T>& b) {
  const T c = geom::computeDistance(a.getCenter(), b.getCenter());
  return doesClip(a, b, c);
}


//! Return the clipping plane distance.
template<typename T>
inline
T
computeClippingPlaneDistance(const Atom<T>& a, const Atom<T>& b) {
  // The distance between the atom's centers.
  const T c = geom::computeDistance(a.getCenter(), b.getCenter());
#ifdef DEBUG_mst
  assert(doesClip(a, b, c));
#endif
  // The atoms should not have the same centers.
  assert(c != 0);
  
  // If the second atom contains the first.
  if (b.getRadius() - a.getRadius() >= c) {
    // The whole atom is clipped.  Return -infinity.
    return - std::numeric_limits<T>::max();
  }
  // Otherwise, the spheres intersect.

  // Compute the clipping plane distance.
  const T dist = (a.getRadius() * a.getRadius() - b.getRadius() * b.getRadius()
		  + c * c) / (2 * c);

  // Check that the plane clips the atom.
  const T bound = ((1.0 + 10.0 * std::numeric_limits<T>::epsilon()) * 
		   a.getRadius());
  assert(-bound <= dist && dist <= bound);

  // Return the clipping plane distance.
  return dist;
}


END_NAMESPACE_MST

#endif
