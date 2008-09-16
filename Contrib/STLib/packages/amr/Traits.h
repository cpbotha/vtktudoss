// -*- C++ -*-

/*!
  \file amr/Traits.h
  \brief The traits for an orthtree.
*/

#if !defined(__amr_Traits_h__)
#define __amr_Traits_h__

#include "SpatialIndexMorton.h"

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_Traits)
#define DEBUG_amr_Traits
#endif

BEGIN_NAMESPACE_AMR

//! The traits for an orthtree.
/*!
  \param _Dimension The dimension of the space.
  \param _MaximumLevel The maximum level in the tree.
  \param _SpatialIndex The spatial index data structure.  This determines
  the level and position of the node.  It also holds a key for storing the
  node in the \c std::map data structure.
  \param _Number The real number type.
*/
template<std::size_t _Dimension, std::size_t _MaximumLevel, 
	 template<std::size_t, std::size_t> class _SpatialIndex = 
	 SpatialIndexMorton,
	 typename _Number = double>
class
Traits {
  //
  // Types and enumerations.
  //
public:

  //! The space dimension and the maximum level.
  enum {Dimension = _Dimension, MaximumLevel = _MaximumLevel};
  //! A multi-index in a multi-array.
  typedef array::Array<int, Dimension> IndexList;
  //! A list if sizes.
  typedef array::Array<std::size_t, Dimension> SizeList;
  //! The spatial index.
  typedef _SpatialIndex<Dimension, MaximumLevel> SpatialIndex;
  //! The number type.
  typedef _Number Number;
  //! A Cartesian point.
  typedef array::Array<Number, Dimension> Point;
  //! The number of orthants = 2^Dimension.
  enum {NumberOfOrthants = SpatialIndex::NumberOfOrthants};

  //
  // Not implemented.
  //
private:

  //! The default constructor is not implemented.
  Traits();
};

END_NAMESPACE_AMR

#endif
