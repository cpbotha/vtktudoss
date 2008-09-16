// -*- C++ -*-

/*!
  \file array/MultiArrayTypes.h
  \brief Array types.
*/

#if !defined(__array_MultiArrayTypes_h__)
#define __array_MultiArrayTypes_h__

#include "IndexTypes.h"

#include <boost/call_traits.hpp>

#include <iterator>

// If we are debugging the whole array package.
#if defined(DEBUG_array) && !defined(DEBUG_array_MultiArrayTypes)
#define DEBUG_array_MultiArrayTypes
#endif

BEGIN_NAMESPACE_ARRAY

//! Array types.
template<typename _T, std::size_t _Dimension>
class
MultiArrayTypes : public IndexTypes<_Dimension> {
  //
  // Enumerations.
  //
public:

  //! The number of dimensions.
  enum {Dimension = _Dimension};

  //
  // Public types.
  //
public:

  // Types for STL compliance.

  //! The element type of the array.  
  typedef _T value_type;

  //! A pointer to an array element.
  typedef value_type* pointer;
  //! A pointer to a constant array element.
  typedef const value_type* const_pointer;

  //! An iterator in the array.
  typedef value_type* iterator;
  //! A iterator on constant elements in the array.
  typedef const value_type* const_iterator;

  //! A reverse iterator in the array.
  typedef std::reverse_iterator<iterator> reverse_iterator;
  //! A reverse iterator on constant elements in the array.
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  //! A reference to an array element.
  typedef value_type& reference;
  //! A reference to a constant array element.
  typedef const value_type& const_reference;

  // Other types.

  //! The parameter type.
  /*! This is used for passing the value type as an argument. */
  typedef typename boost::call_traits<value_type>::param_type Parameter;
};

END_NAMESPACE_ARRAY

#endif
