// -*- C++ -*-

/*! 
  \file rcb.h
  \brief Partitioning with recursive coordinate bisection.
*/

#if !defined(__concurrent_partition_rcb_h__)
#define __concurrent_partition_rcb_h__

#include "../defs.h"

#include "../../ads/array/Array.h"
#include "../../ads/functor/compose.h"
#include "../../ads/functor/Dereference.h"
#include "../../ads/functor/index.h"
#include "../../ads/iterator/IndirectIterator.h"

#include "../../geom/kernel/BBox.h"

// If we are debugging the whole concurrent package.
#if defined(DEBUG_concurrent) && !defined(DEBUG_concurrent_rcb)
#define DEBUG_concurrent_rcb
#endif

BEGIN_NAMESPACE_CONCURRENT

//! Determine an RCB partitioning of the records.
/*!
  \param num_processors is the number of processors.  (Input)
  \param num_records is the number of records.  (Input)
  \param identifiers is an array of length \c num_records that holds the 
  record identifiers.  (Input/Output)
  \param id_partition is an array of length \c (num_processors+1).  The 
  elements are pointers into the \c identifiers array that determine 
  the partition.  (Output)
  \param positions is an array of length \c N*num_records that holds the
  Cartesian positions of the records.  (Input)  In 3-D, the order of the 
  array is
  \code
  pos_0_x pos_0_y pos_0_z 
  pos_1_x pos_1_y pos_1_z 
  ...
  \endcode

  rcb() performs recursive coordinate bisection of data.  (The 
  data is composed of records which have identifiers and positions.)  
  It determines how to divide the records among the processors.  This is
  a sequential algorithm.  It can be run on any processor that has 
  the identifiers and the positions of all the records.

  Template parameters:  \c N is the dimension, it must be a positive integer.
  \c IDType is the identifier type.  For example, the identifier type 
  might be an integer that is the index of records in an array.  Or
  it could be a pointer to a record.  \c NumberType is the floating
  point number type.

  The arrays \c identifiers and \c id_partition determine the partition
  of the records.  \c id_partition is an array of pointers into the
  \c identifiers array.  The function sets the values to determine
  how many records each processor should have.  Processor \c n will 
  hold records with identifiers in the range [ \c id_partion[n] ..
  \c id_partition[n+1]).  The \c identifiers array will be permuted
  so the approprate record identifiers are in that range.  

  If you don't already have identifier arrays and position arrays, you 
  can allocate and free memory for the three arrays with rcb_allocate() 
  and rcb_deallocate().
 */
template<int N, typename IDType, typename NumberType>
void
rcb(const int num_processors, const int num_records,
     IDType* identifiers, IDType** id_partition,
     const NumberType* positions);

//! Allocate memory for an RCB calculation.
/*!
  \param num_processors is the number of processors.
  \param num_records is the number of records.
  \param identifiers will point to an array of length \c num_records.
  \param id_partition will point to an array of length \c (num_processors+1).
  \param positions will point to an array of length \c N*num_records.
  
  It is not necessary to use this function; you can allocate and manage the 
  arrays yourself.  It is provided only for convenience.  

  You can free the arrays allocated here with rcb_deallocate().
 */
template<int N, typename IDType, typename NumberType>
void
rcb_allocate(const int num_processors, const int num_records,
	      IDType*& identifiers, IDType**& id_partition, 
	      NumberType*& positions);

//! Free the memory allocated in rcb_allocate().
template<typename IDType, typename NumberType>
void
rcb_deallocate(IDType*& identifiers, IDType**& id_partition, 
		NumberType*& positions);

END_NAMESPACE_CONCURRENT

#define __partition_rcb_ipp__
#include "rcb.ipp"
#undef __partition_rcb_ipp__

#endif
