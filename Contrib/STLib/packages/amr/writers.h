// -*- C++ -*-

/*!
  \file amr/writers.h
  \brief writers that stores data.
*/

#if !defined(__amr_writers_h__)
#define __amr_writers_h__

#include "Orthtree.h"
#include "PatchDescriptor.h"

#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>

// If we are debugging the whole amr package.
#if defined(DEBUG_amr) && !defined(DEBUG_amr_writers)
#define DEBUG_amr_writers
#endif

BEGIN_NAMESPACE_AMR

//! Write the cell data in ParaView format.
/*!
  \relates Orthtree

  \param name The base name for the ParaView file and the VTK files.
  \param orthtree The orthtree.
  \param fieldDescriptors The vector of field descriptors.
 */
template<typename _Patch, class _Traits>
void
writeCellDataParaview(const std::string& name,
		      const Orthtree<_Patch, _Traits>& orthtree,
		      const PatchDescriptor<_Traits>& patchDescriptor);

END_NAMESPACE_AMR

#define __amr_writers_ipp__
#include "writers.ipp"
#undef __amr_writers_ipp__

#endif
