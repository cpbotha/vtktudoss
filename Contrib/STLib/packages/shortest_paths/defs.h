// -*- C++ -*-

/*! 
  \file shortest_paths/defs.h
  \brief Definitions for the shortest-paths package.
*/

#if !defined(__shortest_paths_defs_h__)
//! Include guard.
#define __shortest_paths_defs_h__

// If we are debugging everything in STLib.
#if defined(DEBUG_stlib) && !defined(DEBUG_shortest_paths)
#define DEBUG_shortest_paths
#endif

//! Begin the shortest_paths namespace.
#define BEGIN_NAMESPACE_SHORTEST_PATHS namespace shortest_paths {
//! End the shortest_paths namespace.
#define END_NAMESPACE_SHORTEST_PATHS }

#endif
