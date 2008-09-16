// -*- C++ -*-

#if !defined(__amr_h__)
#define __amr_h__

#include "amr/Orthtree.h"
#include "amr/Traits.h"

BEGIN_NAMESPACE_AMR

/*! 
  \file amr.h
  \brief Includes the adaptive mesh refinement classes.
*/

//=============================================================================
//=============================================================================
/*!
  \mainpage Adaptive Mesh Refinement

  This package is a C++ templated class library designed to support adaptive
  mesh refinement (AMR) for 
  logically Cartesian, structured meshes. A heirarchy of sub-grids are stored
  in a tree data structure. Each node of the tree is a structured mesh with
  the same index ranges. In N-dimensional space, refinement replaces a single
  mesh with 2<sup>N</sup> meshes; coarsening does the opposite. 

  Notable features include:
  - Support for any space dimension. (The dimension is a template parameter.)
  - Mixed-mode concurrency (both shared-memory and distributed-memory) to
  effectively utilize a wide range of architectures, including clusters
  of multi-core processors.

  This package has a <em>linear orthree</em> called amr::Orthtree .
  Each leaf has an associated <em>spatial index</em> (implemented in 
  amr::SpatialIndexMorton) that encodes its level and position in the tree.
  For distributed-memory applications, use amr::DistributedOrthtree .

  Use this sub-package by including the desired class or by including
  the file amr.h.

  Follow the links below for documentation and results.
  - \ref amr_introduction
  - \ref amr_considerations
  - \ref amr_orthant
  - \ref amr_index
  - \ref amr_container
  - \ref amr_concurrency
  - \ref amr_performance
  - \ref amr_compiling
  - \ref amr_links

  <!-- - \ref amr_examples-->
*/


//=============================================================================
//=============================================================================
/*!
\page amr_introduction Introduction

  A tree data structure that recursively divides a 2-D region into quadrants
  is called a \e quadtree.  (http://en.wikipedia.org/wiki/Quadtree)
  The 3-D equivalent is 
  an \e octree.  The N-D analogue of the 2-D quadrant and the 3-D octant 
  is an \e orthant.  (http://en.wikipedia.org/wiki/Orthant) 
  In N-D space there are \f$2^N\f$ orthants.  A tree that
  recursively divides space into orthants has been called an 
  <em>N-D quadtree</em> and a <em>hyper-octree</em>.  <!--CONTINUE cite-->
  We think the term \e orthtree is a descriptive and succinct name for such
  data structures.  As we have implemented a data structure that is generic
  in the space dimension (the dimension is a template parameter) we will 
  use the term orthtree instead of quadtree or octree.

  There are two primary methods of representing an orthtree: 
  <em>pointer-based</em> and \e linear. With a pointer-based representation,
  one stores the branches as well as the leaves of the tree.  
  Each branch stores \f$2^N\f$ pointers to nodes in the tree.  
  Here we use the term \e node to refer to either a branch or a leaf.
  These pointers the to \e children of each branch are necessary for 
  searching for a leaf.  Usually, each node also stores a pointer to its
  \e parent.  This is convenient for many algorithms, including traversing 
  the nodes.  The position of a leaf is implicitly stored in the branches
  leading from the root to the leaf.

  With the linear representation of an orthtree, only the leaves are stored.
  Instead of implicitly representing the positions of the leaves with branches,
  <em>spatial indices</em> are used to explicitly describe these positions.
  In this context, we call the pair of a spatial index and the leaf a \e node.
  "Linear" refers to the fact that the nodes may be stored in a sequence 
  container (an array, for example).  A specific node may be accessed 
  by searching for its spatial index.
*/


//=============================================================================
//=============================================================================
/*!
\page amr_considerations Considerations

There are a few basic design issues to consider when using AMR. First consider
how one refines a region. There are several possibilities:
- Finer grids are formed from the orthants of a single coarser grid. 
This leads to a regular tree data structure.
- Finer grids are arbitrary rectilinear portions of a single coarser grid.
This leads to an irregular tree data structure.
- Finer grids may be arbitrary rectilinear portions of one or more coarser
grids. In this case the hierarchy of grids could be represented with a graph.
.
We will name these three refinement schemes \e orthant, \e nested, and 
\e arbitrary.

*/

//=============================================================================
//=============================================================================
/*!
\page amr_orthant Orthant Arithmetic


\section amr_orthant_numbering Orthant Numbering

The \f$2^N\f$ orthants are numbered from 0 to \f$2^N - 1\f$.  The following
figure and table show this numbering for the 2-D case of quadrants.

\image html quadrant.jpg "The quadrants."

<table>
<tr> <th> Number <th> y <th> x <th> Binary Number
<tr> <td> 0 <td> - <td> - <td> 00
<tr> <td> 1 <td> - <td> + <td> 01
<tr> <td> 2 <td> + <td> - <td> 10
<tr> <td> 3 <td> + <td> + <td> 11
</table>

Note that the origin is in the lower left corner.  In computer graphics, one 
often places the origin in the \e upper left corner so that row numbering
starts at the top of the screen.  We instead use the convention that is 
more familiar in geometry and physics.

The table below lists octant numbering.

<table>
<tr> <th> Number <th> z <th> y <th> x <th> Binary Number
<tr> <td> 0 <td> - <td> - <td> - <td> 000
<tr> <td> 1 <td> - <td> - <td> + <td> 001
<tr> <td> 2 <td> - <td> + <td> - <td> 010
<tr> <td> 3 <td> - <td> + <td> + <td> 011
<tr> <td> 4 <td> + <td> - <td> - <td> 100
<tr> <td> 5 <td> + <td> - <td> + <td> 101
<tr> <td> 6 <td> + <td> + <td> - <td> 110
<tr> <td> 7 <td> + <td> + <td> + <td> 111
</table>

This is a nifty numbering system because the binary representation of 
the orthant encodes the \e N coordinate directions (0 for negative and 1 for
positive).  Specifically, the \f$n^{\mathrm{th}}\f$ bit (starting with the
least significant bit) encodes the direction for the \f$n^{\mathrm{th}}\f$
coordinate.  Let the coordinates be labeled from 0 to \e N - 1.
If \e i is the orthant number, the we can get the direction bit (0 or 1) for
the \f$n^{\mathrm{th}}\f$ coordinate with a shift and a modulus.
\code
(i >> n) % 2 \endcode
It easy to transform this to a direction sign (\f$\pm 1\f$).
\code
((i >> n) % 2) * 2 - 1 \endcode
By the way, you can write express multiplication and division by a power
of 2 as left or right shifting.  <code>x * 2</code> is the same as
<code>x << 1</code>. Any reasonable compiler will generate the same 
(efficient) code regardless of whether you use the former or latter method.


\section amr_orthant_directions Directions

In N-D space there are \f$2 N\f$ signed coordinate directions.  We label 
these with integers from 0 to \f$2 N - 1\f$.  In the figure and table below
we show the direction numbering in 2-D.

\image html directions.jpg "The direction numbering in 2-D."

<table>
<tr> <th> Number <th> Direction
<tr> <td> 0 <td> -x
<tr> <td> 1 <td> +x
<tr> <td> 2 <td> -y
<tr> <td> 3 <td> +y
</table>

The next table gives the direction numbering in 3-D.
<table>
<tr> <th> Number <th> Direction
<tr> <td> 0 <td> -x
<tr> <td> 1 <td> +x
<tr> <td> 2 <td> -y
<tr> <td> 3 <td> +y
<tr> <td> 4 <td> -z
<tr> <td> 5 <td> +z
</table>

This is a nifty numbering scheme because is easy to extract the coordinate
number and the coordinate direction.  For direction <code>i</code>
the former is <code>i / 2</code> and the latter is <code>i % 2</code>.
(Here we use 0 and 1 to indicate negative and positive.)
*/

//=============================================================================
//=============================================================================
/*!
\page amr_index The Spatial Index

The geometry of a node in an orthtree may be described by its level and 
the \e N index coordinates of its lower corner.  There are various ways of 
storing this information.  One can use an integer to store the level and
\e N integers to store the index coordinates.  Alternatively one can pack
this information into integer or bit array.  For the packed case, one 
can store each index coordinate in contiguous bits, or one can interleave
the bits.  Typically, the method of storage directly translates into 
an ordering of the nodes.  If one packs the information into a single 
integer type, then one uses the standard less than comparison on 
integers to order the nodes.  

There is a standard spatial indexing scheme called the \e quadcode which
stores the 
level and coordinates in a single integer type.  The index coordinates
are interleaved and packed into the most significant bits.  The level is
stored in the less significant bits.  This translates into an ordering of 
the nodes that has some handy properties.  Most importantly, there is 
a spatial coherence to a block of nodes accessed in order.  This is useful
both from an algorithmic perspective (some searching operations are easier)
and for performance (fewer cache misses).  In addition, this ordering can
be used for partitioning the nodes across processors.

\image html zOrdering.jpg "Z-ordering."

While the quadcode is handy for ordering the nodes, a non-interleaved
representation is more useful for manipulating the spatial index.  When
determining the parent, children, or neighbors of a node, one needs to
manipulate the index coordinates and level.  With the quadcode, one extracts
the coordinates and level, manipulates then, and then packs them back
into the quadcode.  Because the bits are interleaved, this unpacking and 
packing is fairly expensive.  (The computational complexity of both packing
and unpacking is linear in the depth of the tree.)  If one uses a 
separated or non-interleaved representation for the spatial index, then
manipulations can be done in constant time.  

For our implementation, we store the interleaved quadcode in a single
integer type and also separately store the index coordinates and level.
We use the quadcode for ordering the nodes and use the separate 
index coordinates and level for spatial index arithmetic.  This makes
the algorithms easier to implement and more efficient, but
roughly doubles the storage requirements for the spatial index.
If one is storing a significant amount of data in each node (a small 
grid, or a spectral element), then this extra storage requirement for the
spatial indices will be neglible.
*/

//=============================================================================
//=============================================================================
/*!
\page amr_container The Node Container

Each node in the orthtree has a spatial index and data.  For a static tree,
one can store the nodes in an array, sorted by the spatial indices.  One
can then find a node in logarithmic time with a binary search.  For a dynamic
tree, storing the nodes in an array may not efficient because inserting or 
deleting a node has linear complexity in the size of the array.  One can
mitigate this by caching the insertions and deletions.  Another drawback
to storing the nodes in an array is that insertions and deletions invalidate
pointers.  If one were storing iterators to adjacent neighbors, this 
information would need to be recomputed following an insertion or deletion.
A different approach (and the one that we have taken) is to store the nodes
in a paired associative container. <!--CONTINUE cite--> We use 
\c std::map .  The node is a \c std::pair with \c first being
a const spatial index and \c second being the data.  A node can still be
found in logarithmic time.  Inserting new nodes does not invalidate any
iterators.  Deleting a node only invalidates iterators to the deleted
node.  Finally, when iterating through the nodes, they are ordered according
to the spatial index.
*/

//=============================================================================
//=============================================================================
/*!
\page amr_concurrency Concurrency

<!-------------------------------------------------------------------------->
\section amr_concurrency_distributed Distributed Memory 

The strategy for achieving distributed-memory concurrency is to partition the
nodes and distribute them across processors.  Each processor has an orthtree
that holds a subset of the nodes.  The domain of each orthtree is the same.
The nodes stored in a particular orthtree cover a portion of its domain.

There are many ways to partition the nodes across processors.  The simplest
method is to use the z-ordering and give each processor a contiguous block
of nodes.  The number of nodes in each processor is approximately equal.

<!-------------------------------------------------------------------------->
\section amr_concurrency_shared Shared Memory 

The strategy for achieving shared-memory concurrency is to partition a 
processor's nodes across threads.  Let there be \e N threads available in 
concurrent sections.  When applying an operation to all nodes, the set of 
nodes is partitioned into \e N subsets using the z-ordering.  Then each thread
acts on its set of nodes.  Here we assume that an operation applied to a node
may alter only that node.  However, it may access the state of other nodes.
We also assume that the operations may be applied in any order.  The threads 
may concurrenty alter the state of their nodes without conflict.

Refinement is more complicated than applying an operation to the nodes
because it involves inserting and erasing nodes.  One could make
insertions and deletions thread-safe by performing these operations in
OpenMP critical sections.  (Recall that only one thread may execute a
given critical section at a time.)  This appoach is easy to implement,
but may not be efficient because of its fine-grained nature.  It is
usually preferable for threads to cache their modifications and
perform them in a single critical section.  Note that determining
which nodes should be split may be done concurrently, each thread
examines its subset of the nodes.  After each thread collects the
nodes that need refinement, it applies splitting operations to those
nodes in a critical section.

The concurrent algorithm for coarsening is similar to that for refinement.
However, the partitioning must be done differently.  Coarsening is done
by repeating sweeps until no further coarsening is required.
In a coarsening sweep, a group of nodes may be merged a single time.
Before each sweep, the nodes are partitioned such that each potential
merging operation involves nodes in only one of the subsets.  The threads
concurrently collect the groups of nodes that need coarsening.  Then
each thread performs the merges in a critical section.

The figure below shows two partitions of a 2-D orthtree.  The first partition
evenly divides the nodes and is suitable for applying an operation to the 
node data or for refinement.  It is not suitable for coarsening because 
a mergeable group of nodes is divided between the two sets.  The second
partition can be used for a concurrent coarsening sweep.  Each set contains one
group of nodes that could be merged.

\image html partition.jpg "Two partitions of an orthtree."
*/



//=============================================================================
//=============================================================================
/*!
\page amr_performance Performance

- \ref amr_performance_AssociativeContainers
- \ref amr_performance_SpatialIndexMorton
- \ref amr_performance_serial

*/

//=============================================================================
//=============================================================================
/*!
\page amr_performance_AssociativeContainers

The container holds <code>size</code> elements which are pairs of 
<code>std::size_t</code>. The key and the value are the same and are the 
integers in the range <code>[0..size)</code>. For searching operations,
the keys are shuffled. The times in the tables below are given in nanoseconds.

Performance for <code>std::map</code>.
\htmlinclude AssociativeContainersMap.txt

Performance for <code>__gnu_cxx::hash_map</code>.
\htmlinclude AssociativeContainersHashMap.txt

Performance for <code>std::tr1::unordered_map</code>.
\htmlinclude AssociativeContainersUnorderedMap.txt
*/

//=============================================================================
//=============================================================================
/*!
\page amr_performance_SpatialIndexMorton Spatial Index Performance

\htmlinclude SpatialIndexMorton3.txt
*/

//=============================================================================
//=============================================================================
/*!
\page amr_performance_serial Serial Performance

<!-------------------------------------------------------------------------->
\section amr_performance_serial_access Node Access

<!-------------------------------------------------------------------------->
\subsection amr_performance_serial_access_ordered Ordered Access

Ordered access of the nodes.  This test is performed on an orthtree with
\f$2^{12} = 4096\f$ nodes.  The time per access is given in nanoseconds.

\htmlinclude OrthtreeAccessOrdered.txt

Ordered access of the nodes.  Next we test access times for a 3-D orthtree
of varying sizes.

\htmlinclude OrthtreeAccessOrdered3.txt

<!-------------------------------------------------------------------------->
\subsection amr_performance_serial_access_random Random Access

Random access of the nodes.  This test is performed on an orthtree with
\f$2^{12} = 4096\f$ nodes.  The time per access is given in nanoseconds.

\htmlinclude OrthtreeAccessRandom.txt

Random access of the nodes.  Next we test access times for a 3-D orthtree
of varying sizes.

\htmlinclude OrthtreeAccessRandom3.txt

<!-------------------------------------------------------------------------->
\section amr_performance_serial_refinement Refinement

Refinement operations.  This test is performed on an orthtree with
up to \f$2^{12} = 4096\f$ nodes.  The time per refinement is given in
nanoseconds.

\htmlinclude OrthtreeRefine.txt

Next we test refinement times for a 3-D orthtree of varying sizes.

\htmlinclude OrthtreeRefine3.txt


<!-------------------------------------------------------------------------->
\section amr_performance_serial_coarsening Coarsening

Coarsening operations.  This test is performed on an orthtree with
up to \f$2^{12} = 4096\f$ nodes.  The time per coarsening operation is given in
nanoseconds.

\htmlinclude OrthtreeCoarsen.txt

Next we test coarsening times for a 3-D orthtree of varying sizes.

\htmlinclude OrthtreeCoarsen3.txt


<!-------------------------------------------------------------------------->
\section amr_performance_serial_neighbors Finding Neighbors

Finding neighbors.  This test is performed on an orthtree with
\f$2^{12} = 4096\f$ nodes.  The time per operation is given in
nanoseconds.

\htmlinclude OrthtreeNeighbors.txt

Next we time the operations for a 3-D orthtree of varying sizes.

\htmlinclude OrthtreeNeighbors3.txt
*/



//=============================================================================
//=============================================================================
/*!
\page amr_compiling Compiling

To use this package, you will need a compiler that supports the C++ TR1 
library. GCC 4.2 and beyond should work. 
*/



//=============================================================================
//=============================================================================
/*!
\page amr_links Links

- http://en.wikipedia.org/wiki/Quadtree
- Visit http://donar.umiacs.umd.edu/quadtree/index.html for spatial index 
demos using Java.

- Chombo http://seesar.lbl.gov/anag/index.html
- Overture https://computation.llnl.gov/casc/Overture/ 
- HyperCLaw http://seesar.lbl.gov/ccse/Software/index.html 
*/

END_NAMESPACE_AMR

#endif
