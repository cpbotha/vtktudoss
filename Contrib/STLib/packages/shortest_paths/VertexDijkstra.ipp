// -*- C++ -*-

#if !defined(__VertexDijkstra_ipp__)
#error This file is an implementation detail of the class VertexDijkstra.
#endif

BEGIN_NAMESPACE_SHORTEST_PATHS

//
// Constructors, Destructor.
//

template <typename WeightType>
inline
VertexDijkstra<WeightType>& 
VertexDijkstra<WeightType>::
operator=( const VertexDijkstra& rhs )
{
  if ( &rhs != this ) {
    base_type::operator=( rhs );
    _adjacent_edges = rhs._adjacent_edges;
    _heap_ptr = rhs._heap_ptr;
  }
  return *this;
}

END_NAMESPACE_SHORTEST_PATHS

// End of file.
