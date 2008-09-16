// -*- C++ -*-

#if !defined(__hj_GridFM_BH_ipp__)
#error This file is an implementation detail of the class GridFM_BH.
#endif

BEGIN_NAMESPACE_HJ

template <int N, typename T, class DifferenceScheme>
inline
void
GridFM_BH<N,T,DifferenceScheme>::
solve( const Number max_solution )
{
  Index i;

  // The heap of labeled unknown grid points.
  ads::PriorityQueueBinaryHeapStoreKeys<const_handle> labeled;

  // Label the neighbors of known grid points.
  {
    const int size = _solution.size();
    for ( int n = 0; n != size; ++n ) {
      _solution.index_to_indices( n, i );
      if ( _scheme.is_initial( i ) ) {
	_scheme.label_neighbors( labeled, i );
      }
    }
  }
      
  // If we are going to solve for all grid points.
  if ( max_solution == 0 ) {
    // All vertices are known when there are no labeled vertices left.
    // Loop while there are labeled vertices left.
    while ( ! labeled.empty() ) {
      // Get a handle to the grid point with minimum solution.
      // Convert the handle to indices.
      _solution.iterator_to_indices( labeled.top(), i );
      // Remove the element at the top of the heap.
      labeled.pop();
      // If the grid point has not already been set.
      if ( _scheme.status()( i ) == LABELED ) {
	// Label the neighbors.
	_scheme.label_neighbors( labeled, i );
      }
    }
  }
  // Else we solve for the grid points around the initial condition.
  else {
    // Loop while there are labeled vertices left and while the solution
    // is less than or equal to max_solution.
    Number soln = max_solution;
    while ( ! labeled.empty() && soln <= max_solution ) {
      // Get a handle to the grid point with minimum solution.
      // Convert the handle to indices.
      _solution.iterator_to_indices( labeled.top(), i );
      soln = _solution( i );
      // Remove the element at the top of the heap.
      labeled.pop();
      // If the grid point has not already been set.
      if ( _scheme.status()( i ) == LABELED ) {
	// Label the neighbors.
	_scheme.label_neighbors( labeled, i );
      }
    }
  }
}

END_NAMESPACE_HJ

// End of file.
