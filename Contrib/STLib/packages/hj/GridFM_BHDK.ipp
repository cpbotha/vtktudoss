// -*- C++ -*-

#if !defined(__hj_GridFM_BHDK_ipp__)
#error This file is an implementation detail of the class GridFM_BHDK.
#endif

BEGIN_NAMESPACE_HJ

template <int N, typename T, class DifferenceScheme>
inline
void
GridFM_BHDK<N,T,DifferenceScheme>::
solve( const Number max_solution )
{
  Index i;

  // The heap of labeled unknown grid points.
  ads::PriorityQueueBinaryHeapArray<const_handle> labeled( _solution );
    
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
      labeled.pop();
      // Label the neighbors.
      _scheme.label_neighbors( labeled, i );
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
      labeled.pop();
      // Label the neighbors.
      _scheme.label_neighbors( labeled, i );
    }
  }
}

END_NAMESPACE_HJ

// End of file.
