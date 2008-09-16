// -*- C++ -*-

#if !defined(__hj_GridMCC_ipp__)
#error This file is an implementation detail of the class GridMCC.
#endif

BEGIN_NAMESPACE_HJ

template <int N, typename T, class DifferenceScheme>
inline
void
GridMCC<N,T,DifferenceScheme>::
solve( const Number max_solution )
{
  Index i;

  // The list of labeled unknown grid points to check during a step.
  std::vector<Number*> labeled;
  // The list of grid points which are added to labeled during a step.
  std::vector<Number*> new_labeled;

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
      
  typename std::vector<Number*>::iterator grid_ptr_iter, labeled_end;
  Number min_unknown = 0;
  ads::LessByHandle<Number*> grid_pt_compare;

  // If we are going to solve for all grid points.
  if ( max_solution == 0 ) {
    // All vertices are known when there are no labeled vertices left.
    // Loop while there are labeled vertices left.
    while ( ! labeled.empty() ) {
      // Find the minimum unknown grid point.
      min_unknown = **(std::min_element( labeled.begin(), labeled.end(), 
					 grid_pt_compare ));

      // Loop through the labeled grid points.
      labeled_end = labeled.end();
      for ( grid_ptr_iter = labeled.begin();
	    grid_ptr_iter != labeled_end;
	    ++grid_ptr_iter ) {
	// Get the indices of the grid point.
	_solution.iterator_to_indices( *grid_ptr_iter, i );
	if ( _solution( i ) <= 
	     _scheme.lower_bound( i, min_unknown ) ) {
	  // Mark the grid point as known and label the neighboring 
	  // grid points.
	  _scheme.label_neighbors( new_labeled, i );
	  // Flag the grid point for deletion from the labeled set.
	  *grid_ptr_iter = 0;
	}
      }
    
      // Remove the grid points that became known.
      labeled.erase( std::remove( labeled.begin(), labeled.end(),
				  static_cast<Number*>(0) ), 
		     labeled.end() );
      // Add the newly labeled vertices.
      labeled.insert( labeled.end(), new_labeled.begin(), 
		      new_labeled.end() );
      new_labeled.clear();
    }
  }
  // Else we solve for the grid points around the initial condition.
  else {
    // Loop while there are labeled vertices left and while the solution
    // is less than or equal to max_solution.
    while ( ! labeled.empty() && min_unknown <= max_solution ) {
      // Find the minimum unknown grid point.
      min_unknown = **(std::min_element( labeled.begin(), labeled.end(), 
					 grid_pt_compare ));

      // Loop through the labeled grid points.
      labeled_end = labeled.end();
      for ( grid_ptr_iter = labeled.begin();
	    grid_ptr_iter != labeled_end;
	    ++grid_ptr_iter ) {
	// Get the indices of the grid point.
	_solution.iterator_to_indices( *grid_ptr_iter, i );
	if ( _solution( i ) <= 
	     _scheme.lower_bound( i, min_unknown ) ) {
	  // Mark the grid point as known and label the neighboring 
	  // grid points.
	  _scheme.label_neighbors( new_labeled, i );
	  // Flag the grid point for deletion from the labeled set.
	  *grid_ptr_iter = 0;
	}
      }
    
      // Remove the grid points that became known.
      labeled.erase( std::remove( labeled.begin(), labeled.end(),
				  static_cast<Number*>(0) ), 
		     labeled.end() );
      // Add the newly labeled vertices.
      labeled.insert( labeled.end(), new_labeled.begin(), 
		      new_labeled.end() );
      new_labeled.clear();
    }
  }
}

END_NAMESPACE_HJ

// End of file.
