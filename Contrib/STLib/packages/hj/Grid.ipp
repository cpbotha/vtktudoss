// -*- C++ -*-

#if !defined(__hj_Grid_ipp__)
#error This file is an implementation detail of the class Grid.
#endif

BEGIN_NAMESPACE_HJ

//
// Constructors
//


template <int N, typename T, class DifferenceScheme>
template <bool A>
inline
Grid<N,T,DifferenceScheme>::
Grid(ads::Array<N,Number,A>& solution, const Number dx, 
     const bool is_concurrent) :
  _solution(solution),
  // For the concurrent algorithm the solution grid has ghost points 
  // made by enlarging the grid by the radius of the finite difference 
  // scheme.  
  _ranges(_solution.lbounds() + (is_concurrent ? 
				   DifferenceScheme::radius() : 0),
	   _solution.ubounds() - (is_concurrent ? 
				   DifferenceScheme::radius() : 0)),
  _dx(dx),
  _scheme(_ranges, _solution, _dx, is_concurrent)
{}

  
template <int N, typename T, class DifferenceScheme>
inline
void
Grid<N,T,DifferenceScheme>::
initialize() {
  // Set the solution to infinity.
  _solution = std::numeric_limits<Number>::max();

  // Initialize the difference scheme.
  _scheme.initialize();
}


template <int N, typename T, class DifferenceScheme>
inline
int
Grid<N,T,DifferenceScheme>::
set_unsigned_initial_condition() {
  // Initialize the difference scheme.
  _scheme.initialize();

  // Loop over the solution array.
  Index i;
  int count = 0;
  const int size = _solution.size();
  for (int n = 0; n != size; ++n) {
    // If the solution is known.
    if (_solution[n] != std::numeric_limits<Number>::max()) {
      _solution.index_to_indices(n, i);
      // If the grid point has an unkown neighbor.
      if (_scheme.has_unknown_neighbor(i)) {
	// Set the status of the solution as known in the initial condition.
	_scheme.set_initial(i);
	++count;
      }
      else {
	// Set the status to KNOWN.
	_scheme.set_known(i);
      }
    }
  }
  return count;
}


template <int N, typename T, class DifferenceScheme>
inline
int
Grid<N,T,DifferenceScheme>::
set_negative_initial_condition() {
  // Initialize the difference scheme.
  _scheme.initialize();

  // Loop over the solution array.
  Index i;
  int initial = 0;
  const int size = _solution.size();
  for (int n = 0; n != size; ++n) {
    // If the solution is finite.
    if (_solution[n] != std::numeric_limits<Number>::max()) {
      _solution.index_to_indices(n, i);
      // If the solution is non-positive.
      if (_solution[n] <= 0 && _scheme.has_unknown_neighbor(i)) {
	// Set the status as known in the initial condition.
	_scheme.set_initial(i);
	++initial;
      }
      else {
	// Set the status to KNOWN.
	_scheme.set_known(i);
      }
      // Reverse the sign of the solution.
      _solution[n] = - _solution[n];
    }
  }
  return initial;
}


template <int N, typename T, class DifferenceScheme>
inline
int
Grid<N,T,DifferenceScheme>::
set_positive_initial_condition() {
  // Loop over the solution array.
  Index i;
  int initial = 0;
  const int size = _solution.size();
  for (int n = 0; n != size; ++n) {
    _solution.index_to_indices(n, i);
    // If the solution is known.
    if (_scheme.is_known(i)) {
#ifdef DEBUG_Grid
      assert(_solution[n] != std::numeric_limits<Number>::max());
#endif
      _solution[n] = - _solution[n];
      if (_solution[n] >= 0 && _scheme.has_unknown_neighbor(i)) {
	_scheme.set_initial(i);
	++initial;
      }
    }
  }
  return initial;
}


template <int N, typename T, class DifferenceScheme>
inline
void
Grid<N,T,DifferenceScheme>::
set_initial(const Index& i, const Number value) {
  _solution(i) = value;
  _scheme.set_initial(i);
}


template <int N, typename T, class DifferenceScheme>
inline
void
Grid<N,T,DifferenceScheme>::
add_source(const ads::FixedArray<2,Number>& x) {
  LOKI_STATIC_CHECK(N == 2, The_dimension_must_be_2);

  assert(_scheme.equation().radius() == 1 || 
	  _scheme.equation().radius() == 2);

  assert(_solution.lbound(0) <= x[0] && x[0] <= _solution.ubound(0) - 1 &&
	  _solution.lbound(1) <= x[1] && x[1] <= _solution.ubound(1) - 1);


  const int radius = _scheme.equation().radius() - 1;
  const Number eps = std::sqrt(std::numeric_limits<Number>::epsilon());
  // Floor.
  const int i_start = int(x[0] - eps) - radius;
  const int j_start = int(x[1] - eps) - radius;
  // Ceiling.
  const int i_stop = int(x[0] + eps) + 1 + radius;
  const int j_stop = int(x[1] + eps) + 1 + radius;
  Index i;
  for (i[1] = j_start; i[1] <= j_stop; ++i[1]) {
    for (i[0] = i_start; i[0] <= i_stop; ++i[0]) {
      if (indices_in_grid(i)) {
	set_initial(i, std::min(_solution(i), index_distance(i, x)));
      }
    }
  }
}


template <int N, typename T, class DifferenceScheme>
inline
void
Grid<N,T,DifferenceScheme>::
add_source(const ads::FixedArray<3,Number>& x) {
  LOKI_STATIC_CHECK(N == 3, The_dimension_must_be_3);

  assert(_scheme.equation().radius() == 1 || 
	  _scheme.equation().radius() == 2);

  assert(_solution.lbound(0) <= x[0] && x[0] <= _solution.ubound(0) - 1 &&
	  _solution.lbound(1) <= x[1] && x[1] <= _solution.ubound(1) - 1 &&
	  _solution.lbound(2) <= x[2] && x[2] <= _solution.ubound(2) - 1);


  const int radius = _scheme.equation().radius() - 1;
  const Number eps = std::sqrt(std::numeric_limits<Number>::epsilon());
  // Floor.
  const int i_start = int(x[0] - eps) - radius;
  const int j_start = int(x[1] - eps) - radius;
  const int k_start = int(x[2] - eps) - radius;
  // Ceiling.
  const int i_stop = int(x[0] + eps) + 1 + radius;
  const int j_stop = int(x[1] + eps) + 1 + radius;
  const int k_stop = int(x[2] + eps) + 1 + radius;
  Index i;
  for (i[2] = k_start; i[2] <= k_stop; ++i[2]) {
    for (i[1] = j_start; i[1] <= j_stop; ++i[1]) {
      for (i[0] = i_start; i[0] <= i_stop; ++i[0]) {
	if (indices_in_grid(i)) {
	  set_initial(i, std::min(_solution(i), index_distance(i, x)));
	}
      }
    }
  }
}


//
// I/O
//

  
template <int N, typename T, class DifferenceScheme>
inline
void
Grid<N,T,DifferenceScheme>::
print_statistics(std::ostream& out) const {
  // Print statistics for the status array.
  _scheme.print_statistics(out);

  T min_known = std::numeric_limits<T>::max();
  T max_known = -std::numeric_limits<T>::max();
  T s;
  int num_known = 0;
  int num_positive = 0;
  int num_nonpositive = 0;
  const int size = _solution.size();
  for (int n = 0; n != size; ++n) {
    s = _solution[n];
    if (s != std::numeric_limits<T>::max()) {
      ++num_known;

      if (s > 0) {
	++num_positive;
      }
      else {
	++num_nonpositive;
      }

      if (s < min_known) {
	min_known = s;
      }
      if (s > max_known) {
	max_known = s;
      }
    }
  }

  out << "Solution array size = " << size << "\n"
      << "  Number known =        " << num_known << "\n"
      << "  Number positive =     " << num_positive << "\n"
      << "  Number non-positive = " << num_nonpositive << "\n"
      << "  Minimum known =       " << min_known << "\n"
      << "  Maximum known =       " << max_known << "\n";
}


template <typename T, bool A>
inline
void
print_solution_array(std::ostream& out, const ads::Array<2,T,A>& solution) {
  // Write the solution.
  out << "Solution:" << '\n';
  for (int j = solution.ubound(1) - 1; j >= solution.lbound(1); --j){
    for (int i = solution.lbound(0); i < solution.ubound(0); ++i){
      out << solution(i, j) << " ";
    }
    out << '\n';
  }
}


template <typename T, bool A>
inline
void
print_solution_array(std::ostream& out, const ads::Array<3,T,A>& solution) {
  // Write the solution.
  out << "Solution:" << '\n';
  for (int k = solution.ubound(2) - 1; k >= solution.lbound(2); --k){
    for (int j = solution.ubound(1) - 1; j >= solution.lbound(1); --j){
      for (int i = solution.lbound(0); i < solution.ubound(0); ++i){
	out << solution(i, j, k) << " ";
      }
      out << '\n';
    }
    out << '\n';
  }
}

template <int N, typename T, class DifferenceScheme>
inline
void
Grid<N,T,DifferenceScheme>::
put(std::ostream& out) const {
  // Write the solution array.
  print_solution_array(out, _solution);
  // Write information from the difference scheme.
  out << _scheme;
}

//
// File I/O
//

template <int N, typename T, class DifferenceScheme>
inline
std::ostream& 
operator<<(std::ostream& out, const Grid<N,T,DifferenceScheme>& x) {
  x.put(out);
  return out;
}

END_NAMESPACE_HJ

// End of file.
