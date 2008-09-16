// -*- C++ -*-

#if !defined(__hj_DiffScheme_ipp__)
#error This file is an implementation detail of the class DiffScheme.
#endif

BEGIN_NAMESPACE_HJ


template <int N, typename T, class Equation>
template <bool A>
inline
DiffScheme<N,T,Equation>::
DiffScheme(const Range& index_ranges, ads::Array<N,Number,A>& solution, 
	    const Number dx, const bool is_concurrent) :
  _index_ranges(index_ranges),
  _solution(solution),
  // For the sequential algorithm, the status array is larger than the 
  // solution array.  It has an extra layer that is one grid point thick.
  // This is because grid points on the edge of the solution array will 
  // call label_neighbors.  For the concurrent algorithm, the status array
  // needs to have two extra grid points around the original solution array.
  // (By "original" I mean the solution array without the ghost points.  
  // This is because the GHOST_ADJACENT points will call label_neighbors.
  _status(Range(_index_ranges.lbounds() - (is_concurrent ? 2 : 1),
		       _index_ranges.ubounds() + (is_concurrent ? 2 : 1))),
  _equation(solution, _status, dx)
{}


template <int N, typename T, class Equation>
inline
void
DiffScheme<N,T,Equation>::
initialize() {
  //
  // Set the status of the grid points.
  //

  // First make everything a void point.
  _status = VOID;

  // Then mark the interior points an unlabeled.

  // Make an array without valid data.
  int* data = 0;
  ads::Array<N,int,false> a(_index_ranges, data);
  // Loop over the indices.
  ads::ArrayIndexIterator<N> begin = a.indices_begin();
  const ads::ArrayIndexIterator<N> end = a.indices_end();
  for (; begin != end; ++begin) {
    _status(*begin) = UNLABELED;
  }
}


template <int N, typename T, class Equation>
template <class Container>
inline
void
DiffScheme<N,T,Equation>::
label(Container& unlabeled_neighbors, const Index& i, const Number value) {
  Status& stat = _status(i);
  Number& soln = _solution(i);

#ifdef DEBUG_DiffScheme
  assert(stat == UNLABELED || stat == LABELED);
#endif

  if (stat == UNLABELED) {
    stat = LABELED;
    soln = value;
    unlabeled_neighbors.push_back(&soln);
  }
  else if (value < soln) {
    soln = value;
  }
}


template <int N, typename T, class Equation>
template <typename HandleType>
inline
void
DiffScheme<N,T,Equation>::
label(ads::PriorityQueueBinaryHeapArray<HandleType>& labeled, 
       const Index& i, const Number value) {
  Status& stat = _status(i);
  Number& soln = _solution(i);

#ifdef DEBUG_DiffScheme
  assert(stat == UNLABELED || stat == LABELED);
#endif

  if (stat == UNLABELED) {
    stat = LABELED;
    soln = value;
    labeled.push(&soln);
  }
  else if (value < soln) {
    soln = value;
    labeled.decrease(&soln);
  }
}


template <int N, typename T, class Equation>
template <typename HandleType>
inline
void
DiffScheme<N,T,Equation>::
label(ads::PriorityQueueBinaryHeapStoreKeys<HandleType>& labeled, 
       const Index& i, const Number value) {
  // CONTINUE: REMOVE
  //std::cerr << "label(labeled, " << i << ", " << value << ")\n";

  Status& stat = _status(i);
  Number& soln = _solution(i);

#ifdef DEBUG_DiffScheme
  assert(stat == UNLABELED || stat == LABELED);
#endif

  if (stat == UNLABELED) {
    stat = LABELED;
    soln = value;
    labeled.push(&soln);
  }
  else if (value < soln) {
    soln = value;
    labeled.push(&soln);
  }
  //std::cerr << "Done label().\n";
}


template <int N, typename T, class Equation>
template <typename HandleType>
inline
void
DiffScheme<N,T,Equation>::
label(ads::PriorityQueueCellArray<HandleType>& labeled, 
       const Index& i, const Number value) {
  Status& stat = _status(i);
  Number& soln = _solution(i);

#ifdef DEBUG_DiffScheme
  assert(stat == UNLABELED || stat == LABELED);
#endif

  if (stat == UNLABELED) {
    stat = LABELED;
    soln = value;
    labeled.push(&soln);
  }
  else if (value < soln) {
    soln = value;
    labeled.push(&soln);
  }
}


template <int N, typename T, class Equation>
inline
void
DiffScheme<N,T,Equation>::
label(int, const Index& i, const Number value) {
  Status& stat = _status(i);
  Number& soln = _solution(i);

#ifdef DEBUG_DiffScheme
  assert(stat == UNLABELED || stat == LABELED);
#endif

  if (stat == UNLABELED) {
    stat = LABELED;
    soln = value;
  }
  else if (value < soln) {
    soln = value;
  }
}


//
// File I/O
//


template <int N, typename T, class Equation>
inline
void
DiffScheme<N,T,Equation>::
print_statistics(std::ostream& out) const {
  int num_known = 0;
  int num_labeled = 0;
  int num_unlabeled = 0;
  int num_void = 0;
  int num_initial = 0;
  Status s;
  // Loop over the status array.
  const int size = _status.size();
  for (int n = 0; n != size; ++n) {
    s = _status[n];
    if (s == KNOWN) {
      ++num_known;
    }
    else if (s == LABELED) {
      ++num_labeled;
    }
    else if (s == UNLABELED) {
      ++num_unlabeled;
    }
    else if (s == VOID) {
      ++num_void;
    }
    else if (s == INITIAL) {
      ++num_initial;
    }
    else {
      assert(false);
    }
  }

  out << "Status array size = " << size << "\n"
      << "  Number KNOWN =     " << num_known << "\n"
      << "  Number LABELED =   " << num_labeled << "\n"
      << "  Number UNLABELED = " << num_unlabeled << "\n"
      << "  Number VOID =      " << num_void << "\n"
      << "  Number INITIAL =   " << num_initial << "\n";
}


inline 
char
status_character(const Status s) {
  if (s == KNOWN) {
    return 'K';
  }
  else if (s == LABELED) {
    return 'L';
  }
  else if (s == UNLABELED) {
    return 'U';
  }
  else if (s == VOID) {
    return 'V';
  }
  else if (s == INITIAL) {
    return 'I';
  }
  else {
    assert(false);
  }
  return ' ';
}


template <bool A>
inline
void
print_status_array(std::ostream& out, const ads::Array<2,Status,A>& status) {
  out << "Status:" << '\n';
  for (int j = status.ubound(1) - 1; j >= status.lbound(1); --j) {
    for (int i = status.lbound(0); i < status.ubound(0); ++i) {
      out << status_character(status(i, j)) << ' ';
    }
    out << '\n';
  }
}


template <bool A>
inline
void
print_status_array(std::ostream& out, const ads::Array<3,Status,A>& status) {
  Status s;
  out << "Status:" << '\n';
  for (int k = status.ubound(2) - 1; k >= status.lbound(2); --k) {
    for (int j = status.ubound(1) - 1; j >= status.lbound(1); --j) {
      for (int i = status.lbound(0); i < status.ubound(0); ++i) {
	out << status_character(status(i, j, k)) << ' ';
      }
      out << '\n';
    }
    out << '\n';
  }
}


template <int N, typename T, class Equation>
inline
void
DiffScheme<N,T,Equation>::
put(std::ostream& out) const {
  print_status_array(out, _status);
  out << _equation;
}


template <int N, typename T, class Equation>
std::ostream& 
operator<<(std::ostream& out, const DiffScheme<N,T,Equation>& x) {
  x.put(out);
  return out;
}


END_NAMESPACE_HJ

// End of file.
