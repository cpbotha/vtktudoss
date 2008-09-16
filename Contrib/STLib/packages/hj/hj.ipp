// -*- C++ -*-

#if !defined(__hj_hj_ipp__)
#error This file is an implementation detail of hj.
#endif

BEGIN_NAMESPACE_HJ


template<int N, typename T, bool A>
inline
void
computeUnsignedDistance(ads::Array<N,T,A>& array, const T dx, 
			const T maximumDistance) {
#ifdef DEBUG_hj
  assert(ads::computeProduct(array.extents()) > 0);
  assert(maximumDistance >= 0);
#endif
  GridMCC< N, T, DiffSchemeAdjDiag< N, T, DistanceAdjDiag1st<N,T> > > 
    grid(array, dx);
  grid.set_unsigned_initial_condition();
  grid.solve(maximumDistance);
}


template<int N, typename T, bool A>
inline
void
computeSignedDistance(ads::Array<N,T,A>& array, const T dx, const T maximumDistance) {
#ifdef DEBUG_hj
  assert(ads::computeProduct(array.extents()) > 0);
  assert(maximumDistance >= 0);
#endif
  // CONTINUE: REMOVE
  GridMCC< N, T, DiffSchemeAdjDiag< N, T, DistanceAdjDiag1st<N,T> > > 
    grid(array, dx);

  std::cerr << "grid.set_negative_initial_condition()\n";
  std::cerr << grid.set_negative_initial_condition() << "\n";
  grid.print_statistics(std::cerr);

  std::cerr << "grid.solve(maximumDistance)\n";
  grid.solve(maximumDistance);
  grid.print_statistics(std::cerr);

  std::cerr << "grid.set_positive_initial_condition()\n";
  std::cerr << grid.set_positive_initial_condition() << "\n";
  grid.print_statistics(std::cerr);

  std::cerr << "grid.solve(maximumDistance)\n";
  grid.solve(maximumDistance);
  grid.print_statistics(std::cerr);
}


template<int N, typename T, bool A>
inline
void
floodFillUnsignedDistance(ads::Array<N,T,A>& array,
			  const T maximumDistance, const T fillValue) {
#ifdef DEBUG_hj
  assert(ads::computeProduct(array.extents()) > 0);
  assert(maximumDistance > 0);
#endif
  // The default value of fillValue.
  if (fillValue == 0) {
    fillValue = maximumDistance;
  }

  typename ads::Array<N,T>::iterator i = array.begin();
  const typename ads::Array<N,T>::iterator end = array.end();
  for (; i != end; ++i) {
    if (*i > maximumDistance) {
      *i = fillValue;
    }
  }
}


template<typename T, bool A>
inline
void
floodFillSignedDistance(ads::Array<2,T,A>& array,
			const T maximumDistance, const T fillValue) {
#ifdef DEBUG_hj
  assert(ads::computeProduct(array.extents()) > 0);
  assert(maximumDistance > 0);
#endif

  // The default value of fillValue.
  if (fillValue == 0) {
    fillValue = maximumDistance;
  }

  //
  // See if the distance is known for any grid points
  //
  bool result = false;
  int sign = 0;
  const int i_begin = array.lbound(0);
  const int i_end = array.ubound(0);
  const int j_begin = array.lbound(1);
  const int j_end = array.ubound(1);
  int i, j;
  for (j = j_begin; !result && j != j_end; ++j) {
    for (i = i_begin; !result && i != i_end; ++i) {
      if (array(i, j) != std::numeric_limits<T>::max()) {
	result = true;
	sign = (array(i, j) > 0) ? 1 : -1;
      }
    }
  }

  //
  // If there are any points in a known distance.
  //
  if (result) {
    int ysign = sign;

    //
    // Flood fill the distance with +- far_away.
    //
    for (j = j_begin; j != j_end; ++j) {
      if (array(0, j) != std::numeric_limits<T>::max()) {
	ysign = (array(0, j) > 0) ? 1 : -1;
      }
      sign = ysign;
      for (i = i_begin; i != i_end; ++i) {
	if (array(i, j) != std::numeric_limits<T>::max()) {
	  sign = (array(i, j) > 0) ? 1 : -1;
	}
	else {
	  // Set the distance to +- far_away.
	  array(i, j) = sign * fillValue;
	}
      }
    }
  } // end if (result)
  else {
    // Set the distance to +fillValue.
    array = fillValue;
  }
}


template<typename T, bool A>
inline
void
floodFillSignedDistance(ads::Array<3,T,A>& array,
			const T maximumDistance, const T fillValue) {
#ifdef DEBUG_hj
  assert(ads::computeProduct(array.extents()) > 0);
  assert(maximumDistance > 0);
#endif

  // The default value of fillValue.
  if (fillValue == 0) {
    fillValue = maximumDistance;
  }

  //
  // See if the distance is known for any grid points
  //
  bool result = false;
  int sign = 0;
  const int i_begin = array.lbound(0);
  const int i_end = array.ubound(0);
  const int j_begin = array.lbound(1);
  const int j_end = array.ubound(1);
  const int k_begin = array.lbound(2);
  const int k_end = array.ubound(2);
  int i, j, k;
  for (k = k_begin; !result && k != k_end; ++k) {
    for (j = j_begin; !result && j != j_end; ++j) {
      for (i = i_begin; !result && i != i_end; ++i) {
	if (array(i,j,k) != std::numeric_limits<T>::max()) {
	  result = true;
	  sign = (array(i,j,k) > 0) ? 1 : -1;
	}
      }
    }
  }

  //
  // If there are any points in a known distance.
  //
  if (result) {
    int ysign = sign, zsign = sign;

    //
    // Flood fill the distance with +- far_away.
    //
    for (k = k_begin; k != k_end; ++k) {
      if (array(0,0,k) != std::numeric_limits<T>::max()) {
	zsign = (array(0,0,k) > 0) ? 1 : -1;
      }
      ysign = zsign;
      for (j = j_begin; j != j_end; ++j) {
	if (array(0,j,k) != std::numeric_limits<T>::max()) {
	  ysign = (array(0,j,k) > 0) ? 1 : -1;
	}
	sign = ysign;
	for (i = i_begin; i != i_end; ++i) {
	  if (array(i,j,k) != std::numeric_limits<T>::max()) {
	    sign = (array(i,j,k) > 0) ? 1 : -1;
	  }
	  else {
	    // Set the distance to +-fillValue.
	    array(i,j,k) = sign * fillValue;
	  }
	}
      }
    }
  } // end if (result)
  else {
    // Set the distance to +far_away.
    array = fillValue;
  }
}


namespace neighbor {

  const int adj_x[4] = {1, 0, -1, 0};
  const int adj_y[4] = {0, 1, 0, -1};

  const int diag_x[4] = {1, -1, -1, 1};
  const int diag_y[4] = {1, 1, -1, -1};

  const int any_x[8] = {1, 1, 0, -1, -1, -1, 0, 1};
  const int any_y[8] = {0, 1, 1, 1, 0, -1, -1, -1};
}


template<typename T, bool A>
inline
bool
is_in_narrow_band(const ads::Array<2,T,A>& array, const int index) {
#ifdef DEBUG_hj
  assert(array[index] != std::numeric_limits<T>::max());
#endif

  int i, j;
  array.index_to_indices(index, i, j);
  int x, y;
  for (int n = 0; n != 8; ++n) {
    x = i + neighbor::any_x[n];
    y = j + neighbor::any_y[n];
    if (array.ranges().is_in(x, y) && 
	array(i, j) * array(x, y) <= 0) {
      return true;
    }
  }
  return false;
}


template<typename T, bool A>
inline
T
max_derivative(const ads::Array<2,T,A>& array, const T dx, 
	       const int index) {
  const T a = array[index];
#ifdef DEBUG_hj
  assert(a != std::numeric_limits<T>::max());
  assert(dx > 0);
#endif
  static const T sqrt_2_inverse = 1.0 / std::sqrt(T(2));
    
  int i, j;
  array.index_to_indices(index, i, j);
  T d, delta = 0;
  int x, y;
  // The four adjacent and four diagonal directions.
  for (int n = 0; n != 4; ++n) {
    x = i + neighbor::adj_x[n];
    y = j + neighbor::adj_y[n];
    if (array.ranges().is_in(x, y)) {
      d = std::abs(array(x, y) - a);
      if (d > delta) {
	delta = d;
      }
    }
    x = i + neighbor::diag_x[n];
    y = j + neighbor::diag_y[n];
    if (array.ranges().is_in(x, y)) {
      d = std::abs(array(x, y) - a) * sqrt_2_inverse;
      if (d > delta) {
	delta = d;
      }
    }
  }
  // Change the difference to a derivative.
  delta /= dx;
    
  if (delta == 0) {
    // This should occur only if this point and all its known neighbors are
    // zero.  In this case, set the derivative to unity to avoid division
    // by zero.
    delta = 1;
  }
  return delta;
}


namespace neighbor {

  const int num_adjacent = 6;
  const int adjacent_x[6] = { 0, 0,-1, 1, 0, 0};
  const int adjacent_y[6] = { 0,-1, 0, 0, 1, 0};
  const int adjacent_z[6] = {-1, 0, 0, 0, 0, 1};

  const int num_diagonal = 12;
  const int diagonal_x[12] = { 0,-1, 1, 0,-1, 1,-1, 1, 0,-1, 1, 0};
  const int diagonal_y[12] = {-1, 0, 0, 1,-1,-1, 1, 1,-1, 0, 0, 1};
  const int diagonal_z[12] = {-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1};

  const int num_corner = 8;
  const int corner_x[8] = {-1, 1,-1, 1,-1, 1,-1, 1};
  const int corner_y[8] = {-1,-1, 1, 1,-1,-1, 1, 1};
  const int corner_z[8] = {-1,-1,-1,-1, 1, 1, 1, 1};

  const int num_neighbor = 26;
  const int neighbor_x[26] = {-1, 0, 1,-1, 0, 1,-1, 0, 1,
			      -1, 0, 1,-1, 1,-1, 0, 1,
			      -1, 0, 1,-1, 0, 1,-1, 0, 1};
  const int neighbor_y[26] = {-1,-1,-1, 0, 0, 0, 1, 1, 1,
			      -1,-1,-1, 0, 0, 1, 1, 1,
			      -1,-1,-1, 0, 0, 0, 1, 1, 1};
  const int neighbor_z[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,
			      0, 0, 0, 0, 0, 0, 0, 0,
			      1, 1, 1, 1, 1, 1, 1, 1, 1};
}


template<typename T, bool A>
inline
bool
is_in_narrow_band(const ads::Array<3,T,A>& array, 
		  const int index) {
  const T a = array[index];
#ifdef DEBUG_hj
  assert(a != std::numeric_limits<T>::max());
#endif
  int i, j, k;
  array.index_to_indices(index, i, j, k);
  int x, y, z;
  for (int n = 0; n != neighbor::num_neighbor; ++n) {
    x = i + neighbor::neighbor_x[n];
    y = j + neighbor::neighbor_y[n];
    z = k + neighbor::neighbor_z[n];
    if (array.ranges().is_in(x, y, z) && a * array(x, y, z) <= 0) {
      return true;
    }
  }
  return false;
}


template<typename T, bool A>
inline
T
max_derivative(const ads::Array<3,T,A>& array, const T dx, 
	       const int index) {
  const T a = array[index];
#ifdef DEBUG_hj
  assert(a != std::numeric_limits<T>::max());
  assert(dx > 0);
#endif
  static const T sqrt_2_inverse = 1.0 / std::sqrt(T(2));
  static const T sqrt_3_inverse = 1.0 / std::sqrt(T(3));
    
  int i, j, k;
  array.index_to_indices(index, i, j, k);
  T d, delta = 0;
  int x, y, z;

  // The adjacent directions.
  for (int n = 0; n != neighbor::num_adjacent; ++n) {
    x = i + neighbor::adjacent_x[n];
    y = j + neighbor::adjacent_y[n];
    z = k + neighbor::adjacent_z[n];
    if (array.ranges().is_in(x, y, z)) {
      d = std::abs(array(x, y, z) - a);
      if (d > delta) {
	delta = d;
      }
    }
  }

  // The diagonal directions.
  for (int n = 0; n != neighbor::num_diagonal; ++n) {
    x = i + neighbor::diagonal_x[n];
    y = j + neighbor::diagonal_y[n];
    z = k + neighbor::diagonal_z[n];
    if (array.ranges().is_in(x, y, z)) {
      d = std::abs(array(x, y, z) - a) * sqrt_2_inverse;
      if (d > delta) {
	delta = d;
      }
    }
  }

  // The corner directions.
  for (int n = 0; n != neighbor::num_corner; ++n) {
    x = i + neighbor::corner_x[n];
    y = j + neighbor::corner_y[n];
    z = k + neighbor::corner_z[n];
    if (array.ranges().is_in(x, y, z)) {
      d = std::abs(array(x, y, z) - a) * sqrt_3_inverse;
      if (d > delta) {
	delta = d;
      }
    }
  }

  // Change the difference to a derivative.
  delta /= dx;
    
  if (delta == 0) {
    // This should occur only if this point an all its known neighbors are
    // zero.  In this case, set the derivative to unity to avoid division
    // by zero.
    delta = 1;
  }
  return delta;
}


template<int N, typename T, bool A>
inline
void
convertLevelSetToSignedDistance(ads::Array<N,T,A>& array, const T dx,
				const T isoValue, const T maximumDistance, 
				const T fillValue) {
#ifdef DEBUG_hj
  assert(ads::computeProduct(array.extents()) > 0);
  assert(maximumDistance >= 0);
#endif

  // The default value of fillValue.
  if (fillValue == 0) {
    fillValue = maximumDistance;
  }

  // Adjust the field of values so the iso-curve is zero.
  if (isoValue != 0) {
    array -= isoValue;
  }
    
  {
    // The grid points in the narrow band around the zero iso-curve.
    std::vector<int> narrow_band;
    std::vector<T> narrow_band_value;

    // Loop over the array, looking for grid points that neighbor the zero
    // iso-curve.
    const int size = array.size();
    for (int n = 0; n != size; ++n) {
      if (is_in_narrow_band(array, n)) {
	narrow_band.push_back(n);
	narrow_band_value.push_back(array[n] / 
				    max_derivative(array, dx, n));
      }
    }

    // Set all the grid points to infinity.
    array = std::numeric_limits<T>::max();

    // Make the points in the narrow band the initial condition.
    const int nbv_size = narrow_band_value.size();
    for (int n = 0; n != nbv_size; ++n) {
      array[ narrow_band[n] ] = narrow_band_value[n];
    }
  }

  // Compute the signed distance from the initial condition.
  computeSignedDistance(array, dx, maximumDistance);

  // Flood fill the signed distance if necessary.
  if (maximumDistance != 0) {
    floodFillSignedDistance(array, maximumDistance, fillValue);
  }
}










template<typename T, typename F, bool A1, bool A2>
inline
void
advectConstantIntoNegativeDistance(ads::Array<3,F,A1>& field,
				   const geom::RegularGrid<3,T>& grid,
				   const ads::Array<3,T,A2>& distance, 
				   const T maximumDistance,
				   const F defaultValue) {
  // The point type.
  typedef ads::FixedArray<3,T> Point;
  // The field array type.
  typedef ads::Array<3,F,A1> FieldArray;
  // An array of numbers.
  typedef ads::Array<3,T,A2> NumberArray;
  // A multi-index into a 3-D array.
  typedef typename NumberArray::index_type Index;

  assert(field.extents() == grid.getExtents());
  assert(distance.extents() == grid.getExtents());

  // The grid spacing.
  const Point delta = grid.getDelta();
  const Point deltaInverse(1.0 / delta[0], 1.0 / delta[1], 1.0 / delta[2]);

  // If the distance is known.
  ads::Array<3,bool> 
    is_known(ads::Array<3,bool>::
	     range_type(distance.lbounds() - 1, distance.ubounds() + 1)); 
  is_known = false;
  
  //
  // Find the grid points in the ghost fluid region.
  // Make a vector of iterators on these points.
  // Find the grid points with known distances.
  // Set the field values far inside the solid to the default value.
  //
  typedef typename NumberArray::const_iterator grid_point_const_iterator;
  typedef std::vector< grid_point_const_iterator > grid_point_container;
  grid_point_container ghost_points;
  {
    int i, j, k;
    const typename NumberArray::const_iterator iter_end = distance.end();
    for (typename NumberArray::const_iterator iter = distance.begin();
	 iter != iter_end; ++iter) {
      // If the distance is known.
      if (-maximumDistance < *iter && *iter < maximumDistance) {
	distance.iterator_to_indices(iter, i, j, k);
	is_known(i, j, k) = true;
	// If this point is in the ghost fluid region.
	if (*iter < 0) {
	  ghost_points.push_back(iter);
	}
      }
    }
  }

  //
  // Set the field values for points far inside the ghost fluid region.
  //
  {
    // Loop over the distance array and the field array.
    const typename NumberArray::const_iterator dist_iter_end = distance.end();
    typename NumberArray::const_iterator dist_iter = distance.begin();
    typename FieldArray::iterator field_iter = field.begin();
    for (; dist_iter != dist_iter_end; ++dist_iter, ++field_iter) {
      // Negative distance, not close to the boundary.
      if (*dist_iter <= - maximumDistance) {
	*field_iter = defaultValue;
      }
    }
  }


  //
  // Sort the grid points in the ghost fluid region by distance.
  //
  std::sort(ghost_points.begin(), ghost_points.end(),
	    ads::greater_handle<typename NumberArray::const_iterator>());

  //
  // Compute the upwind directions.
  //
  std::vector< Index > upwind;
  upwind.reserve(ghost_points.size());
  std::vector< Point > gradient;
  gradient.reserve(ghost_points.size());
  {
    int i, j, k, in, jn, kn, ip, jp, kp;
    int f[3]; // offset.
    bool is_known_n, is_known_p;
    T d, dn = 0, dp = 0;
    Index direction;
    Point grad;
    T mag;
    typename grid_point_container::const_iterator 
      iter_end = ghost_points.end();
    for (typename grid_point_container::const_iterator 
	   iter = ghost_points.begin(); 
	 iter != iter_end; ++iter) {
      distance.iterator_to_indices(*iter, i, j, k);

      for (int n = 0; n != 3; ++n) {
	f[0] = f[1] = f[2] = 0;
	f[n] = 1;
	in = i-f[0]; jn = j-f[1]; kn = k-f[2];
	ip = i+f[0]; jp = j+f[1]; kp = k+f[2];
	d = distance(i, j, k);
	if (is_known_n = is_known(in, jn, kn)) {
	  dn = distance(in, jn, kn);
	}
	if (is_known_p = is_known(ip, jp, kp)) {
	  dp = distance(ip, jp, kp);
	}

	direction[n] = 0;
	grad[n] = 0;
	if (is_known_n && is_known_p) {
	  if (dn > dp) {
	    if (dn > d) {
	      direction[n] = -1;
	      grad[n] = (dn - d) * deltaInverse[n];
	    }
	  }
	  else {
	    if (dp > d) {
	      direction[n] = 1;
	      grad[n] = (dp - d) * deltaInverse[n];
	    }
	  }
	}
	else if (is_known_n && dn > d) {
	  direction[n] = -1;
	  grad[n] = (dn - d) * deltaInverse[n];
	}
	else if (is_known_p && dp > d) {
	  direction[n] = 1;
	  grad[n] = (dp - d) * deltaInverse[n];
	}
      }

      // CONTINUE REMOVE
      /*
	if (direction == Index(0, 0, 0)) {
	std::cout << "direction = " << direction << '\n'
	<< "grad = " << grad << '\n'
	<< "distance(" << i << ", " << j << ", " << k << ") = "
	<< distance(i, j, k) << '\n';
	}
      */

      upwind.push_back(direction);
      // Normalize the gradient.
      mag = geom::computeMagnitude(grad);
      if (mag != 0) {
	grad /= mag;
      }
      gradient.push_back(grad);
    }
  }

  //
  // Apply the difference scheme.
  //
  T f[3], g[3], d[3];
  int i[3], m;
  Index direction;
  const Index zero_direction(0, 0, 0);
  Point grad;

  //
  // Loop over the grid points in the ghost fluid region.
  //
  const int sz = ghost_points.size();
  for(int n = 0; n != sz; ++n) {
    // The grid indices.
    distance.iterator_to_indices(ghost_points[n], i[0], i[1], i[2]);
    // The upwind direction.
    direction = upwind[n];

    // If there are any upwind directions.
    if (direction != zero_direction) {
      // For each direction.
      for (m = 0; m != 3; ++m) {
	if (direction[m]) {
	  i[m] += direction[m];
	  f[m] = field(i[0], i[1], i[2]);
	  i[m] -= direction[m];
	  d[m] = delta[m];
	  g[m] = gradient[n][m];
	}
	else {
	  f[m] = 0;
	  d[m] = 1;
	  g[m] = 0;
	}
      }

      assert(! (g[0] == 0 && g[1] == 0 && g[2] == 0));

      // A first order, upwind, finite difference scheme to solve:
      // (grad field) . (grad distance) = 0
      field(i[0], i[1], i[2]) = (f[0] * g[0] * d[1] * d[2] + 
				 f[1] * g[1] * d[0] * d[2] + 
				 f[2] * g[2] * d[0] * d[1]) / 
	(g[0] * d[1] * d[2] + d[0] * g[1] * d[2] + d[0] * d[1] * g[2]);
    }
    // Else there are no upwind coordinate directions.
    else {
      field(i[0], i[1], i[2]) = defaultValue;
    }
  }
}

END_NAMESPACE_HJ

// End of file.
