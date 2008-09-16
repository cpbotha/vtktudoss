// -*- C++ -*-

#if !defined(__partition_BspTree_ipp__)
#error This file is an implementation detail of BspTree.
#endif

BEGIN_NAMESPACE_CONCURRENT

//---------------------------------------------------------------------------
// 1-D
//---------------------------------------------------------------------------

template<typename _T>
inline
void
partitionRegularGridWithBspTree
(const ads::Array<1, _T>& costs, ads::Array<1, int>* identifiers,
 const _T totalCost, geom::SemiOpenInterval<1, int> indexRange,
 int identifiersBegin,  int identifiersEnd) {
  // If we are down to a single piece or there are no elements.
  if (identifiersEnd - identifiersBegin <= 1 || indexRange.isEmpty()) {
    // Do nothing.
    return;
  }

  const int identifiersMiddle = (identifiersBegin + identifiersEnd) / 2;
  const _T splittingCost = totalCost * (identifiersMiddle - identifiersBegin) /
    _T(identifiersEnd - identifiersBegin);
  
  // Find the splitting index.
  const int begin = indexRange.getLowerCorner()[0];
  const int end = indexRange.getUpperCorner()[0];
  _T cost = 0;
  int splittingIndex = begin;
  for (; splittingIndex != end; ++splittingIndex) {
    // If we should put this element in the lower partition.
    if (cost + costs(splittingIndex) / 2 < splittingCost) {
      // Place the element in the lower partition.
      (*identifiers)(splittingIndex) = identifiersMiddle - 1;
      cost += costs(splittingIndex);
    }
    else {
      // Stop partitioning.
      break;
    }
  }

  // Recurse.
  const int upper = indexRange.getUpperCorner()[0];
  indexRange.setUpperCoordinate(0, splittingIndex);
  partitionRegularGridWithBspTree(costs, identifiers, cost, indexRange,
				  identifiersBegin, identifiersMiddle);
  indexRange.setLowerCoordinate(0, splittingIndex);
  indexRange.setUpperCoordinate(0, upper);
  partitionRegularGridWithBspTree(costs, identifiers, totalCost - cost, 
				  indexRange, identifiersMiddle,
				  identifiersEnd);
}



template<typename _T>
inline
void
partitionRegularGridWithBspTree(const ads::Array<1, _T>& costs,
				ads::Array<1, int>* identifiers,
				const int numberOfPartitions) {
  assert(costs.ranges() == identifiers->ranges());
  // Initialize by putting all of the elements in the final partition.
  *identifiers = numberOfPartitions - 1;
  geom::SemiOpenInterval<1, int> indexRange(costs.lbounds(), costs.ubounds());
  partitionRegularGridWithBspTree(costs, identifiers,
				  ads::computeSum(costs),
				  indexRange, 0, numberOfPartitions);
  
}



//---------------------------------------------------------------------------
// 2-D
//---------------------------------------------------------------------------

inline
double
predictBestSplitting(const int n, double x, double y, double remains) {
  if (n == 1) {
    // No communication costs.
    return 0;
  }
  if (y > x) {
    std::swap(x, y);
  }
  if (y > remains) {
    return 2 * y;
  }
  double value;
  double minValue = std::numeric_limits<double>::max();
  for (int m = n / 2; m != 0; --m) {
    // The left part is no more expensive than the right, so the left can
    // use up to half the remaining value.
    value = y + predictBestSplitting(m, m * x / n, y, (remains - y) / 2);
    value += predictBestSplitting(n - m, (n - m) * x / n, y, remains - value);
    if (value < minValue) {
      minValue = value;
    }
  }
  return minValue;
}

inline
double
predictBestSplitting(int* splittingIndex, const int n, double x, double y) {
  assert(n > 1);
  if (y > x) {
    std::swap(x, y);
  }
  double value;
  double minValue = std::numeric_limits<double>::max();
  for (int m = n / 2; m != 0; --m) {
    value = y + predictBestSplitting(m, m * x / n, y, (minValue - y) / 2);
    value += predictBestSplitting(n - m, (n - m) * x / n, y, minValue - value);
    if (value < minValue) {
      minValue = value;
      *splittingIndex = m;
    }
  }
  return minValue;
}


inline
int
predictBestSplitting(const int n, const double x, const double y) {
  int splittingIndex = 0;
  predictBestSplitting(&splittingIndex, n, x, y);
  return splittingIndex;
}


// Reduce the index range to tightly contain the specified value.
inline
void
tightenIndexRange(const ads::Array<2, int>& identifiers, 
		  geom::SemiOpenInterval<2, int>* indexRange, 
		  const int value) {
  ads::FixedArray<2, int> index;
  bool hasValue;

  // For each dimension.
  for (int d = 0; d != 2; ++d) {
    const int i = d;
    const int j = (d + 1) % 2;
    // Lower.
    do {
      hasValue = false;
      index[i] = indexRange->getLowerCorner()[i];
      for (index[j] = indexRange->getLowerCorner()[j]; 
	   index[j] != indexRange->getUpperCorner()[j] && ! hasValue; 
	   ++index[j]) {
	if (identifiers(index) == value) {
	  hasValue = true;
	}
      }
      if (! hasValue) {
	indexRange->setLowerCoordinate(i, index[i] + 1);
      }
      if (indexRange->isEmpty()) {
	return;
      }
    } while (! hasValue);

    // Upper.
    do {
      hasValue = false;
      index[i] = indexRange->getUpperCorner()[i] - 1;
      for (index[j] = indexRange->getLowerCorner()[j]; 
	   index[j] != indexRange->getUpperCorner()[j] && ! hasValue; 
	   ++index[j]) {
	if (identifiers(index) == value) {
	  hasValue = true;
	}
      }
      if (! hasValue) {
	indexRange->setUpperCoordinate(i, index[i]);
      }
      if (indexRange->isEmpty()) {
	return;
      }
    } while (! hasValue);
  }
}




template<typename _T>
inline
void
partitionRegularGridWithBspTree
(const ads::Array<2, _T>& costs, ads::Array<2, int>* identifiers,
 const _T totalCost, geom::SemiOpenInterval<2, int> indexRange,
 const int identifiersBegin, const int identifiersEnd,
 const int predictionThreshhold) {
  // If we are down to a single piece or there are no elements.
  if (identifiersEnd - identifiersBegin <= 1 || indexRange.isEmpty()) {
    // Do nothing.
    return;
  }
  
  // See if we can tighten the index range.
  tightenIndexRange(*identifiers, &indexRange, identifiersEnd - 1);
  // If there are one or fewer elements.
  if (indexRange.computeContent() <= 1) {
    return;
  }

  // There are at least two elements to partition.

  // Choose the splitting dimension.
  const ads::FixedArray<2, int> extents = indexRange.getUpperCorner() - 
    indexRange.getLowerCorner();
  const int i = extents.max_index();
  const int j = (i + 1) % 2;

  // CONTINUE:
  //const int identifiersMiddle = (identifiersBegin + identifiersEnd) / 2;
  int identifiersMiddle;
  if (identifiersEnd - identifiersBegin <= predictionThreshhold) {
    identifiersMiddle = identifiersBegin + 
      predictBestSplitting(identifiersEnd - identifiersBegin,
			   extents[0], extents[1]);
  }
  else {
    identifiersMiddle = (identifiersBegin + identifiersEnd) / 2;
  }

  const _T splittingCost = totalCost * (identifiersMiddle - identifiersBegin) /
    _T(identifiersEnd - identifiersBegin);

  // Split.
  int splittingIndex = indexRange.getUpperCorner()[i];
  _T cost = 0;
  bool stop = false;
  ads::FixedArray<2, int> index;
  const ads::FixedArray<2, int>
    begin(indexRange.getLowerCorner()[j], indexRange.getUpperCorner()[j] - 1),
    end(indexRange.getUpperCorner()[j], indexRange.getLowerCorner()[j] - 1),
    stride(1, -1);
  bool flip = ((*identifiers)(indexRange.getLowerCorner()) != 
	       identifiersEnd - 1);
  for (index[i] = indexRange.getLowerCorner()[i]; 
       index[i] != indexRange.getUpperCorner()[i] && ! stop; 
       ++index[i]) {
    flip = ! flip;
    for (index[j] = begin[flip]; index[j] != end[flip] && ! stop; 
	 index[j] += stride[flip]) {
      // If this is an element that we are partitioning.
      if ((*identifiers)(index) == identifiersEnd - 1) {
	// If we should put this element in the lower partition.
	if (cost + costs(index) / 2 < splittingCost) {
	  // Place the element in the lower partition.
	  (*identifiers)(index) = identifiersMiddle - 1;
	  cost += costs(index);
	}
	else {
	  // Record the splitting index.
	  splittingIndex = index[i];
	  // Stop partitioning.  Break out of the loops.
	  stop = true;
	}
      }
    }
  }
  assert(splittingIndex != indexRange.getUpperCorner()[i]);

  // Recurse.
  const int upper = indexRange.getUpperCorner()[i];
  indexRange.setUpperCoordinate(i, splittingIndex + 1);
  partitionRegularGridWithBspTree(costs, identifiers, cost, indexRange,
				  identifiersBegin, identifiersMiddle,
				  predictionThreshhold);
  indexRange.setLowerCoordinate(i, splittingIndex);
  indexRange.setUpperCoordinate(i, upper);
  partitionRegularGridWithBspTree(costs, identifiers, totalCost - cost, 
				  indexRange, identifiersMiddle,
				  identifiersEnd, predictionThreshhold);
}


#if 0
template<typename _T>
inline
void
partitionRegularGridWithBspTree
(const ads::Array<2, _T>& costs, ads::Array<2, int>* identifiers,
 const _T totalCost, geom::SemiOpenInterval<2, int> indexRange,
 int identifiersBegin,  int identifiersEnd) {
  // If we are down to a single piece or there are no elements.
  if (identifiersEnd - identifiersBegin <= 1 || indexRange.isEmpty()) {
    // Do nothing.
    return;
  }

  const int identifiersMiddle = (identifiersBegin + identifiersEnd) / 2;
  const _T splittingCost = totalCost * (identifiersMiddle - identifiersBegin) /
    _T(identifiersEnd - identifiersBegin);

  // See if we can tighten the index range.
  tightenIndexRange(*identifiers, &indexRange, identifiersEnd - 1);
  // If there are one or fewer elements.
  if (indexRange.computeContent() <= 1) {
    return;
  }

  // There are at least two elements to partition.

  // Choose the splitting dimension.
  const ads::FixedArray<2, int> extents = indexRange.getUpperCorner() - 
    indexRange.getLowerCorner();
  const int i = extents.max_index();
  const int j = (i + 1) % 2;

  // Split.
  int splittingIndex = indexRange.getUpperCorner()[i];
  _T cost = 0;
  bool stop = false;
  ads::FixedArray<2, int> index;
  for (index[i] = indexRange.getLowerCorner()[i]; 
       index[i] != indexRange.getUpperCorner()[i] && ! stop; 
       ++index[i]) {
    for (index[j] = indexRange.getLowerCorner()[j]; 
	 index[j] != indexRange.getUpperCorner()[j] && ! stop; 
	 ++index[j]) {
      // If this is an element that we are partitioning.
      if ((*identifiers)(index) == identifiersEnd - 1) {
	// If we should put this element in the lower partition.
	if (cost + costs(index) / 2 < splittingCost) {
	  // Place the element in the lower partition.
	  (*identifiers)(index) = identifiersMiddle - 1;
	  cost += costs(index);
	}
	else {
	  // Record the splitting index.
	  splittingIndex = index[i];
	  // Stop partitioning.  Break out of the loops.
	  stop = true;
	}
      }
    }
  }
  assert(splittingIndex != indexRange.getUpperCorner()[i]);

  // Recurse.
  const int upper = indexRange.getUpperCorner()[i];
  indexRange.setUpperCoordinate(i, splittingIndex + 1);
  partitionRegularGridWithBspTree(costs, identifiers, cost, indexRange,
				  identifiersBegin, identifiersMiddle);
  indexRange.setLowerCoordinate(i, splittingIndex);
  indexRange.setUpperCoordinate(i, upper);
  partitionRegularGridWithBspTree(costs, identifiers, totalCost - cost, 
				  indexRange, identifiersMiddle,
				  identifiersEnd);
}
#endif

template<typename _T>
inline
void
partitionRegularGridWithBspTree(const ads::Array<2, _T>& costs,
				ads::Array<2, int>* identifiers,
				const int numberOfPartitions,
				const int predictionThreshhold) {
  assert(costs.ranges() == identifiers->ranges());
  // Initialize by putting all of the elements in the final partition.
  *identifiers = numberOfPartitions - 1;
  geom::SemiOpenInterval<2, int> indexRange(costs.lbounds(), costs.ubounds());
  partitionRegularGridWithBspTree(costs, identifiers,
				  ads::computeSum(costs),
				  indexRange, 0, numberOfPartitions,
				  predictionThreshhold);
}



//---------------------------------------------------------------------------
// 3-D
//---------------------------------------------------------------------------


// Reduce the index range to tightly contain the specified value.
inline
void
tightenIndexRange(const ads::Array<3, int>& identifiers, 
		  geom::SemiOpenInterval<3, int>* indexRange, 
		  const int value) {
  ads::FixedArray<3, int> index;
  bool hasValue;

  // For each dimension.
  for (int d = 0; d != 3; ++d) {
    const int i = d;
    const int j = (d + 1) % 3;
    const int k = (d + 2) % 3;
    // Lower.
    do {
      hasValue = false;
      index[i] = indexRange->getLowerCorner()[i];
      for (index[j] = indexRange->getLowerCorner()[j]; 
	   index[j] != indexRange->getUpperCorner()[j] && ! hasValue; 
	   ++index[j]) {
	for (index[k] = indexRange->getLowerCorner()[k]; 
	     index[k] != indexRange->getUpperCorner()[k] && ! hasValue;
	     ++index[k]) {
	  if (identifiers(index) == value) {
	    hasValue = true;
	  }
	}
      }
      if (! hasValue) {
	indexRange->setLowerCoordinate(i, index[i] + 1);
      }
      if (indexRange->isEmpty()) {
	return;
      }
    } while (! hasValue);

    // Upper.
    do {
      hasValue = false;
      index[i] = indexRange->getUpperCorner()[i] - 1;
      for (index[j] = indexRange->getLowerCorner()[j]; 
	   index[j] != indexRange->getUpperCorner()[j] && ! hasValue; 
	   ++index[j]) {
	for (index[k] = indexRange->getLowerCorner()[k]; 
	     index[k] != indexRange->getUpperCorner()[k] && ! hasValue;
	     ++index[k]) {
	  if (identifiers(index) == value) {
	    hasValue = true;
	  }
	}
      }
      if (! hasValue) {
	indexRange->setUpperCoordinate(i, index[i]);
      }
      if (indexRange->isEmpty()) {
	return;
      }
    } while (! hasValue);
  }
}


template<typename _T>
inline
void
partitionRegularGridWithBspTree
(const ads::Array<3, _T>& costs, ads::Array<3, int>* identifiers,
 const _T totalCost, geom::SemiOpenInterval<3, int> indexRange,
 int identifiersBegin,  int identifiersEnd) {
  // If we are down to a single piece or there are no elements.
  if (identifiersEnd - identifiersBegin <= 1 || indexRange.isEmpty()) {
    // Do nothing.
    return;
  }

  const int identifiersMiddle = (identifiersBegin + identifiersEnd) / 2;
  const _T splittingCost = totalCost * (identifiersMiddle - identifiersBegin) /
    _T(identifiersEnd - identifiersBegin);

  // See if we can tighten the index range.
  tightenIndexRange(*identifiers, &indexRange, identifiersEnd - 1);
  // If there are one or fewer elements.
  if (indexRange.computeContent() <= 1) {
    return;
  }

  // There are at least two elements to partition.

  // Choose the splitting dimension.  Assign i the index with the longest 
  // extent.
  const ads::FixedArray<3, int> extents = indexRange.getUpperCorner() - 
    indexRange.getLowerCorner();
  const int i = extents.max_index();
  const int j = (i + 1) % 3;
  const int k = (i + 2) % 3;
#if 0
  // CONTINUE: Does this help?
  // Assign k the index with the shortest extent.
  int j = (i + 1) % 3;
  int k = (i + 2) % 3;
  if (extents[j] < extents[k]) {
    std::swap(j, k);
  }
#endif

  // Split.
  int splittingIndex = indexRange.getUpperCorner()[i];
  _T cost = 0;
  bool stop = false;
  ads::FixedArray<3, int> index;
  for (index[i] = indexRange.getLowerCorner()[i]; 
       index[i] != indexRange.getUpperCorner()[i] && ! stop; 
       ++index[i]) {
    for (index[j] = indexRange.getLowerCorner()[j]; 
	 index[j] != indexRange.getUpperCorner()[j] && ! stop; 
	 ++index[j]) {
      for (index[k] = indexRange.getLowerCorner()[k]; 
	   index[k] != indexRange.getUpperCorner()[k] && ! stop;
	   ++index[k]) {
	// If this is an element that we are partitioning.
	if ((*identifiers)(index) == identifiersEnd - 1) {
	  // If we should put this element in the lower partition.
	  if (cost + costs(index) / 2 < splittingCost) {
	    // Place the element in the lower partition.
	    (*identifiers)(index) = identifiersMiddle - 1;
	    cost += costs(index);
	  }
	  else {
	    // Record the splitting index.
	    splittingIndex = index[i];
	    // Stop partitioning.  Break out of the loops.
	    stop = true;
	  }
	}
      }
    }
  }
  assert(splittingIndex != indexRange.getUpperCorner()[i]);

  // Recurse.
  const int upper = indexRange.getUpperCorner()[i];
  indexRange.setUpperCoordinate(i, splittingIndex + 1);
  partitionRegularGridWithBspTree(costs, identifiers, cost, indexRange,
				  identifiersBegin, identifiersMiddle);
  indexRange.setLowerCoordinate(i, splittingIndex);
  indexRange.setUpperCoordinate(i, upper);
  partitionRegularGridWithBspTree(costs, identifiers, totalCost - cost, 
				  indexRange, identifiersMiddle,
				  identifiersEnd);
}



template<typename _T>
inline
void
partitionRegularGridWithBspTree(const ads::Array<3, _T>& costs,
				ads::Array<3, int>* identifiers,
				const int numberOfPartitions) {
  assert(costs.ranges() == identifiers->ranges());
  // Initialize by putting all of the elements in the final partition.
  *identifiers = numberOfPartitions - 1;
  geom::SemiOpenInterval<3, int> indexRange(costs.lbounds(), costs.ubounds());
  partitionRegularGridWithBspTree(costs, identifiers,
				  ads::computeSum(costs),
				  indexRange, 0, numberOfPartitions);
}

END_NAMESPACE_CONCURRENT
