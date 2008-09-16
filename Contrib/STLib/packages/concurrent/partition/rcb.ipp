// -*- C++ -*-

#if !defined(__partition_rcb_ipp__)
#error This file is an implementation detail of rcb.
#endif

BEGIN_NAMESPACE_CONCURRENT

template<int N, typename id_type, typename number_type>
inline
void
rcb_split(id_type** partition_begin,id_type** partition_end,
	  const ads::FixedArray<N,number_type>** records_begin,
	  const ads::FixedArray<N,number_type>** records_end) {
  // If there is only a single processor, there is no need to split.
  const int num_processors = partition_end - partition_begin;
  if (num_processors == 1) {
    return;
  }

  //
  // Determine the number of processors and records in the left and 
  // right branches.
  //
  const int num_processors_left = num_processors / 2;
  id_type** partition_mid = partition_begin + num_processors_left;
  const int num_records = records_end - records_begin;
  const int num_records_left =
    int(number_type(num_records) * 
	number_type(num_processors_left) / 
	number_type(num_processors));
  *partition_mid = *partition_begin + num_records_left;
  
  //
  // Determine the splitting dimension.
  //
  geom::BBox<N,number_type> bbox;
  bbox.bound(ads::constructIndirectIterator(records_begin),
	     ads::constructIndirectIterator(records_end));
  ads::FixedArray<N,number_type> extents = bbox.getUpperCorner() - 
    bbox.getLowerCorner();
  const int split_dimension = extents.max_index();
  
  //
  // The comparison functor for splitting the records.
  //
  // Dereference the pointer.
  typedef ads::Dereference< const ads::FixedArray<N>* > Deref;
  Deref deref;
  // Index the point.
  typedef ads::IndexConstObject< ads::FixedArray<N> > Index;
  Index ind;
  // Index in the splitting dimension.
  typedef std::binder2nd< Index > IndexSD;
  IndexSD indsd(ind, split_dimension);
  // Dereference; then index in the splitting dimension.
  typedef ads::unary_compose_unary_unary< IndexSD, Deref > ISDD;
  ISDD isdd(indsd, deref);
  // Compare two handles to points by their splitting dimension coordinate.
  typedef ads::binary_compose_binary_unary< std::less<double>, ISDD, ISDD > 
    Comp;
  std::less<double> less_than;
  Comp comp(less_than, isdd, isdd);
  
  //
  // Split the records.
  //
  const ads::FixedArray<N,number_type>** records_mid = 
    records_begin + num_records_left;
  std::nth_element(records_begin, records_mid, records_end, comp);
  
  //
  // Recurse.
  //
  rcb_split(partition_begin, partition_mid, records_begin, records_mid);
  rcb_split(partition_mid, partition_end, records_mid, records_end);
}


template<int N, typename id_type, typename number_type>
inline
void
rcb(const int num_processors, const int num_records,
    id_type* identifiers, id_type** id_partition,
    const number_type* positions) {
  typedef ads::FixedArray<N,number_type> Point;

  assert(num_processors >= 1);

  // Assign the first and last partition pointers.
  id_partition[0] = identifiers;
  id_partition[num_processors] = identifiers + num_records;

  // Make an array of pointers to the positions.
  ads::Array<1,const Point*> pos_ptrs(num_records);
  {
    const Point* p = reinterpret_cast< const Point* >(positions);
    typename ads::Array<1,const Point*>::iterator i = pos_ptrs.begin();
    const typename ads::Array<1,const Point*>::iterator i_end = pos_ptrs.end();
    while (i != i_end) {
      *i++ = p++;
    }
  }

  // Determine the partition.
  rcb_split(id_partition, id_partition + num_processors,
	     pos_ptrs.begin(), pos_ptrs.end());

  // Re-arrange the identifiers.
  ads::Array<1,id_type> tmp(identifiers, identifiers + num_records);
  const Point** p = pos_ptrs.data();
  const Point* pos = reinterpret_cast< const Point* >(positions);
  id_type* i = identifiers;
  id_type* i_end = identifiers + num_records;
  while (i != i_end) {
    *i++ = tmp[*p++ - pos];
  }
}

template<int N, typename id_type, typename number_type>
inline
void
rcb_allocate(const int num_processors, const int num_records,
	     id_type*& identifiers, id_type**& id_partition, 
	     number_type*& positions) {
  identifiers = new id_type[num_records];
  id_partition = new id_type*[num_processors + 1];
  positions = new number_type[N * num_records];
}

template<typename id_type, typename number_type>
inline
void
rcb_deallocate(id_type*& identifiers, id_type**& id_partition, 
	       number_type*& positions) {
  delete[] identifiers;
  delete[] id_partition;
  delete[] positions;
}

END_NAMESPACE_CONCURRENT
