// -*- C++ -*-

#if !defined(__amr_writers_ipp__)
#error This file is an implementation detail.
#endif

BEGIN_NAMESPACE_AMR

template<typename _T, std::size_t _N>
inline
void
writeElements(std::ostream& out, const array::Array<_T, _N>& vector) {
  out << vector[0];
  for (std::size_t n = 1; n != _N; ++n) {
    out << " " << vector[n];
  }
}

template<typename _T, std::size_t _N>
inline
void
writeElementsInterlaced(std::ostream& out, const array::Array<_T, _N>& a,
			const array::Array<_T, _N>& b) {
  out << a[0] << " " << b[0];
  for (std::size_t n = 1; n != _N; ++n) {
    out << " " << a[n] << " " << b[n];
  }
}

template<typename _Patch, class _Traits>
inline
void
writeCellDataVtkXml(std::ostream& out,
		    const Orthtree<_Patch, _Traits>& orthtree,
		    typename Orthtree<_Patch, _Traits>::const_iterator node,
		    const PatchDescriptor<_Traits>& patchDescriptor) {
  typedef typename _Traits::SpatialIndex SpatialIndex;
  typedef typename _Traits::IndexList IndexList;

  typedef Orthtree<_Patch, _Traits> Orthtree;
  typedef typename Orthtree::Point Point;
  typedef typename Orthtree::const_iterator const_iterator;

  typedef typename _Patch::PatchData PatchData;
  typedef typename PatchData::ArrayConstView ArrayConstView;

  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"ImageData\">\n";

  // Exclude the ghost cells.
  ArrayConstView array(node->second.getPatchData().getInteriorArray());

  // The spatial index.
  const SpatialIndex& spatialIndex = node->first;
  // Begin RectilinearGrid
  out << "<ImageData WholeExtent=\"";
  writeElementsInterlaced(out,
			  Point(spatialIndex.getCoordinates() * 
				array.extents()),
			  Point((spatialIndex.getCoordinates() + 1) *
				array.extents()));
  out << "\" Origin=\"";
  writeElements(out, orthtree.getLowerCorner());
  out << "\" Spacing=\"";
  writeElements(out, Point(orthtree.getExtents(spatialIndex) / 
			   array.extents()));
  out << "\">\n";

  // Begin Piece.
  out << "<Piece Extent=\"";
  writeElementsInterlaced(out,
			  Point(spatialIndex.getCoordinates() * 
				array.extents()),
			  Point((spatialIndex.getCoordinates() + 1) *
				array.extents()));
  out << "\">\n";

  // PointData.
  out << "<PointData></PointData>\n";

  // Begin CellData.
  out << "<CellData>\n";

  const FieldDescriptor& field = patchDescriptor.getFieldDescriptor();
  // Begin DataArray.
  out << "<DataArray type=\"Float64\" Name=\"" << field.getName()
      << "\" NumberOfComponents=\"" << field.getNumberOfComponents()
      << "\" format=\"ascii\">\n";
  // Loop over the array.
  IndexList begin = array.bases();
  IndexList end = array.bases() + array.extents();
  IndexList i;
  for (i[2] = begin[2]; i[2] != end[2]; ++i[2]) {
    for (i[1] = begin[1]; i[1] != end[1]; ++i[1]) {
      for (i[0] = begin[0]; i[0] != end[0]; ++i[0]) {
	// Write the components of this field on a line.
	for (int n = 0; n != field.getNumberOfComponents(); ++n) {
	  out << array(i)[n] << " ";
	}
	out << '\n';
      }
    }
  }
  // End DataArray.
  out << "</DataArray>\n";

  // End CellData.
  out << "</CellData>\n";

  // End Piece.
  out << "</Piece>\n";
  out << "</ImageData>\n";
  out << "</VTKFile>\n";
}

template<typename _Patch, class _Traits>
inline
void
writeCellDataParaview(const std::string& name,
		      const Orthtree<_Patch, _Traits>& orthtree,
		      const PatchDescriptor<_Traits>& patchDescriptor) {
  typedef Orthtree<_Patch, _Traits> Orthtree;
  typedef typename Orthtree::const_iterator const_iterator;

  // Open the ParaView file.
  std::string paraviewName = name;
  paraviewName += ".pvd";
  std::ofstream paraviewFile(paraviewName.c_str());
  paraviewFile << "<?xml version=\"1.0\"?>\n"
	       << "<VTKFile type=\"Collection\">\n"
	       << "<Collection>\n";
  // For each node.
  for (const_iterator node = orthtree.begin(); node != orthtree.end(); ++node) {
    std::string vtkName = name;
    std::ostringstream oss;
    oss << node->first.getCode();
    vtkName += oss.str();
    vtkName += ".vti";

    paraviewFile << "<DataSet part=\"" << int(node->first.getLevel())
		 << "\" file=\"" << vtkName << "\"/>\n";

    std::ofstream vtkFile(vtkName.c_str());
    writeCellDataVtkXml(vtkFile, orthtree, node, patchDescriptor);
  }
  paraviewFile << "</Collection>\n";
  paraviewFile << "</VTKFile>\n";
}

END_NAMESPACE_AMR
