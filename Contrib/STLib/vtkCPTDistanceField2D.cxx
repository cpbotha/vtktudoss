// Stef Busking, 080530
// VTK wrapper class for Mauch's STLib CPT

#include "vtkCPTDistanceField2D.h"

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <cpt.h>

#include <cassert>


vtkStandardNewMacro(vtkCPTDistanceField2D);

//----------------------------------------------------------------------------
// Constructor
vtkCPTDistanceField2D::vtkCPTDistanceField2D()
 : MaximumDistance(1.0),
   Padding(0.0),
   UnusedAxis(2)
{
  Dimensions[0] = Dimensions[1] = 64;
  Dimensions[2] = 1;
  Domain[0] = Domain[1] = Domain[2] = Domain[3] = 0;
}

//----------------------------------------------------------------------------
// Destructor
vtkCPTDistanceField2D::~vtkCPTDistanceField2D()
{
}

//----------------------------------------------------------------------------
// PrintSelf
void vtkCPTDistanceField2D::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Maximum Distance: " << this->MaximumDistance << "\n";
  os << indent << "Dimensions: (" << this->Dimensions[0] << ", "
     << this->Dimensions[1] << ", "
     << this->Dimensions[2] << ")\n";
  os << indent << "Padding: " << this->Padding << "\n";
}

//----------------------------------------------------------------------------
// RequestInformation
int vtkCPTDistanceField2D::RequestInformation(vtkInformation *vtkNotUsed(request),
                                            vtkInformationVector **inputVector,
                                            vtkInformationVector *outputVector)
{
  // Get info about the first input and output
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Set scalar type
  vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_DOUBLE, 1);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               0, this->Dimensions[0]-1,
               0, this->Dimensions[1]-1,
               0, this->Dimensions[2]-1);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!input)
  {
    vtkErrorMacro(<<"Input not present or incorrect type!");
    return 1;
  }

  // TODO: Figure out how to get an updated source here so we can get the
  // correct bounds since we can't call input->Update() anymore since VTK 6.
  // To handle this for now, the source needs to be updated manually or data
  // must be passed in with SetInputData().

  // Get padded bounds
  double bounds[6];
  input->GetBounds(bounds);
  if (!PadBounds(bounds))
  {
    // The bounds were invalid
    return 1;
  }

  // Set origin and spacing for the distance field
  double origin[3];
  double spacing[3];
  for (int i = 0; i < 3; ++i)
  {
    origin[i] = bounds[2 * i];
    spacing[i] = (bounds[2 * i + 1] - bounds[2 * i]) / (Dimensions[i] - 1);
  }
  outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
  outInfo->Set(vtkDataObject::SPACING(), spacing, 3);

  return 1;
}

//----------------------------------------------------------------------------
// RequestData
int vtkCPTDistanceField2D::RequestData(vtkInformation *vtkNotUsed(request),
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector)
{
  // Get the information objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and output data objects
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *output = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDebugMacro(<<"Starting CPT filter...");

  // Check validity of network
  if (!input || !output)
  {
    vtkErrorMacro(<<"Input / output not present or incorrect type!");
    return 1;
  }
  vtkIdType numPoints = input->GetNumberOfPoints();
  vtkIdType numCells = input->GetNumberOfCells();
  if (numPoints < 1 || numCells < 1)
  {
    vtkErrorMacro(<<"No mesh data!");
    return 1;
  }

  // Check validity of parameters
  if (MaximumDistance < 0)
  {
    vtkErrorMacro(<<"MaximumDistance can not be negative!");
    return 1;
  }
  if (Dimensions[0] < 0 || Dimensions[1] < 0 || Dimensions[2] < 0)
  {
    vtkErrorMacro(<<"Invalid Dimensions!");
    return 1;
  }
  // Calculate bounds
  double bounds[6];
  input->GetBounds(bounds);
  if (!PadBounds(bounds))
  {
    // The bounds were invalid
    return 1;
  }
  // CPT domain uses a different ordering...
  double domain[4];
  if (Domain[0] != 0 || Domain[1] != 0 || Domain[2] != 0 || Domain[3] != 0) {
    domain[0] = Domain[0];
    domain[1] = Domain[1];
    domain[2] = Domain[2];
    domain[3] = Domain[3];
  }
  else {
    switch (UnusedAxis) {
    case 0:
      domain[0] = bounds[2];
      domain[1] = bounds[4];
      domain[2] = bounds[3];
      domain[3] = bounds[5];
      break;
    case 1:
      domain[0] = bounds[0];
      domain[1] = bounds[4];
      domain[2] = bounds[1];
      domain[3] = bounds[5];
      break;
    case 2:
      domain[0] = bounds[0];
      domain[1] = bounds[2];
      domain[2] = bounds[1];
      domain[3] = bounds[3];
      break;
    }
  }

  vtkDebugMacro(<<"Converting to BRep...");

  // Triangulate the input mesh first
  vtkSmartPointer<vtkTriangleFilter> triangulate =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triangulate->SetInputData(input);
  triangulate->Update();
  vtkPolyData *mesh = triangulate->GetOutput();

  // Extract the BRep from the mesh
  vtkDebugMacro(<<"Collecting points...");
  vtkIdType numVerts = mesh->GetNumberOfPoints();
  double *verts = new double[2 * numVerts];
  int axis_indices[2];
  switch (UnusedAxis) {
    case 0:
    // Drop the x value
    axis_indices[0] = 1;
    axis_indices[1] = 2;
    break;
    case 1:
    axis_indices[0] = 0;
    // Drop the y value
    axis_indices[1] = 2;
    break;
    case 2:
    axis_indices[0] = 0;
    axis_indices[1] = 1;
    // Drop the z value
    break;
  }
  for (vtkIdType i = 0; i < numVerts; ++i)
  {
    double v[3];
    // Copy the point
    mesh->GetPoint(i, v);
    verts[i * 2] = v[axis_indices[0]];
    verts[i * 2 + 1] = v[axis_indices[1]];
  }
  vtkDebugMacro(<<"Collecting faces...");
  vtkIdType numFaces = mesh->GetNumberOfCells();
  unsigned int *faces = new unsigned int[2 * numFaces];
  for (vtkIdType i = 0; i < numFaces; ++i)
  {
    // Get the cell
    vtkCell *cell = mesh->GetCell(i);
    // Ignore non-line cells
    if (cell->GetCellType() == VTK_LINE)
    {
      // Get and copy the cell points
      vtkIdList *cellPoints = cell->GetPointIds();
      assert(cellPoints->GetNumberOfIds() == 2);
      for (int j = 0; j < 2; ++j)
      {
        faces[i * 2 + j] = cellPoints->GetId(j);
      }
    }
  }

  vtkDebugMacro(<<"Preparing the CPT module...");

  // Allocate the field
  // VTK6+ prototype: vtkDataObject, vtkInformation
  AllocateOutputData(output, outInfo);

  // Set up the CPT module
  cpt::State<2, double> state;
  state.setParameters(domain, this->MaximumDistance);
  state.setLattice(this->Dimensions, domain);
  // We use a single grid covering the entire lattice
  int zeroes[2] = {0, 0};
  state.insertGrid(zeroes, this->Dimensions,
                   static_cast<double*>(output->GetScalarPointer()), 0, 0, 0);
  //state.setBRepWithNoClipping(numVerts, verts, numFaces, faces);
  state.setBRep(numVerts, verts, numFaces, faces);

  // Now compute the distance field
  vtkDebugMacro(<<"Computing CPT...");
  state.computeClosestPointTransform();

  //state.floodFillDetermineSign(this->MaximumDistance);
  state.floodFillAtBoundary(this->MaximumDistance);

  // Check the grids (see cpt3.cc example)
  //if (!state.areGridsValid())
  //{
  //	vtkErrorMacro(<<"Warning! Computed distance field is not valid!");
  // Not sure if we should exit here...
  //}

  vtkDebugMacro(<<"Done!");

  // Clean up
  delete [] verts;
  delete [] faces;

  return 1;
}

//----------------------------------------------------------------------------
// FillInputPortInformation
int vtkCPTDistanceField2D::FillInputPortInformation(int vtkNotUsed(port),
                                                  vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

//----------------------------------------------------------------------------
// PadBounds
bool vtkCPTDistanceField2D::PadBounds(double *bounds)
{
  if (bounds[0] > bounds[1] || bounds[2] > bounds[3] || bounds[4] > bounds[5])
  {
    // TODO: this seems to occur sometimes, e.g., when a reader has not yet updated
    vtkErrorMacro(<<"Input polydata has invalid bounds: (" << bounds[0] << ", " << bounds[2] << ", " << bounds[4]<<
                  ") to (" << bounds[1] << ", " << bounds[3] << ", " << bounds[5] << ")");
    return false;
  }
  for (int i = 0; i < 3; ++i)
  {
    bounds[2 * i] -= Padding;
    bounds[2 * i + 1] += Padding;
  }
  vtkDebugMacro(<<"Bounds after padding: (" << bounds[0] << ", " << bounds[2] << ", " << bounds[4]<<
                ") to (" << bounds[1] << ", " << bounds[3] << ", " << bounds[5] << ")");

  return true;
}
