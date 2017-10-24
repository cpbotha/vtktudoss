#pragma once
// Stef Busking, 080530
// VTK wrapper class for Mauch's STLib CPT

#include <vtkImageAlgorithm.h>

class VTK_EXPORT vtkCPTDistanceField2D : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(vtkCPTDistanceField2D, vtkImageAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  static vtkCPTDistanceField2D *New();

  //int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

  // Maximum distance up to which the field should be computed
  vtkGetMacro(MaximumDistance, double);
  vtkSetMacro(MaximumDistance, double);

  // The distance field will have the polydata's bounds, padded by this amount
  vtkGetMacro(Padding, double);
  vtkSetMacro(Padding, double);

  // Resolution of the distance field
  // The first two dimensions are used no matter what the unused axis is
  vtkGetVector3Macro(Dimensions, int);
  vtkSetVector3Macro(Dimensions, int);

  // The index of the axis that isn't used
  // 0 - X, 1 - Y, 2 - Z
  vtkGetMacro(UnusedAxis, int);
  vtkSetMacro(UnusedAxis, int);

  // Manually set the domain of the distance field
  // Note: If this is not set, the bounds of the input
  // polydata is used to set the domain
  vtkSetVector4Macro(Domain, double);
  vtkGetVector4Macro(Domain, double);

protected:
  vtkCPTDistanceField2D();
  ~vtkCPTDistanceField2D();

  virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  bool PadBounds(double *bounds);

  double MaximumDistance;
  double Padding;
  int Dimensions[3];
  double Bounds[6];
  int UnusedAxis;
  double Domain[4];

private:
  vtkCPTDistanceField2D(const vtkCPTDistanceField2D&);        // Not implemented.
  void operator=(const vtkCPTDistanceField2D&);        // Not implemented.
};
