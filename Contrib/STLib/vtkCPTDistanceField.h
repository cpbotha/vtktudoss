#pragma once
// Stef Busking, 080530
// VTK wrapper class for Mauch's STLib CPT

#include <vtkImageAlgorithm.h>

class VTK_EXPORT vtkCPTDistanceField : public vtkImageAlgorithm
{
public:
	vtkTypeRevisionMacro(vtkCPTDistanceField, vtkImageAlgorithm);
	void PrintSelf(ostream &os, vtkIndent indent);

	static vtkCPTDistanceField *New();

	//int ProcessRequest(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

	// Maximum distance up to which the field should be computed
	vtkGetMacro(MaximumDistance, double);
	vtkSetMacro(MaximumDistance, double);

	// The distance field will have the polydata's bounds, padded by this amount
	vtkGetMacro(Padding, double);
	vtkSetMacro(Padding, double);

	// Resolution of the distance field
	vtkGetVector3Macro(Dimensions, int);
	vtkSetVector3Macro(Dimensions, int);

protected:
	vtkCPTDistanceField();
	~vtkCPTDistanceField();

	virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	virtual int FillInputPortInformation(int port, vtkInformation* info);

	bool PadBounds(double *bounds);

	double MaximumDistance;
	double Padding;
	int Dimensions[3];

private:
	vtkCPTDistanceField(const vtkCPTDistanceField&);  // Not implemented.
	void operator=(const vtkCPTDistanceField&);  // Not implemented.
};
