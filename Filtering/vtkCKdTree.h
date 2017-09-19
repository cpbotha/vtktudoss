/*=========================================================================
vtkCKdTree was created by Peter Krekel and he cannot remember why he made it.
HOWEVER, it seems as if some of the functionality of the vtkKdTree
was not wrapped for python, so the C probably means something like Crappy.
=========================================================================*/

#include <iostream>
#include <cmath>
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkKdTree.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkPointData.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPolyDataAlgorithm.h"

#ifndef __vtkCollisionDetection_h
#define __vtkCollisionDetection_h

class vtkDataArray;
class vtkMatrix4x4;

class VTK_EXPORT vtkCKdTree : public vtkPolyDataAlgorithm
{
public:
	vtkTypeRevisionMacro(vtkCKdTree,vtkPolyDataAlgorithm);
	
	
	vtkGetObjectMacro(ptidtypearray, vtkIdTypeArray);
	vtkGetMacro(Boxsize, double);
	vtkSetMacro(Boxsize, double);
	static vtkCKdTree *New();

	
	void BuildKdTree();
	void QueryRegion(double point[3]);

private:
	vtkIdTypeArray *ptidtypearray;
	vtkPolyData *input_first;
	vtkKdTree *kdtree1;
	double Boxsize;
	

protected:
    vtkCKdTree(void);
	~vtkCKdTree(void);
};


#endif




