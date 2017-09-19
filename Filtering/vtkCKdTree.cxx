/*=========================================================================
vtkCKdTree was created by Peter Krekel and he cannot remember why he made it.
HOWEVER, it seems as if some of the functionality of the vtkKdTree
was not wrapped for python, so the C probably means something like Crappy.
=========================================================================*/
#include "vtkCKdTree.h"

vtkCxxRevisionMacro(vtkCKdTree, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkCKdTree);



// Constructor
vtkCKdTree::vtkCKdTree()
{
	this->Boxsize = 5;
}

// Destructor
vtkCKdTree::~vtkCKdTree()
{
}



void vtkCKdTree::BuildKdTree()
{
	this->kdtree1 = vtkKdTree::New();	
	//this->kdtree1->SetDataSet(input_first);
	//vtkInformation *in_first_Info = this->GetInput()->GetInformationObject(0);
	input_first = vtkPolyData::SafeDownCast(this->GetInput());
		//in_first_Info->Get(vtkDataObject::DATA_OBJECT()));
	//printf("%d", input_first->GetNumberOfPoints());
	this->kdtree1->BuildLocatorFromPoints(input_first->GetPoints());
}

void vtkCKdTree::QueryRegion(double point[3])
{
	//int *return_idlist = idlist;
	double area[6];
	area[0] = point[0] - this->Boxsize;
	area[1] = point[0] + this->Boxsize;
	area[2] = point[1] - this->Boxsize;
	area[3] = point[1] + this->Boxsize;
	area[4] = point[2] - this->Boxsize;
	area[5] = point[2] + this->Boxsize;
	this->ptidtypearray = vtkIdTypeArray::New();
	this->kdtree1->FindPointsInArea(area, this->ptidtypearray);
	//return_idlist = ptidtypearray->GetPointer(6);
	
}


