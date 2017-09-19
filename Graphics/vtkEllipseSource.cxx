/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkEllipseSource.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkEllipseSource.h"

#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"



vtkCxxRevisionMacro(vtkEllipseSource, "$Revision: 1.7 $");
vtkStandardNewMacro(vtkEllipseSource);

vtkEllipseSource::vtkEllipseSource()
{
  this->SemiMajorAxisLength = 1.0;
  this->SemiMinorAxisLength = 0.5;
  this->NumberOfSteps = 16;
  
  this->SetNumberOfInputPorts(0);
}

int vtkEllipseSource::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the ouptut
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  double interval = 2 * vtkMath::Pi() / this->NumberOfSteps;
  double alpha, sinalpha, cosalpha;
  double pt[3];
  pt[2] = 0.0;

  vtkPoints *new_points;
  new_points = vtkPoints::New();
  new_points->Allocate(this->NumberOfSteps);

  vtkCellArray *new_lines;
  new_lines = vtkCellArray::New();
  new_lines->Allocate(this->NumberOfSteps);

  vtkIdType prev_pt=0, cur_pt=0, pts[2];
  
  // calculate points
  for (int i = 0; i < this->NumberOfSteps; i++)
    {
    alpha = i * interval; // current angle
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);
/*
    pt[0] = this->SemiMajorAxisLength * cosalpha -
      this->SemiMinorAxisLength * sinalpha;
    pt[1] = this->SemiMajorAxisLength * cosalpha +
      this->SemiMinorAxisLength * sinalpha;
*/
    pt[0] = this->SemiMajorAxisLength * cosalpha;
    pt[1] = this->SemiMinorAxisLength * sinalpha;


    prev_pt = cur_pt;
    cur_pt = new_points->InsertNextPoint(pt);

    if (i > 0)
      {
      pts[0] = prev_pt;
      pts[1] = cur_pt;
      new_lines->InsertNextCell(2, pts);
      }
				
    }

	// connect up the last segment
	pts[0] = pts[1];
	pts[1] = 0;
	new_lines->InsertNextCell(2, pts);

  output->SetPoints(new_points);
  new_points->Delete();

  output->SetLines(new_lines);
  new_lines->Delete();

  return 1;
}

void vtkEllipseSource::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NumberOfSteps: " << this->NumberOfSteps << "\n";
  os << indent << "SemiMajorAxisLength: " << this->SemiMajorAxisLength << "\n";
  os << indent << "SemiMajorAxisLength: " << this->SemiMinorAxisLength << "\n";
}
