/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkArrowSource.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkEllipseSource - Draws an ellipse.
// .SECTION Description
// Draws an ellipse (by definition 2D) centred at 0,0,0; minor axis up
// to 0,1,0, major axis across to 1,0,0, so its a circle to start off
// with.  You can set the major and minor axis lengths, but you have to
// use a vtkTransformPolyDataFilter if you want to rotate and/or
// translate the ellipse.
//
// .SECTION Thanks
// Charl P. Botha c.p.botha@tudelft.nl for creating and contributing
// this class.
//
// .SECTION See Also
// vtkTransformPolyData

#ifndef __vtkEllipseSource_h
#define __vtkEllipseSource_h

#include "vtkPolyDataAlgorithm.h"

class VTK_EXPORT vtkEllipseSource : public vtkPolyDataAlgorithm
{
public:
  // Description
  // Construct default ellipse source.
  static vtkEllipseSource *New();

  vtkTypeRevisionMacro(vtkEllipseSource,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the length of the ellipse major and minor axes
  vtkSetMacro(SemiMajorAxisLength,double);
  vtkGetMacro(SemiMajorAxisLength,double);
  vtkSetMacro(SemiMinorAxisLength,double);
  vtkGetMacro(SemiMinorAxisLength,double);
    
  // Description:
  // Set the resolution of ellipse polyline.  This is the number of
  // points / lines on the ellipse arc.
  vtkSetClampMacro(NumberOfSteps,int,2,1024);
  vtkGetMacro(NumberOfSteps,int);

protected:
  vtkEllipseSource();
  ~vtkEllipseSource() {};

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  double SemiMajorAxisLength;
  double SemiMinorAxisLength;
  int NumberOfSteps;

private:
  vtkEllipseSource(const vtkEllipseSource&); // Not implemented.
  void operator=(const vtkEllipseSource&); // Not implemented.
};

#endif


