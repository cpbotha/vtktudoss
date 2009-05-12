/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFixedPointVolumeRayCastMapper.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCVFixedPointVolumeRayCastMapper - A fixed point mapper for volumes
// .SECTION Description


// .SECTION see also
// vtkVolumeMapper

#ifndef __vtkCVFixedPointVolumeRayCastMapper_h
#define __vtkCVFixedPointVolumeRayCastMapper_h

#include "vtkFixedPointVolumeRayCastMapper.h"

class vtkFixedPointVolumeRayCastCVHelper;

//BTX
// Forward declaration needed for use by friend declaration below.
VTK_THREAD_RETURN_TYPE CVFixedPointVolumeRayCastMapper_CastRays( void *arg );
//ETX

class VTK_EXPORT vtkCVFixedPointVolumeRayCastMapper : public vtkFixedPointVolumeRayCastMapper
{
public:
  static vtkCVFixedPointVolumeRayCastMapper *New();
  vtkTypeRevisionMacro(vtkCVFixedPointVolumeRayCastMapper,vtkFixedPointVolumeRayCastMapper);
  void PrintSelf( ostream& os, vtkIndent indent );


  void RenderSubVolume();

protected:
  vtkCVFixedPointVolumeRayCastMapper();
  ~vtkCVFixedPointVolumeRayCastMapper();

  
  friend VTK_THREAD_RETURN_TYPE CVFixedPointVolumeRayCastMapper_CastRays( void *arg );

  vtkFixedPointVolumeRayCastCVHelper *CVHelper;
  

private:
  vtkCVFixedPointVolumeRayCastMapper(const vtkCVFixedPointVolumeRayCastMapper&);  // Not implemented.
  void operator=(const vtkCVFixedPointVolumeRayCastMapper&);  // Not implemented.
};




#endif
