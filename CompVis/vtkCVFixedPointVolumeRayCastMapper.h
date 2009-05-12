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

//BTX
  // Description:
  // WARNING: INTERNAL METHOD - NOT INTENDED FOR GENERAL USE
  // Initialize rendering for this volume.
  // cpbotha: I've had to override this JUST so that my RenderSubVolume() is called,
  // as VTK doesn't have that as virtual.
  void Render( vtkRenderer *, vtkVolume * );
//ETX

  // Description:
  // Activate comparative visualisation mode.  When 0, this acts just like a normal
  // vtkFixedPointVRCM.  With CompVis 1, the first comparison mode is activated.
  vtkSetMacro(CompVisMode, int);
  vtkGetMacro(CompVisMode, int);

  vtkGetObjectMacro( CVHelper, vtkFixedPointVolumeRayCastCVHelper );

protected:
  vtkCVFixedPointVolumeRayCastMapper();
  ~vtkCVFixedPointVolumeRayCastMapper();

  
  friend VTK_THREAD_RETURN_TYPE CVFixedPointVolumeRayCastMapper_CastRays( void *arg );

  vtkFixedPointVolumeRayCastCVHelper *CVHelper;

  int CompVisMode;
  

private:
  vtkCVFixedPointVolumeRayCastMapper(const vtkCVFixedPointVolumeRayCastMapper&);  // Not implemented.
  void operator=(const vtkCVFixedPointVolumeRayCastMapper&);  // Not implemented.
};




#endif
