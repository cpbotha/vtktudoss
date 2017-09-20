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
// .NAME vtkCVFixedPointVolumeRayCastMapper - A fixed point mapper for the
// comparative visualisation of multi-component volumes.
// .SECTION Description
//
// I derived this class and also created the vtkFixedPointVolumeRayCastCVHelper
// to test ideas for multi-component comparative visualisation of volumes.  This
// specific experiment was not that successful, primarily due to speed issues
// that I should have expected. :)  It does work, but too slow: 128x128x80 x 3 
// volume gets rendered at about 2 to 3 seconds per frame on my core duo 1 2GHz
// laptop.  That's too slow to get for example animation going at any reasonable rate.
//
// I'm keeping this class in here, because it does show a clear example of how
// to stuff at the ray level with multiple TFs and multiple components.
//
// -- cpbotha

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
  vtkTypeMacro(vtkCVFixedPointVolumeRayCastMapper,vtkFixedPointVolumeRayCastMapper);
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
