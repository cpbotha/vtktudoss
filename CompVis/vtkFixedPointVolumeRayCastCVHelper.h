/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFixedPointVolumeRayCastCVHelper.h,v $
  Language:  C++
  Date:      $Date: 2005-05-04 14:13:58 $
  Version:   $Revision: 1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkFixedPointVolumeRayCastCVHelper - A helper that generates composite images for the volume ray cast mapper
// .SECTION Description
//
// .SECTION see also
// vtkFixedPointVolumeRayCastMapper

#ifndef __vtkFixedPointVolumeRayCastCVHelper_h
#define __vtkFixedPointVolumeRayCastCVHelper_h

#include "vtkFixedPointVolumeRayCastHelper.h"

class vtkCVFixedPointVolumeRayCastMapper;
class vtkVolume;

//BTX  
#define VTKKWRCHelper_InterpolateSingleScalarComponent( VAL, IDX )               \
    {                                                                                   \
    VAL[IDX] =                                                                         \
    (0x7fff + ((A[IDX]*((0x4000 + w1Xw1Y*w1Z)>>VTKKW_FP_SHIFT)) +                      \
               (B[IDX]*((0x4000 + w2Xw1Y*w1Z)>>VTKKW_FP_SHIFT)) +                      \
               (C[IDX]*((0x4000 + w1Xw2Y*w1Z)>>VTKKW_FP_SHIFT)) +                      \
               (D[IDX]*((0x4000 + w2Xw2Y*w1Z)>>VTKKW_FP_SHIFT)) +                      \
               (E[IDX]*((0x4000 + w1Xw1Y*w2Z)>>VTKKW_FP_SHIFT)) +                      \
               (F[IDX]*((0x4000 + w2Xw1Y*w2Z)>>VTKKW_FP_SHIFT)) +                      \
               (G[IDX]*((0x4000 + w1Xw2Y*w2Z)>>VTKKW_FP_SHIFT)) +                      \
               (H[IDX]*((0x4000 + w2Xw2Y*w2Z)>>VTKKW_FP_SHIFT)))) >> VTKKW_FP_SHIFT;   \
    }                                                                                   \
//ETX 

class VTK_EXPORT vtkFixedPointVolumeRayCastCVHelper : public vtkFixedPointVolumeRayCastHelper
{
public:
  static vtkFixedPointVolumeRayCastCVHelper *New();
  vtkTypeRevisionMacro(vtkFixedPointVolumeRayCastCVHelper,vtkFixedPointVolumeRayCastHelper);
  void PrintSelf( ostream& os, vtkIndent indent );

  virtual void  GenerateImage( int threadID, 
                               int threadCount,
                               vtkVolume *vol,
                               vtkCVFixedPointVolumeRayCastMapper *mapper);

protected:
  vtkFixedPointVolumeRayCastCVHelper();
  ~vtkFixedPointVolumeRayCastCVHelper();

private:
  vtkFixedPointVolumeRayCastCVHelper(const vtkFixedPointVolumeRayCastCVHelper&);  // Not implemented.
  void operator=(const vtkFixedPointVolumeRayCastCVHelper&);  // Not implemented.
};

#endif


