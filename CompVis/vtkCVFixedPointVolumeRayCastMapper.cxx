/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFixedPointVolumeRayCastMapper.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCVFixedPointVolumeRayCastMapper.h"

#include "vtkCommand.h"

#include "vtkFixedPointVolumeRayCastCVHelper.h"

#include "vtkFixedPointVolumeRayCastCompositeGOHelper.h"
#include "vtkFixedPointVolumeRayCastCompositeGOShadeHelper.h"
#include "vtkFixedPointVolumeRayCastCompositeHelper.h"
#include "vtkFixedPointVolumeRayCastCompositeShadeHelper.h"
#include "vtkFixedPointVolumeRayCastMIPHelper.h"
#include "vtkMultiThreader.h"
#include "vtkObjectFactory.h"
#include <math.h>

vtkCxxRevisionMacro(vtkCVFixedPointVolumeRayCastMapper, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkCVFixedPointVolumeRayCastMapper); 



// Construct a new vtkFixedPointVolumeRayCastMapper with default values
vtkCVFixedPointVolumeRayCastMapper::vtkCVFixedPointVolumeRayCastMapper()
{
  this->CVHelper   = vtkFixedPointVolumeRayCastCVHelper::New();
}

// Destruct a vtkFixedPointVolumeRayCastMapper - clean up any memory used
vtkCVFixedPointVolumeRayCastMapper::~vtkCVFixedPointVolumeRayCastMapper()
{

  this->CVHelper->Delete();

}

// This is the render method for the subvolume
void vtkCVFixedPointVolumeRayCastMapper::RenderSubVolume()
{
  // Set the number of threads to use for ray casting,
  // then set the execution method and do it.
  this->InvokeEvent( vtkCommand::VolumeMapperRenderStartEvent, NULL );
  this->Threader->SetSingleMethod( CVFixedPointVolumeRayCastMapper_CastRays,
                                   (void *)this);
  this->Threader->SingleMethodExecute();
  this->InvokeEvent( vtkCommand::VolumeMapperRenderEndEvent, NULL );
}


VTK_THREAD_RETURN_TYPE CVFixedPointVolumeRayCastMapper_CastRays( void *arg )
{
  // Get the info out of the input structure
  int threadID    = ((vtkMultiThreader::ThreadInfo *)(arg))->ThreadID;
  int threadCount = ((vtkMultiThreader::ThreadInfo *)(arg))->NumberOfThreads;

  vtkCVFixedPointVolumeRayCastMapper *me = (vtkCVFixedPointVolumeRayCastMapper *)(((vtkMultiThreader::ThreadInfo *)arg)->UserData);
  
  if ( !me )
    {
    vtkGenericWarningMacro("Irrecoverable error: no mapper specified");
    return VTK_THREAD_RETURN_VALUE;
    }

  vtkVolume *vol = me->GetVolume();
  
  if ( me->GetBlendMode() == vtkVolumeMapper::MAXIMUM_INTENSITY_BLEND ||
       me->GetBlendMode() == vtkVolumeMapper::MINIMUM_INTENSITY_BLEND )
    {
    me->GetMIPHelper()->GenerateImage( threadID, threadCount, vol, me );
    }
  else
    {
    if ( me->GetShadingRequired() == 0 )
      {
      if ( me->GetGradientOpacityRequired() == 0 )
        {
        me->GetCompositeHelper()->GenerateImage( threadID, threadCount, vol, me );
        }
      else
        {
        me->GetCompositeGOHelper()->GenerateImage( threadID, threadCount, vol, me );
        }
      }
    else
      {
      if ( me->GetGradientOpacityRequired() == 0 )
        {
        me->GetCompositeShadeHelper()->GenerateImage( threadID, threadCount, vol, me );
        }
      else
        {
        me->GetCompositeGOShadeHelper()->GenerateImage( threadID, threadCount, vol, me );
        }
      }
    }
  
  return VTK_THREAD_RETURN_VALUE;
}


// Print method for vtkFixedPointVolumeRayCastMapper
void vtkCVFixedPointVolumeRayCastMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}
