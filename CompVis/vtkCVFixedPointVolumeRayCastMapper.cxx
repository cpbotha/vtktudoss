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
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkTimerLog.h"
#include <math.h>

vtkCxxRevisionMacro(vtkCVFixedPointVolumeRayCastMapper, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkCVFixedPointVolumeRayCastMapper); 



// Construct a new vtkFixedPointVolumeRayCastMapper with default values
vtkCVFixedPointVolumeRayCastMapper::vtkCVFixedPointVolumeRayCastMapper()
{
  this->CVHelper   = vtkFixedPointVolumeRayCastCVHelper::New();
  this->CompVisMode = 1;
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

void vtkCVFixedPointVolumeRayCastMapper::Render( vtkRenderer *ren, vtkVolume *vol )
{
  this->Timer->StartTimer();
  
  // Since we are passing in a value of 0 for the multiRender flag
  // (this is a single render pass - not part of a multipass AMR render)
  // then we know the origin, spacing, and extent values will not
  // be used so just initialize everything to 0. No need to check
  // the return value of the PerImageInitialization method - since this
  // is not a multirender it will always return 1.
  double dummyOrigin[3]  = {0.0, 0.0, 0.0};
  double dummySpacing[3] = {0.0, 0.0, 0.0};
  int dummyExtent[6] = {0, 0, 0, 0, 0, 0};
  this->PerImageInitialization( ren, vol, 0,
                                dummyOrigin,
                                dummySpacing,
                                dummyExtent );

  this->PerVolumeInitialization( ren, vol );
  
  vtkRenderWindow *renWin=ren->GetRenderWindow();
  
  if ( renWin && renWin->CheckAbortStatus() )
    {
    this->AbortRender();
    return;
    }

  this->PerSubVolumeInitialization( ren, vol, 0 );
  if ( renWin && renWin->CheckAbortStatus() )
    {
    this->AbortRender();
    return;
    }
  
  this->RenderSubVolume();

  if ( renWin && renWin->CheckAbortStatus() )
    {
    this->AbortRender();
    return;
    }

  this->DisplayRenderedImage( ren, vol );
  
  this->Timer->StopTimer();
  this->TimeToDraw = this->Timer->GetElapsedTime();
  // If we've increased the sample distance, account for that in the stored time. Since we
  // don't get linear performance improvement, use a factor of .66
  this->StoreRenderTime( ren, vol, 
                         this->TimeToDraw * 
                         this->ImageSampleDistance * 
                         this->ImageSampleDistance *
                         ( 1.0 + 0.66*
                           (this->SampleDistance - this->OldSampleDistance) / 
                           this->OldSampleDistance ) );
  
  this->SampleDistance = this->OldSampleDistance;
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
          if (me->GetCompVisMode() > 0)
          {
            me->GetCVHelper()->GenerateImage(threadID, threadCount, vol, me);
          }
          else
          {
          me->GetCompositeShadeHelper()->GenerateImage( threadID, threadCount, vol, me );
          }
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
