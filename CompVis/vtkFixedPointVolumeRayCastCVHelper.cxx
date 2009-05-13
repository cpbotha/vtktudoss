/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFixedPointVolumeRayCastCVHelper.cxx,v $
  Language:  C++
  Date:      $Date: 2007-11-08 16:51:06 $
  Version:   $Revision: 1.9 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkFixedPointVolumeRayCastCVHelper.h"

#include "vtkImageData.h"
#include "vtkCommand.h"
#include "vtkCVFixedPointVolumeRayCastMapper.h"
//#include "vtkFixedPointVolumeRayCastHelper.h"
#include "vtkObjectFactory.h"
#include "vtkRenderWindow.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkFixedPointRayCastImage.h"
#include "vtkDataArray.h"

#include <math.h>

vtkCxxRevisionMacro(vtkFixedPointVolumeRayCastCVHelper, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkFixedPointVolumeRayCastCVHelper);

// Construct a new vtkFixedPointVolumeRayCastCVHelper with default values
vtkFixedPointVolumeRayCastCVHelper::vtkFixedPointVolumeRayCastCVHelper()
{
}

// Destruct a vtkFixedPointVolumeRayCastCVHelper - clean up any memory used
vtkFixedPointVolumeRayCastCVHelper::~vtkFixedPointVolumeRayCastCVHelper()
{
}


// This method is used when the interpolation type is nearest neighbor and
// the data has one component and scale == 1.0 and shift == 0.0. In the inner 
// loop we get the data value as an unsigned short, and use this index to 
// lookup a color and opacity for this sample. We then composite this into 
// the color computed so far along the ray, and check if we can terminate at
// this point (if the accumulated opacity is higher than some threshold).
// Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageOneSimpleNN( T *data, 
                                                   int threadID,
                                                   int threadCount,
                                                   vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                   vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();

  for ( k = 0; k < numSteps; k++ )
    {    
    if ( k )
      {
      VTKKWRCHelper_MoveToNextSampleShadeNN();    
      }
    
    VTKKWRCHelper_SpaceLeapCheck();    
    VTKKWRCHelper_CroppingCheckNN( pos );
    
    unsigned short val    = static_cast<unsigned short>(((*dptr)));
    VTKKWRCHelper_LookupColorUS( colorTable[0], scalarOpacityTable[0], val, tmp );
    if (tmp[3])
      {
      unsigned short normal = *dirPtr;
      VTKKWRCHelper_LookupShading( diffuseShadingTable[0], specularShadingTable[0], normal, tmp );    
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
      }
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has one component. In the inner loop we get the data value as
// an unsigned short using the scale/shift, and use this index to lookup
// a color and opacity for this sample. We then composite this into the
// color computed so far along the ray, and check if we can terminate at
// this point (if the accumulated opacity is higher than some threshold).
// Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageOneNN( T *data, 
                                             int threadID,
                                             int threadCount,
                                             vtkCVFixedPointVolumeRayCastMapper *mapper,
                                             vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();
  
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      VTKKWRCHelper_MoveToNextSampleShadeNN();    
      }
    
    VTKKWRCHelper_SpaceLeapCheck();    
    VTKKWRCHelper_CroppingCheckNN( pos );
    
    unsigned short val    = static_cast<unsigned short>(((*dptr) + shift[0])*scale[0]);
    VTKKWRCHelper_LookupColorUS( colorTable[0], scalarOpacityTable[0], val, tmp );
    if (tmp[3])
      {
      unsigned short normal = *dirPtr;
      VTKKWRCHelper_LookupShading( diffuseShadingTable[0], specularShadingTable[0], normal, tmp );
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
      }
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}


// This method is used when the interpolation type is nearest neighbor and
// the data has two components which are not considered independent. In the
// inner loop we compute the two unsigned short index values from the data
// values (using the scale/shift). We use the first index to lookup a color,
// and we use the second index to look up the opacity. We then composite 
// the color into the color computed so far along this ray, and check to 
// see if we can terminate here (if the opacity accumulated exceed some
// threshold). Finally we move to the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageTwoDependentNN( T *data, 
                                                      int threadID,
                                                      int threadCount,
                                                      vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                      vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();
  
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      VTKKWRCHelper_MoveToNextSampleShadeNN();    
      }
    
    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckNN( pos );
    
    unsigned short val[2];
    val[0] = static_cast<unsigned short>(((*(dptr  )) + shift[0])*scale[0]);
    val[1] = static_cast<unsigned short>(((*(dptr+1)) + shift[1])*scale[1]);
    
    tmp[3] = scalarOpacityTable[0][val[1]];
    if (tmp[3])
      {
      tmp[0] = static_cast<unsigned short>
        ((colorTable[0][3*val[0]  ]*tmp[3] + 0x7fff)>>(VTKKW_FP_SHIFT));
      tmp[1] = static_cast<unsigned short>
        ((colorTable[0][3*val[0]+1]*tmp[3] + 0x7fff)>>(VTKKW_FP_SHIFT));
      tmp[2] = static_cast<unsigned short>
        ((colorTable[0][3*val[0]+2]*tmp[3] + 0x7fff)>>(VTKKW_FP_SHIFT));
    
      unsigned short normal = *dirPtr;    
      VTKKWRCHelper_LookupShading( diffuseShadingTable[0], specularShadingTable[0], normal, tmp );
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
      }
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has four components which are not considered independent . This
// means that the first three components directly represent color, and this
// data must be of unsigned char type. In the inner loop we directly access
// the four data values (no scale/shift is needed). The first three are the
// color of this sample and the fourth is used to look up an opacity in the
// scalar opacity transfer function. We then composite this color into the
// color we have accumulated so far along the ray, and check if we can 
// terminate here (if our accumulated opacity has exceed some threshold). 
// Finally we move onto the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageFourDependentNN( T *data, 
                                                       int threadID,
                                                       int threadCount,
                                                       vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                       vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeOneNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
  VTKKWRCHelper_SpaceLeapSetup();

  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      VTKKWRCHelper_MoveToNextSampleShadeNN();    
      }
    
    VTKKWRCHelper_SpaceLeapCheck();    
    VTKKWRCHelper_CroppingCheckNN( pos );
    
    unsigned short val[4];
    val[0] = *(dptr  );
    val[1] = *(dptr+1);
    val[2] = *(dptr+2);
    val[3] = static_cast<unsigned short>(((*(dptr+3)) + shift[3])*scale[3]);
    
    tmp[3] = scalarOpacityTable[0][val[3]];
    if ( tmp[3] )
      {
      tmp[0] = (val[0]*tmp[3]+0x7f)>>(8);
      tmp[1] = (val[1]*tmp[3]+0x7f)>>(8);
      tmp[2] = (val[2]*tmp[3]+0x7f)>>(8);
      
      unsigned short normal   = *dirPtr;    
      VTKKWRCHelper_LookupShading( diffuseShadingTable[0], specularShadingTable[0], normal, tmp );
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
      }
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is nearest neighbor and
// the data has more than one component and the components are considered to
// be independent. In the inner loop we access each component value, using
// the scale/shift to turn the data value into an unsigned short index. We
// then lookup the color/opacity for each component and combine them according
// to the weighting value for each component. We composite this resulting
// color into the color already accumulated for this ray, and we check
// wether we can terminate here (if the accumulated opacity exceeds some
// threshold). Finally we increment to the next sample on the ray.
//
// TODO: short circuit calculations when opacity is 0
template <class T>
void vtkFixedPointCVHelperGenerateImageIndependentNN( T *data, 
                                                     int threadID,
                                                     int threadCount,
                                                     vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                     vtkVolume *vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartShadeNN();
  VTKKWRCHelper_InitializeCompositeMultiNN();
  VTKKWRCHelper_InitializeCompositeShadeNN();
      
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      VTKKWRCHelper_MoveToNextSampleShadeNN();    
      }
    
    VTKKWRCHelper_CroppingCheckNN( pos );
    
    unsigned short normal[4];
    for ( c = 0; c < components; c++ )
      {
      val[c]    = static_cast<unsigned short>(((*(dptr+c)) + shift[c])*scale[c]);
      normal[c] = *(dirPtr+c);
      }
    
    VTKKWRCHelper_LookupAndCombineIndependentColorsShadeUS( colorTable, scalarOpacityTable,
                                                            diffuseShadingTable,
                                                            specularShadingTable,
                                                            val, normal, weights, 
                                                            components, tmp );        

    if ( tmp[3] )
      {
      VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
      }
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear and the data 
// has one component and scale = 1.0 and shift = 0.0. In the inner loop we 
// get the data value for the eight cell corners (if we have changed cells)
// as an unsigned short (the range must be right and we don't need the 
// scale/shift). We compute our weights within the cell according to our 
// fractional position within the cell, apply trilinear interpolation to
// compute the index, and use this index to lookup a color and opacity for 
// this sample. We then composite this into the color computed so far along 
// the ray, and check if we can terminate at this point (if the accumulated 
// opacity is higher than some threshold). Finally we move on to the next 
// sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageOneSimpleTrilin( T *data, 
                                                       int threadID,
                                                       int threadCount,
                                                       vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                       vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();

  int needToSampleDirection = 0;
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      mapper->FixedPointIncrement( pos, dir );
      }
    
    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin( pos );
    
    mapper->ShiftVectorDown( pos, spos );
    if ( spos[0] != oldSPos[0] ||
         spos[1] != oldSPos[1] ||
         spos[2] != oldSPos[2] )
      {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];
      
      dptr = data + spos[0]*inc[0] + spos[1]*inc[1] + spos[2]*inc[2];
      VTKKWRCHelper_GetCellScalarValuesSimple( dptr );
      dirPtrABCD = gradientDir[spos[2]  ] + spos[0]*dInc[0] + spos[1]*dInc[1];
      dirPtrEFGH = gradientDir[spos[2]+1] + spos[0]*dInc[0] + spos[1]*dInc[1];
      needToSampleDirection = 1;
      }
          
    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);
    
    VTKKWRCHelper_LookupColorUS( colorTable[0], scalarOpacityTable[0], val, tmp );
    if ( needToSampleDirection )
      {
      VTKKWRCHelper_GetCellDirectionValues( dirPtrABCD, dirPtrEFGH );
      needToSampleDirection = 0;
      }
    
    VTKKWRCHelper_InterpolateShading( diffuseShadingTable[0], specularShadingTable[0], tmp ); 
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );    
    }
      
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear and the data 
// has one component and scale != 1.0 or shift != 0.0. In the inner loop we 
// get the data value for the eight cell corners (if we have changed cells)
// as an unsigned short (we use the scale/shift to ensure the correct range). 
// We compute our weights within the cell according to our fractional position 
// within the cell, apply trilinear interpolation to compute the index, and use 
// this index to lookup a color and opacity for this sample. We then composite 
// this into the color computed so far along the ray, and check if we can 
// terminate at this point (if the accumulated opacity is higher than some 
// threshold). Finally we move on to the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageOneTrilin( T *data, 
                                                 int threadID,
                                                 int threadCount,
                                                 vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                 vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeOneTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();
        
  int needToSampleDirection = 0;
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      mapper->FixedPointIncrement( pos, dir );
      }

    VTKKWRCHelper_SpaceLeapCheck();    
    VTKKWRCHelper_CroppingCheckTrilin( pos );
    
    mapper->ShiftVectorDown( pos, spos );    
    if ( spos[0] != oldSPos[0] ||
         spos[1] != oldSPos[1] ||
         spos[2] != oldSPos[2] )
      {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];
      
      
      dptr = data + spos[0]*inc[0] + spos[1]*inc[1] + spos[2]*inc[2];
      VTKKWRCHelper_GetCellScalarValues( dptr, scale[0], shift[0] );
      dirPtrABCD = gradientDir[spos[2]  ] + spos[0]*dInc[0] + spos[1]*dInc[1];
      dirPtrEFGH = gradientDir[spos[2]+1] + spos[0]*dInc[0] + spos[1]*dInc[1];
      needToSampleDirection = 1;
      }
    
    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalar(val);
    
    VTKKWRCHelper_LookupColorUS( colorTable[0], scalarOpacityTable[0], val, tmp );
    if ( needToSampleDirection )
      {
      VTKKWRCHelper_GetCellDirectionValues( dirPtrABCD, dirPtrEFGH );
      needToSampleDirection = 0;
      }
    VTKKWRCHelper_InterpolateShading( diffuseShadingTable[0], specularShadingTable[0], tmp ); 
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );    
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}


// This method is used when the interpolation type is linear, the data has 
// two components and the components are not considered independent. In the
// inner loop we get the data value for the eight cell corners (if we have 
// changed cells) for both components as an unsigned shorts (we use the 
// scale/shift to ensure the correct range). We compute our weights within 
// the cell according to our fractional position within the cell, and apply 
// trilinear interpolation to compute the two index value. We use the first 
// index to lookup a color and the second to look up an opacity for this sample. 
// We then composite this into the color computed so far along the ray, and 
// check if we can terminate at this point (if the accumulated opacity is 
// higher than some threshold). Finally we move on to the next sample along 
// the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageTwoDependentTrilin( T *data, 
                                                          int threadID,
                                                          int threadCount,
                                                          vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                          vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();
    
  int needToSampleDirection = 0;
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      mapper->FixedPointIncrement( pos, dir );
      }
    
    VTKKWRCHelper_SpaceLeapCheck();    
    VTKKWRCHelper_CroppingCheckTrilin( pos );
    
    mapper->ShiftVectorDown( pos, spos );
    if ( spos[0] != oldSPos[0] ||
         spos[1] != oldSPos[1] ||
         spos[2] != oldSPos[2] )
      {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];
      
      dptr = data + spos[0]*inc[0] + spos[1]*inc[1] + spos[2]*inc[2];
      VTKKWRCHelper_GetCellComponentScalarValues( dptr, 0, scale[0], shift[0] );
      
      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues( dptr, 1, scale[1], shift[1] );
      
      dirPtrABCD = gradientDir[spos[2]  ] + spos[0]*dInc[0] + spos[1]*dInc[1];
      dirPtrEFGH = gradientDir[spos[2]+1] + spos[0]*dInc[0] + spos[1]*dInc[1];
      needToSampleDirection = 1;
      }
    
    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent( val, c, 2 );
    
    tmp[3] = scalarOpacityTable[0][val[1]];
    if ( !tmp[3] )
      {
      continue;
      }
    
    if ( needToSampleDirection )
      {
      VTKKWRCHelper_GetCellDirectionValues( dirPtrABCD, dirPtrEFGH );
      needToSampleDirection = 0;
      }
    
    tmp[0] = static_cast<unsigned short>
      ((colorTable[0][3*val[0]  ]*tmp[3] + 0x7fff)>>(VTKKW_FP_SHIFT));
    tmp[1] = static_cast<unsigned short>
      ((colorTable[0][3*val[0]+1]*tmp[3] + 0x7fff)>>(VTKKW_FP_SHIFT));
    tmp[2] = static_cast<unsigned short>
      ((colorTable[0][3*val[0]+2]*tmp[3] + 0x7fff)>>(VTKKW_FP_SHIFT));
    
    VTKKWRCHelper_InterpolateShading( diffuseShadingTable[0], specularShadingTable[0], tmp ); 
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
    
    }

  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}



// This method is used when the interpolation type is linear, the data has 
// four components and the components are not considered independent. In the
// inner loop we get the data value for the eight cell corners (if we have 
// changed cells) for all components as an unsigned shorts (we don't have to
// use the scale/shift because only unsigned char data is supported for four
// component data when the components are not independent). We compute our 
// weights within the cell according to our fractional position within the cell, 
// and apply trilinear interpolation to compute a value for each component. We 
// use the first three directly as the color of the sample, and the fourth is
// used to look up an opacity for this sample. We then composite this into the 
// color computed so far along the ray, and check if we can terminate at this 
// point (if the accumulated opacity is higher than some threshold). Finally we 
// move on to the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageFourDependentTrilin( T *data, 
                                                           int threadID,
                                                           int threadCount,
                                                           vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                           vtkVolume *vol)
{
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeOneShadeTrilin();
  VTKKWRCHelper_SpaceLeapSetup();
    
  int needToSampleDirection = 0;
  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      mapper->FixedPointIncrement( pos, dir );
      }

    VTKKWRCHelper_SpaceLeapCheck();
    VTKKWRCHelper_CroppingCheckTrilin( pos );
    
    mapper->ShiftVectorDown( pos, spos );
    if ( spos[0] != oldSPos[0] ||
         spos[1] != oldSPos[1] ||
         spos[2] != oldSPos[2] )
      {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];
      
      dptr = data + spos[0]*inc[0] + spos[1]*inc[1] + spos[2]*inc[2];
      VTKKWRCHelper_GetCellComponentRawScalarValues( dptr, 0 );
      
      dptr++;
      VTKKWRCHelper_GetCellComponentRawScalarValues( dptr, 1 );
      
      dptr++;
      VTKKWRCHelper_GetCellComponentRawScalarValues( dptr, 2 );
      
      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues( dptr, 3, scale[3], shift[3] );
      
      dirPtrABCD = gradientDir[spos[2]  ] + spos[0]*dInc[0] + spos[1]*dInc[1];
      dirPtrEFGH = gradientDir[spos[2]+1] + spos[0]*dInc[0] + spos[1]*dInc[1];
      needToSampleDirection = 1;
      }
    
    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent( val, c, components );
    
    tmp[3] = scalarOpacityTable[0][val[3]];
    if ( !tmp[3] )
      {
      continue;
      }
    
    if ( needToSampleDirection )
      {
      VTKKWRCHelper_GetCellDirectionValues( dirPtrABCD, dirPtrEFGH );
      needToSampleDirection = 0;
      }
    
    tmp[0] = (val[0]*tmp[3]+0x7f)>>8;
    tmp[1] = (val[1]*tmp[3]+0x7f)>>8;
    tmp[2] = (val[2]*tmp[3]+0x7f)>>8;
    
    VTKKWRCHelper_InterpolateShading( diffuseShadingTable[0], specularShadingTable[0], tmp ); 
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}

// This method is used when the interpolation type is linear, the data has 
// more than one component and the components are considered independent. In 
// the inner loop we get the data value for the eight cell corners (if we have 
// changed cells) for all components as an unsigned shorts (we have to use the
// scale/shift to ensure that we obtained unsigned short indices) We compute our 
// weights within the cell according to our fractional position within the cell, 
// and apply trilinear interpolation to compute a value for each component. We 
// look up a color/opacity for each component and blend them according to the 
// component weights. We then composite this resulting color into the 
// color computed so far along the ray, and check if we can terminate at this 
// point (if the accumulated opacity is higher than some threshold). Finally we 
// move on to the next sample along the ray.
template <class T>
void vtkFixedPointCVHelperGenerateImageIndependentTrilin( T *data, 
                                                         int threadID,
                                                         int threadCount,
                                                         vtkCVFixedPointVolumeRayCastMapper *mapper,
                                                         vtkVolume *vol)
{
  VTKKWRCHelper_InitializeWeights();
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();
  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeMultiShadeTrilin();

  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      mapper->FixedPointIncrement( pos, dir );
      }
    
    VTKKWRCHelper_CroppingCheckTrilin( pos );
    
    mapper->ShiftVectorDown( pos, spos );
    if ( spos[0] != oldSPos[0] ||
         spos[1] != oldSPos[1] ||
         spos[2] != oldSPos[2] )
      {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];
      
      dptr = data + spos[0]*inc[0] + spos[1]*inc[1] + spos[2]*inc[2];
      VTKKWRCHelper_GetCellComponentScalarValues( dptr, 0, scale[0], shift[0] );
      
      dptr++;
      VTKKWRCHelper_GetCellComponentScalarValues( dptr, 1, scale[1], shift[1] );
      
      if ( components > 2 )
        {
        dptr++;
        VTKKWRCHelper_GetCellComponentScalarValues( dptr, 2, scale[2], shift[2] );
        if ( components > 3 )
          {
          dptr++;
          VTKKWRCHelper_GetCellComponentScalarValues( dptr, 3, scale[3], shift[3] );
          }
        }
      
      
      dirPtrABCD = gradientDir[spos[2]  ] + spos[0]*dInc[0] + spos[1]*dInc[1];
      dirPtrEFGH = gradientDir[spos[2]+1] + spos[0]*dInc[0] + spos[1]*dInc[1];
      VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 0 ); 

      dirPtrABCD++;
      dirPtrEFGH++;
      VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 1 ); 
      
      if ( components > 2 )
        {
        dirPtrABCD++;
        dirPtrEFGH++;
        VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 2 ); 
        if ( components > 3 )
          {
          dirPtrABCD++;
          dirPtrEFGH++;
          VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 3 ); 
          }
        }
      }
    
    VTKKWRCHelper_ComputeWeights(pos);
    VTKKWRCHelper_InterpolateScalarComponent( val, c, components );
    
    VTKKWRCHelper_LookupAndCombineIndependentColorsInterpolateShadeUS( 
      colorTable, scalarOpacityTable, diffuseShadingTable,
      specularShadingTable, val, weights, components, tmp );  
    
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
    }
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}


// volume render 3-component short data
template <class T>
void vtkFixedPointCVHelperGenerateImageCV1( T *data, 
                                            int threadID,
                                            int threadCount,
                                            vtkCVFixedPointVolumeRayCastMapper *mapper,
                                            vtkVolume *vol)
{
  VTKKWRCHelper_InitializeWeights();

  // this defines c
  VTKKWRCHelper_InitializationAndLoopStartShadeTrilin();

  VTKKWRCHelper_InitializeCompositeMultiTrilin();
  VTKKWRCHelper_InitializeCompositeMultiShadeTrilin();

  for ( k = 0; k < numSteps; k++ )
    {
    if ( k )
      {
      mapper->FixedPointIncrement( pos, dir );
      }
    
    VTKKWRCHelper_CroppingCheckTrilin( pos );
    
    mapper->ShiftVectorDown( pos, spos );
    if ( spos[0] != oldSPos[0] ||
         spos[1] != oldSPos[1] ||
         spos[2] != oldSPos[2] )
      {
      oldSPos[0] = spos[0];
      oldSPos[1] = spos[1];
      oldSPos[2] = spos[2];

      // address of first component
      dptr = data + spos[0]*inc[0] + spos[1]*inc[1] + spos[2]*inc[2];

      // third component (distance) - we just use nearest neighbour in this case.
      val[2] = static_cast<short>(((*(dptr+2)) + shift[2])*scale[2]);

      dirPtrABCD = gradientDir[spos[2]  ] + spos[0]*dInc[0] + spos[1]*dInc[1];
      dirPtrEFGH = gradientDir[spos[2]+1] + spos[0]*dInc[0] + spos[1]*dInc[1];

      if (val[2] < 32)
        {
        // cell corner scalar values for first component
        VTKKWRCHelper_GetCellComponentScalarValues( dptr, 0, scale[0], shift[0] );
      
        // cell corner gradient values for first component

        VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 0 ); 

        }
      else
        {
        // second component
        dptr++;
        VTKKWRCHelper_GetCellComponentScalarValues( dptr, 1, scale[1], shift[1] );

        // second component
        dirPtrABCD++;
        dirPtrEFGH++;
        VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 1 ); 
        }
      
      // third component (distance) - won't need gradients, just using as mask
      //dirPtrABCD++;
      //dirPtrEFGH++;
      //VTKKWRCHelper_GetCellComponentDirectionValues( dirPtrABCD, dirPtrEFGH, 2 ); 


      } // if (spos[0] != ...

    
    
    // calculate weights for combining corner values
    VTKKWRCHelper_ComputeWeights(pos);

    // use weights and corner values to interpolate scalar values
    // for all N components (3d param; we only need 2).  Use c (2nd param) as a loop control var,
    // stuff interpolated values in 'val' (1st param) array.
    //VTKKWRCHelper_InterpolateScalarComponent( val, c, 2 );

    // first component inside focus, second component outside
    if (val[2] < 32)
      VTKKWRCHelper_InterpolateSingleScalarComponent( val, 0 )
    else
      VTKKWRCHelper_InterpolateSingleScalarComponent( val, 1);
    
    //VTKKWRCHelper_LookupAndCombineIndependentColorsInterpolateShadeUS( 
    //  colorTable, scalarOpacityTable, diffuseShadingTable,
    //  specularShadingTable, val, weights, components, tmp );  
    // cpbotha: call above is usually used, we're going to replace it
    // with inline code and MO-DI-FY!

#define COLORTABLE colorTable
#define SOTABLE scalarOpacityTable
#define DTABLE diffuseShadingTable
#define STABLE specularShadingTable
#define SCALAR val
#define WEIGHTS weights
//#define COMPONENTS components
#define COMPONENTS 2
#define COLOR tmp


  unsigned int _tmp[4] = {0,0,0,0};                                                                                             \
  unsigned short _alpha[4] = {0,0,0,0};  
  int _idx;
  
  if (val[2] < 32)
    _idx = 0;
  else
    _idx = 1;


    _alpha[_idx] = static_cast<unsigned short>(SOTABLE[_idx][SCALAR[_idx]]);                                      \
    
                                                                                                                                \
  if ( !_alpha[_idx] ) {continue;}    
  
    
      _tmp[0] = static_cast<unsigned short>(((COLORTABLE[_idx][3*SCALAR[_idx]  ])*_alpha[_idx] + 0x7fff)>>(VTKKW_FP_SHIFT));   \
      _tmp[1] = static_cast<unsigned short>(((COLORTABLE[_idx][3*SCALAR[_idx]+1])*_alpha[_idx] + 0x7fff)>>(VTKKW_FP_SHIFT));   \
      _tmp[2] = static_cast<unsigned short>(((COLORTABLE[_idx][3*SCALAR[_idx]+2])*_alpha[_idx] + 0x7fff)>>(VTKKW_FP_SHIFT));   \
      _tmp[3] = _alpha[_idx];                                                                                                  \
      VTKKWRCHelper_InterpolateShadingComponent( DTABLE, STABLE, _tmp, _idx );                                                 \

  if (!_tmp[3]) {continue;}                                                                                                     \
  COLOR[0] = (_tmp[0]>32767)?(32767):(_tmp[0]);                                                                                 \
  COLOR[1] = (_tmp[1]>32767)?(32767):(_tmp[1]);                                                                                 \
  COLOR[2] = (_tmp[2]>32767)?(32767):(_tmp[2]);                                                                                 \
  COLOR[3] = (_tmp[3]>32767)?(32767):(_tmp[3]);

    // end of MO-DI-FY!

    
    VTKKWRCHelper_CompositeColorAndCheckEarlyTermination( color, tmp, remainingOpacity );
    } // for (k = ...
  
  VTKKWRCHelper_SetPixelColor( imagePtr, color, remainingOpacity );  
  VTKKWRCHelper_IncrementAndLoopEnd();
}



void vtkFixedPointVolumeRayCastCVHelper::GenerateImage(
  int threadID,
  int threadCount,
  vtkVolume *vol,
  vtkCVFixedPointVolumeRayCastMapper *mapper )
{
  void *data     = mapper->GetCurrentScalars()->GetVoidPointer(0);
  int scalarType = mapper->GetCurrentScalars()->GetDataType();

  // Nearest Neighbor interpolate
  if ( mapper->ShouldUseNearestNeighborInterpolation( vol ) )
    {
    // One component data
    if ( mapper->GetCurrentScalars()->GetNumberOfComponents() == 1 )
      {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if ( mapper->GetTableScale()[0] == 1.0 &&
           mapper->GetTableShift()[0] == 0.0 )
        {
        switch ( scalarType )
          {
          vtkTemplateMacro( 
            vtkFixedPointCVHelperGenerateImageOneSimpleNN(
              static_cast<VTK_TT *>(data),
              threadID, threadCount, mapper, vol) );
          }
        }
      else
        {
        switch ( scalarType )
          {
          vtkTemplateMacro( 
            vtkFixedPointCVHelperGenerateImageOneNN(
              static_cast<VTK_TT *>(data),
              threadID, threadCount, mapper, vol) );
          }
        }
      }
    // More that one independent components
    else if ( vol->GetProperty()->GetIndependentComponents() )
      {
      switch ( scalarType )
        {
        vtkTemplateMacro( 
          vtkFixedPointCVHelperGenerateImageIndependentNN(
            static_cast<VTK_TT *>(data),
            threadID, threadCount, mapper, vol) );
        }
      }
    // Dependent (color) components
    else
      {
      // Two components - the first specifies color (through a lookup table)
      // and the second specified opacity (through a lookup table)
      if ( mapper->GetCurrentScalars()->GetNumberOfComponents() == 2 )
        {
        switch ( scalarType )
          {
          vtkTemplateMacro( 
            vtkFixedPointCVHelperGenerateImageTwoDependentNN(
              static_cast<VTK_TT *>(data),
              threadID, threadCount, mapper, vol) );
          }
        }
      // Four components - they must be unsigned char, the first three directly
      // specify color and the fourth specifies opacity (through a lookup
      // table)
      else
        {
        if ( scalarType == VTK_UNSIGNED_CHAR )
          {
          vtkFixedPointCVHelperGenerateImageFourDependentNN( 
            static_cast<unsigned char *>(data), threadID, threadCount, mapper,
            vol );
          }
        else
          {
          vtkErrorMacro("Four component dependent data must be unsigned char");
          }
        }
      }
    }
  // Trilinear Interpolation
  else
    {
    // One component
    if ( mapper->GetCurrentScalars()->GetNumberOfComponents() == 1 )
      {
      // Scale == 1.0 and shift == 0.0 - simple case (faster)
      if ( mapper->GetTableScale()[0] == 1.0 &&
           mapper->GetTableShift()[0] == 0.0 )
        {
        switch ( scalarType )
          {
          vtkTemplateMacro( 
            vtkFixedPointCVHelperGenerateImageOneSimpleTrilin(
              static_cast<VTK_TT *>(data),
              threadID, threadCount, mapper, vol) );
          }
        }
      // Scale != 1.0 or shift != 0.0 - must apply scale/shift in inner loop
      else
        {
        switch ( scalarType )
          {
          vtkTemplateMacro( 
            vtkFixedPointCVHelperGenerateImageOneTrilin(
              static_cast<VTK_TT *>(data),
              threadID, threadCount, mapper, vol) );
          }
        }
      }
    // Indepedent components (more than one)
    else if ( vol->GetProperty()->GetIndependentComponents() )
      {
        // cpbotha: we're going to integrate here...
        // logic in the base vtkFPVRCMapper dictates that it prepares as
        // many TFs / LUTs as there are components.  We are expecting
        // three-component datasets, so we'll just have to add a dummy¡
        // third transfer function.  One could also consider overriding
        // UpdateColorTable() and its caller PerVolumeInitialization()
        //
        // yay! seems we DON'T have to add a dummy.
        if (mapper->GetCompVisMode() > 0)
          {
          // for now, we only do signed shorts
          if ( scalarType == VTK_SHORT && mapper->GetCurrentScalars()->GetNumberOfComponents() >= 3)
            {
              vtkFixedPointCVHelperGenerateImageCV1(
                static_cast<short*>(data),
                threadID, threadCount, mapper, vol);
            }
          else
            {
              // work-around for Python threading bug
              if (threadID == 0)
                vtkErrorMacro("CompVis mode data must have at least 3 short components.");
            }
          }
        else
          {
          switch ( scalarType )
            {
            vtkTemplateMacro( 
              vtkFixedPointCVHelperGenerateImageIndependentTrilin(
                static_cast<VTK_TT *>(data),
                threadID, threadCount, mapper, vol) );
            }
          }
      }
    // Dependent components
    else
      {
      // Two components - the first specifies color (through a lookup table)
      // and the second specified opacity (through a lookup table)
      if ( mapper->GetCurrentScalars()->GetNumberOfComponents() == 2 )
        {
        switch ( scalarType )
          {
          vtkTemplateMacro( 
            vtkFixedPointCVHelperGenerateImageTwoDependentTrilin(
              static_cast<VTK_TT *>(data),
              threadID, threadCount, mapper, vol) );
          }
        }
      // Four components - they must be unsigned char, the first three directly
      // specify color and the fourth specifies opacity (through a lookup
      // table)
      else
        {
        if ( scalarType == VTK_UNSIGNED_CHAR )
          {
          vtkFixedPointCVHelperGenerateImageFourDependentTrilin( 
            static_cast<unsigned char *>(data), threadID, threadCount, mapper,
            vol );
          }
        else
          {
          vtkErrorMacro("Four component dependent data must be unsigned char");
          }
        }
      }
    }
}

// Print method for vtkFixedPointVolumeRayCastCVHelper
void vtkFixedPointVolumeRayCastCVHelper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

