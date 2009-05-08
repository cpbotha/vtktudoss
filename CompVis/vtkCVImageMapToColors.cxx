/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCVImageMapToColors.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCVImageMapToColors.h"

#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkScalarsToColors.h"
#include "vtkPointData.h"

vtkCxxRevisionMacro(vtkCVImageMapToColors, "$Revision: 1.0 $")
vtkStandardNewMacro(vtkCVImageMapToColors);
vtkCxxSetObjectMacro(vtkCVImageMapToColors,LookupTable2,vtkScalarsToColors);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkCVImageMapToColors::vtkCVImageMapToColors()
{
  this->LookupTable2 = NULL;
}

//----------------------------------------------------------------------------
vtkCVImageMapToColors::~vtkCVImageMapToColors()
{
  if (this->LookupTable2 != NULL)
    {
    this->LookupTable2->UnRegister(this);
    }


}

//----------------------------------------------------------------------------
unsigned long vtkCVImageMapToColors::GetMTime()
{
  unsigned long t1, t2;

  t1 = this->Superclass::GetMTime();
  if (this->LookupTable2)
    {
    t2 = this->LookupTable2->GetMTime();
    if (t2 > t1)
      {
      t1 = t2;
      }
    }
  return t1;
}


//----------------------------------------------------------------------------
// This non-templated function executes the filter for any type of data.

void vtkCVImageMapToColorsExecute(vtkCVImageMapToColors *self,
                                vtkImageData *inData, void *inPtr,
                                vtkImageData *outData,
                                unsigned char *outPtr,
                                int outExt[6], int id)
{
  int idxY, idxZ;
  int extX, extY, extZ;
  vtkIdType inIncX, inIncY, inIncZ;
  vtkIdType outIncX, outIncY, outIncZ;
  unsigned long count = 0;
  unsigned long target;
  int dataType = inData->GetScalarType();
  // size of scalar type in bytes
  int scalarSize = inData->GetScalarSize();
  int numberOfComponents,numberOfOutputComponents,outputFormat;
  int rowLength;
  vtkScalarsToColors *lookupTable = self->GetLookupTable();
  unsigned char *outPtr1;
  void *inPtr1;

  // find the region to loop over
  extX = outExt[1] - outExt[0] + 1;
  extY = outExt[3] - outExt[2] + 1;
  extZ = outExt[5] - outExt[4] + 1;

  target = static_cast<unsigned long>(extZ*extY/50.0);
  target++;

  // Get increments to march through data, taking into account number of components
  inData->GetContinuousIncrements(outExt, inIncX, inIncY, inIncZ);
  // because we are using void * and char * we must take care
  // of the scalar size in the increments
  inIncY *= scalarSize;
  inIncZ *= scalarSize;
  outData->GetContinuousIncrements(outExt, outIncX, outIncY, outIncZ);
  numberOfComponents = inData->GetNumberOfScalarComponents();
  numberOfOutputComponents = outData->GetNumberOfScalarComponents();
  outputFormat = self->GetOutputFormat();
  rowLength = extX*scalarSize*numberOfComponents;

  if (dataType != VTK_SHORT)
  {
    vtkErrorWithObjectMacro(self, <<"Input data needs to be type short.");
    return;
  }

  if (outputFormat != VTK_RGBA && outputFormat != VTK_RGB)
  {
    vtkErrorWithObjectMacro(self, <<"Output format needs to be RGB or RGBA.");
    return;
  }

  // Loop through output pixels
  outPtr1 = outPtr;
  inPtr1 = static_cast<void *>(
    static_cast<char *>(inPtr) + self->GetActiveComponent()*scalarSize);

  double c1[3], c2[3];

  for (idxZ = 0; idxZ < extZ; idxZ++)
    {
    for (idxY = 0; !self->AbortExecute && idxY < extY; idxY++)
      {
      if (!id)
        {
        if (!(count%target))
          {
          self->UpdateProgress(count/(50.0*target));
          }
        count++;
        }

      // inPtr1 is now a pointer to the first instance of the first active component for this row
      // outPtr1 is a pointer to the output RGB(A) data

      //lookupTable->MapScalarsThroughTable2(inPtr1,outPtr1,
      //                                     dataType,extX,numberOfComponents,
      //                                     outputFormat);

      unsigned char *outPtr3 = outPtr1;
      short *inPtr3 = static_cast<short*>(inPtr1);
      for (int i = 0; i < extX; i++)
      {
        // lookup first component at ...
        lookupTable->GetColor(static_cast<double>(*inPtr3), c1);
        outPtr3[0] = 0; // c1[0] * 255.0;
        outPtr3[1] = c1[1] * 255.0;
        outPtr3[2] = c1[2] * 255.0;
        if (outputFormat = VTK_RGBA) outPtr3[3] = 255;
        outPtr3 += numberOfOutputComponents;
        inPtr3 += numberOfComponents;
        
      }

      outPtr1 += outIncY + extX*numberOfOutputComponents;
      inPtr1 = static_cast<void *>(
        static_cast<char *>(inPtr1) + inIncY + rowLength);
      }
    outPtr1 += outIncZ;
    inPtr1 = static_cast<void *>(static_cast<char *>(inPtr1) + inIncZ);
    }
}

//----------------------------------------------------------------------------
// This method is passed a input and output data, and executes the filter
// algorithm to fill the output from the input.

void vtkCVImageMapToColors::ThreadedRequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *vtkNotUsed(outputVector),
  vtkImageData ***inData,
  vtkImageData **outData,
  int outExt[6], int id)
{
  
  
  
  void *inPtr = inData[0][0]->GetScalarPointerForExtent(outExt);
  void *outPtr = outData[0]->GetScalarPointerForExtent(outExt);

  vtkCVImageMapToColorsExecute(this, inData[0][0], inPtr,
                             outData[0], static_cast<unsigned char *>(outPtr),
                             outExt, id);
}

//----------------------------------------------------------------------------
void vtkCVImageMapToColors::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);


  os << indent << "LookupTable2: ";
  if (this->LookupTable2)
    {
    this->LookupTable2->PrintSelf(os << endl,indent.GetNextIndent());
    }
  else
    {
    os << "(none)\n";
    }
}





