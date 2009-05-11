#ifndef __vtkCVImageMapToColors_h
#define __vtkCVImageMapToColors_h


#include "vtkImageMapToColors.h"

class vtkScalarsToColors;

class VTK_EXPORT vtkCVImageMapToColors : public vtkImageMapToColors
{
public:
  static vtkCVImageMapToColors *New();
  vtkTypeRevisionMacro(vtkCVImageMapToColors,vtkImageMapToColors);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the lookup table.
  virtual void SetLookupTable2(vtkScalarsToColors*);
  vtkGetObjectMacro(LookupTable2,vtkScalarsToColors);

  // Description:
  // We need to check the modified time of the lookup table too.
  virtual unsigned long GetMTime();

  vtkSetMacro(ConfidenceThreshold, double);
  vtkGetMacro(ConfidenceThreshold, double);


protected:
  vtkCVImageMapToColors();
  ~vtkCVImageMapToColors();

  virtual int RequestData(vtkInformation *request,
                          vtkInformationVector **inputVector,
                          vtkInformationVector *outputVector);

  void ThreadedRequestData(vtkInformation *request,
                           vtkInformationVector **inputVector,
                           vtkInformationVector *outputVector,
                           vtkImageData ***inData, vtkImageData **outData,
                           int extent[6], int id);

  vtkScalarsToColors *LookupTable2;

  double ConfidenceThreshold;

private:
  vtkCVImageMapToColors(const vtkCVImageMapToColors&);  // Not implemented.
  void operator=(const vtkCVImageMapToColors&);  // Not implemented.
};

#endif







