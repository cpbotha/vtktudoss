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

  // Description:
  // if the confidence distance (3d component) is greater than this,
  // then we're in the context region, otherwise we're in the focus
  // region.
  vtkSetMacro(ConfidenceThreshold, double);
  vtkGetMacro(ConfidenceThreshold, double);

  // Description:
  // Specify what is shown in the focus area.
  vtkSetMacro(FocusTarget, int);
  vtkGetMacro(FocusTarget, int);
  void SetFocusTargetToC0() { this->FocusTarget = 0; };
  void SetFocusTargetToC1() { this->FocusTarget = 1; };
  void SetFocusTargetToMinC1() { this->FocusTarget = 2; };
  void SetFocusTargetToDiff() {this->FocusTarget = 3; };
  void SetFocusTargetToMagicLens() {this->FocusTarget = 4; };

  // Description:
  // Specify what is shown in the context area.
  vtkSetMacro(ContextTarget, int);
  vtkGetMacro(ContextTarget, int);
  void SetContextTargetToC0() { this->ContextTarget = 0; };
  void SetContextTargetToC1() { this->ContextTarget = 1; };
  void SetContextTargetToMinC1() { this->ContextTarget = 2; };
  void SetContextTargetToDiff() {this->ContextTarget = 3; };

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
  int FocusTarget;
  int ContextTarget;

private:
  vtkCVImageMapToColors(const vtkCVImageMapToColors&);  // Not implemented.
  void operator=(const vtkCVImageMapToColors&);  // Not implemented.
};

#endif







