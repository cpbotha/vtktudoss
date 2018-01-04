/*=========================================================================
   vtkCustomWidget
   =========================================================================*/

#include "vtk3DWidget.h"

class vtkTransform;

#ifndef __vtkCustomWidget_h
#define __vtkCustomWidget_h

class VTK_EXPORT vtkCustomWidget : public vtk3DWidget
{
public:
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkCustomWidget,vtk3DWidget);
  static vtkCustomWidget *New();
  virtual void PlaceWidget(double bds[6]);

protected:
  vtkCustomWidget(void);
  ~vtkCustomWidget(void);
};


#endif
