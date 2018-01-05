/*=========================================================================
   vtkCustomWidget
   =========================================================================*/

#include "vtk3DWidget.h"

class vtkActor;
class vtkCellPicker;
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

  // Override methods from parent
  virtual void SetEnabled(int);

protected:
  vtkCustomWidget(void);
  ~vtkCustomWidget(void);

  enum WidgetState
  {
    Start=0,
    Moving
  };

  // Handles the events
  static void ProcessEvents(vtkObject* object,
                            unsigned long event,
                            void* clientdata,
                            void* calldata);

  vtkActor* TranslationHandle;

  vtkCellPicker* HandlePicker;

  WidgetState State;

private:
  void CreateTranslationHandle();

  // ProcessEvents() dispatches to these methods.
  void OnLeftButtonDown();
  void OnLeftButtonUp();
  void OnMouseMove();
};


#endif
