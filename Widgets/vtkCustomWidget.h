/*=========================================================================
   vtkCustomWidget
   =========================================================================*/

#ifndef __vtkCustomWidget_h
#define __vtkCustomWidget_h

#include "vtkAbstractWidget.h"

class vtkActor;
class vtkCellPicker;
class vtkTransform;

class VTK_EXPORT vtkCustomWidget : public vtkAbstractWidget
{
public:
  void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkCustomWidget, vtkAbstractWidget);
  static vtkCustomWidget *New();
  virtual void CreateDefaultRepresentation() VTK_OVERRIDE;

  // Override methods from parent
  virtual void SetEnabled(int);

  // Public member variables
  double Center[3];

  // Getters and setters
  vtkSetVectorMacro(Center, double, 3);
  vtkGetVectorMacro(Center, double, 3);

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

  void SubtractVectors(double vector1[3],
                       double vector2[3],
                       double* resultvector);
  void AddVectors(double vector1[3],
                  double vector2[3],
                  double* resultvector);
  void CopyVector(double inputVector[3],
                  double outputVector[3]);

  // Used to store the center when the user starts to drag the
  // translation handle.
  double StartCenter[3];
  double StartPick[3];
};


#endif
