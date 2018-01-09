/*=========================================================================
   vtkProsthesisWidget
   =========================================================================*/

#ifndef __vtkProsthesisWidget_h
#define __vtkProsthesisWidget_h

#include "vtkAbstractWidget.h"

class vtkCustomRepresentation;

class VTK_EXPORT vtkProsthesisWidget : public vtkAbstractWidget
{
public:
  // Instantiate the class
  static vtkProsthesisWidget *New();

  // Standard methods for a VTK class
  vtkTypeMacro(vtkProsthesisWidget,vtkAbstractWidget);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  void SetRepresentation(vtkCustomRepresentation *r)
    {this->Superclass::SetWidgetRepresentation(reinterpret_cast<vtkWidgetRepresentation*>(r));}

  // Create the default widget representation if one is not set
  void CreateDefaultRepresentation() VTK_OVERRIDE;  

protected:
  vtkProsthesisWidget();
  ~vtkProsthesisWidget() VTK_OVERRIDE;

  // The states the widget can be in
  enum _WidgetState {Start=0, Active};
  int WidgetState;

  // Action callback functions
  static void SelectAction(vtkAbstractWidget*);
  static void EndSelectAction(vtkAbstractWidget*);
  static void MoveAction(vtkAbstractWidget*);

private:
  vtkProsthesisWidget(const vtkProsthesisWidget&) VTK_DELETE_FUNCTION;
  void operator=(const vtkProsthesisWidget&) VTK_DELETE_FUNCTION;
};

#endif
