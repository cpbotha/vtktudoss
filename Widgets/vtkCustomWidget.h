/*=========================================================================
   vtkCustomWidget
   =========================================================================*/

#ifndef __vtkCustomWidget_h
#define __vtkCustomWidget_h

#include "vtkAbstractWidget.h"

class vtkCustomRepresentation;

class VTK_EXPORT vtkCustomWidget : public vtkAbstractWidget
{
public:
  // Instantiate the class
  static vtkCustomWidget *New();

  // Standard methods for a VTK class
  vtkTypeMacro(vtkCustomWidget,vtkAbstractWidget);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  void SetRepresentation(vtkCustomRepresentation *r)
    {this->Superclass::SetWidgetRepresentation(reinterpret_cast<vtkWidgetRepresentation*>(r));}

  // Create the default widget representation if one is not set
  void CreateDefaultRepresentation() VTK_OVERRIDE;  

protected:
  vtkCustomWidget();
  ~vtkCustomWidget() VTK_OVERRIDE;

  // The states the widget can be in
  enum _WidgetState {Start=0, Active};
  int WidgetState;

  // Action callback functions
  static void SelectAction(vtkAbstractWidget*);
  static void EndSelectAction(vtkAbstractWidget*);
  static void MoveAction(vtkAbstractWidget*);

private:
  vtkCustomWidget(const vtkCustomWidget&) VTK_DELETE_FUNCTION;
  void operator=(const vtkCustomWidget&) VTK_DELETE_FUNCTION;
};

#endif
