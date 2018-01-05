#include "vtkCustomWidget.h"

#include "vtkSphereSource.h"
#include "vtkTransform.h"

// we're using snprintf to convert from integer to string representation
#ifdef _MSC_VER  //all MS compilers define this (version)
     #define snprintf _snprintf
#endif


vtkStandardNewMacro(vtkCustomWidget);

// Constructor
vtkCustomWidget::vtkCustomWidget()
{
}

// Destructor
vtkCustomWidget::~vtkCustomWidget()
{
}

void vtkCustomWidget::PlaceWidget(double bds[6])
{
}

void vtkCustomWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkCustomWidget::SetEnabled(int enabled)
{
  if (!this->Interactor)
  {
    vtkErrorMacro(<<"The interactor must be set prior to enabling/disabling widget");
    return;
  }

  if (enabled)
  {
    vtkDebugMacro(<<"Enabling widget")
    if (this->Enabled)
    {
      vtkDebugMacro(<<"Widget already enabled")
      return;
    }

    if (!this->CurrentRenderer)
    {
      this->SetCurrentRenderer((vtkRenderer *)(this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()));
      if (this->CurrentRenderer == NULL)
      {
        vtkDebugMacro(<<"Widget has no renderer")
        return;
      }
    }

    // Set enabled to true
    this->Enabled = 1;

    this->InvokeEvent(vtkCommand::EnableEvent, NULL);
  }
  else
  {
    vtkDebugMacro(<<"Disabling widget")
    if (!this->Enabled)
    {
      vtkDebugMacro(<<"Widget already disabled")
      return;
    }

    // Set enabled to false
    this->Enabled = 0;

    this->InvokeEvent(vtkCommand::DisableEvent, NULL);
    this->SetCurrentRenderer(NULL);
  }

  this->Interactor->Render();
}
