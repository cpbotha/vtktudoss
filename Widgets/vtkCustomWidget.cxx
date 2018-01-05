#include "vtkCustomWidget.h"

#include "vtkActor.h"
#include "vtkCallbackCommand.h"
#include "vtkCellPicker.h"
#include "vtkDataSetMapper.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRendererCollection.h"
#include "vtkTransform.h"

// we're using snprintf to convert from integer to string representation
#ifdef _MSC_VER  //all MS compilers define this (version)
     #define snprintf _snprintf
#endif


vtkStandardNewMacro(vtkCustomWidget);

// Constructor
vtkCustomWidget::vtkCustomWidget() :
  TranslationHandle(vtkActor::New()),
  HandlePicker(vtkCellPicker::New()),
  State(vtkCustomWidget::Start)
{
  this->EventCallbackCommand->SetCallback(vtkCustomWidget::ProcessEvents);

  this->HandlePicker->SetTolerance(0.001);
  this->HandlePicker->AddPickList(this->TranslationHandle);
  this->HandlePicker->PickFromListOn();
}

// Destructor
vtkCustomWidget::~vtkCustomWidget()
{
  // Instantiated in constructor
  this->TranslationHandle->Delete();
}

void vtkCustomWidget::PlaceWidget(double bds[6])
{
  this->CreateTranslationHandle();
}

void vtkCustomWidget::CreateTranslationHandle()
{
  vtkCellArray *cells = vtkCellArray::New();
  cells->InsertNextCell(4);
  cells->InsertCellPoint(0);
  cells->InsertCellPoint(1);
  cells->InsertCellPoint(2);
  cells->InsertCellPoint(3);

  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(4);
  points->InsertPoint(0, 1, 1, 0);
  points->InsertPoint(1, -1, 1, 0);
  points->InsertPoint(2, -1, -1, 0);
  points->InsertPoint(3, 1, -1, 0);

  vtkPolyData* polyData = vtkPolyData::New();
  polyData->SetPolys(cells);
  polyData->SetPoints(points);
  cells->Delete();
  points->Delete();

  vtkDataSetMapper* mapper = vtkDataSetMapper::New();
  mapper->SetInputData(polyData);
  polyData->Delete();

  vtkProperty* orangeProperty = vtkProperty::New();
  orangeProperty->SetColor(1.0, 0.5, 0.0);
  orangeProperty->SetAmbient(1);
  orangeProperty->SetDiffuse(0);
  this->TranslationHandle->SetMapper(mapper);
  this->TranslationHandle->SetProperty(orangeProperty);
  mapper->Delete();
  orangeProperty->Delete();
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

    // Listen to the following events
    this->Interactor->AddObserver(vtkCommand::MouseMoveEvent,
                                  this->EventCallbackCommand, this->Priority);
    this->Interactor->AddObserver(vtkCommand::LeftButtonPressEvent,
                                  this->EventCallbackCommand, this->Priority);
    this->Interactor->AddObserver(vtkCommand::LeftButtonReleaseEvent,
                                  this->EventCallbackCommand, this->Priority);

    this->CurrentRenderer->AddActor(this->TranslationHandle);

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

    // Don't listen for events any more
    this->Interactor->RemoveObserver(this->EventCallbackCommand);

    this->CurrentRenderer->RemoveActor(this->TranslationHandle);

    this->InvokeEvent(vtkCommand::DisableEvent, NULL);
    this->SetCurrentRenderer(NULL);
  }

  this->Interactor->Render();
}

void vtkCustomWidget::ProcessEvents(vtkObject* object,
                                    unsigned long event,
                                    void* clientdata,
                                    void* calldata)
{
  vtkCustomWidget* self = reinterpret_cast<vtkCustomWidget*>(clientdata);

  //okay, let's do the right thing
  switch(event)
  {
  case vtkCommand::LeftButtonPressEvent:
    self->OnLeftButtonDown();
    break;
  case vtkCommand::LeftButtonReleaseEvent:
    self->OnLeftButtonUp();
    break;
  case vtkCommand::MouseMoveEvent:
    self->OnMouseMove();
    break;
  }
}

void vtkCustomWidget::OnLeftButtonDown()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->HandlePicker->Pick(x, y, 0.0, this->CurrentRenderer);
  vtkAssemblyPath* path = this->HandlePicker->GetPath();

  // Return if the handle wasn't clicked
  if ( path == NULL )
  {
    this->State = vtkCustomWidget::Start;
    return;
  }

  this->State = vtkCustomWidget::Moving;
  this->EventCallbackCommand->SetAbortFlag(1);
}

void vtkCustomWidget::OnLeftButtonUp()
{
  if (this->State != vtkCustomWidget::Moving) {
    return;
  }

  this->State = vtkCustomWidget::Start;
  this->EventCallbackCommand->SetAbortFlag(1);
  this->InvokeEvent(vtkCommand::EndInteractionEvent, NULL);
}

void vtkCustomWidget::OnMouseMove()
{
  if (this->State != vtkCustomWidget::Moving) {
    return;
  }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
}
