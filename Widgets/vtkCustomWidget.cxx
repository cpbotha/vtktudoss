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
  this->Center[0] = this->Center[1] = this->Center[2] = 0.0;
  this->CopyVector(this->Center, this->StartCenter);
  this->CopyVector(this->Center, this->StartPick);
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

void vtkCustomWidget::CreateDefaultRepresentation()
{
  this->CreateTranslationHandle();
}

vtkPolyData* vtkCustomWidget::CreateTranslationHandlePolydata()
{
  vtkCellArray* cells = vtkCellArray::New();
  cells->InsertNextCell(4);
  cells->InsertCellPoint(0);
  cells->InsertCellPoint(1);
  cells->InsertCellPoint(2);
  cells->InsertCellPoint(3);
  cells->InsertNextCell(4);
  cells->InsertCellPoint(0);
  cells->InsertCellPoint(3);
  cells->InsertCellPoint(4);
  cells->InsertCellPoint(5);
  cells->InsertNextCell(4);
  cells->InsertCellPoint(0);
  cells->InsertCellPoint(5);
  cells->InsertCellPoint(6);
  cells->InsertCellPoint(7);
  cells->InsertNextCell(4);
  cells->InsertCellPoint(0);
  cells->InsertCellPoint(7);
  cells->InsertCellPoint(8);
  cells->InsertCellPoint(1);

  double scale = this->Scale;

  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(8);
  points->InsertPoint(0, scale * 0,     scale * 0,    0);
  points->InsertPoint(1, scale * 0.2,   scale * 0.2,  0);
  points->InsertPoint(2, scale * 0,     scale * 1,    0);
  points->InsertPoint(3, scale * -0.2,  scale * 0.2,  0);
  points->InsertPoint(4, scale * -1,    scale * 0,    0);
  points->InsertPoint(5, scale * -0.2,  scale * -0.2, 0);
  points->InsertPoint(6, scale * 0,     scale * -1,   0);
  points->InsertPoint(7, scale * 0.2,   scale * -0.2, 0);
  points->InsertPoint(8, scale * 1,     scale * 0,    0);

  vtkPolyData* polyData = vtkPolyData::New();
  polyData->SetPolys(cells);
  polyData->SetPoints(points);
  cells->Delete();
  points->Delete();
  return polyData;
}

void vtkCustomWidget::CreateTranslationHandle()
{
  vtkDataSetMapper* mapper = vtkDataSetMapper::New();
  vtkPolyData* polyData = this->CreateTranslationHandlePolydata();
  mapper->SetInputData(polyData);
  polyData->Delete();

  vtkProperty* orangeProperty = vtkProperty::New();
  orangeProperty->SetColor(1.0, 0.5, 0.0);
  orangeProperty->SetAmbient(1);
  orangeProperty->SetDiffuse(0);
  this->TranslationHandle->SetMapper(mapper);
  this->TranslationHandle->SetProperty(orangeProperty);
  this->TranslationHandle->SetPosition(this->Center);
  mapper->Delete();
  orangeProperty->Delete();
}

void vtkCustomWidget::SetScale(double value)
{
  this->Scale = value;
  if (this->TranslationHandle != NULL) {
    vtkDataSetMapper* mapper = vtkDataSetMapper::New();
    vtkPolyData* polyData = this->CreateTranslationHandlePolydata();
    mapper->SetInputData(polyData);
    polyData->Delete();
    this->TranslationHandle->SetMapper(mapper);
  }
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
  this->HandlePicker->GetPickPosition(this->StartPick);
  this->CopyVector(this->Center, this->StartCenter);
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

  double x = (double)this->Interactor->GetEventPosition()[0];
  double y = (double)this->Interactor->GetEventPosition()[1];

  double lastPos[4];
  this->ComputeWorldToDisplay(this->StartPick[0],
                              this->StartPick[1],
                              this->StartPick[2],
                              lastPos);

  // Compute the two points defining the motion
  double currentPickPosition[4];
  this->ComputeDisplayToWorld(x, y, lastPos[2], currentPickPosition);

  double offset[3];
  // Get the offset between the pick position when dragging started and the
  // current pick position.
  this->SubtractVectors(currentPickPosition, this->StartPick, offset);
  // Set the center to the center when dragging started plus the offset of
  // the current pick position. This is done so the dragged object is always
  // under the mouse and doesn't drift because of rounding errors.
  this->AddVectors(this->StartCenter, offset, this->Center);
  this->TranslationHandle->SetPosition(this->Center);

  this->EventCallbackCommand->SetAbortFlag(1);
  this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
  this->Interactor->Render();
}

void vtkCustomWidget::SubtractVectors(double vector1[3], double vector2[3],
                                      double* resultvector)
{
  for (int i = 0; i < 3; i++)
    resultvector[i] = vector1[i] - vector2[i];
}

void vtkCustomWidget::AddVectors(double vector1[3], double vector2[3],
                                 double* resultvector)
{
  for (int i = 0; i < 3; i++)
    resultvector[i] = vector1[i] + vector2[i];
}

void vtkCustomWidget::CopyVector(double inputVector[3], double outputVector[3])
{
  for (int i = 0; i < 3; i++)
    outputVector[i] = inputVector[i];
}
