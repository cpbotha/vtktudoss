#include "vtkProsthesisWidget.h"

#include "vtkProsthesisRepresentation.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkObjectFactory.h" // Remember to always include this one, otherwise you get a error in vtkStandardNewMacro
#include "vtkWidgetCallbackMapper.h"
#include "vtkEvent.h"
#include "vtkWidgetEvent.h"
#include "vtkRenderer.h"


vtkStandardNewMacro(vtkProsthesisWidget);

// Constructor
vtkProsthesisWidget::vtkProsthesisWidget()
{
  this->WidgetState = vtkProsthesisWidget::Start;
  // Define widget events
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonPressEvent,
                                          vtkEvent::NoModifier,
                                          0, 0, NULL,
                                          vtkWidgetEvent::Select,
                                          this, vtkProsthesisWidget::SelectAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonReleaseEvent,
                                          vtkEvent::NoModifier,
                                          0, 0, NULL,
                                          vtkWidgetEvent::EndSelect,
                                          this, vtkProsthesisWidget::EndSelectAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                          vtkWidgetEvent::Move,
                                          this, vtkProsthesisWidget::MoveAction);
}

// Destructor
vtkProsthesisWidget::~vtkProsthesisWidget()
{
}

void vtkProsthesisWidget::CreateDefaultRepresentation()
{
  if (!this->WidgetRep)
  {
    this->WidgetRep = vtkProsthesisRepresentation::New();
  }
}

void vtkProsthesisWidget::SelectAction(vtkAbstractWidget* w)
{
  // We are in a static method, cast to ourself
  vtkProsthesisWidget* self = reinterpret_cast<vtkProsthesisWidget*>(w);

  // Get the event position
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!self->CurrentRenderer || !self->CurrentRenderer->IsInViewport(X,Y))
  {
    self->WidgetState = vtkProsthesisWidget::Start;
    return;
  }

  self->WidgetState = vtkProsthesisWidget::Active;
  // Begin the widget interaction which has the side effect of setting the
  // interaction state.
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->StartWidgetInteraction(e);
  int interactionState = self->WidgetRep->GetInteractionState();
  if (interactionState == vtkProsthesisRepresentation::Outside)
  {
    return;
  }
  else if (interactionState == vtkProsthesisRepresentation::UpClick)
  {
    self->InvokeEvent(vtkCommand::UserEvent, NULL);
  }
  else if (interactionState == vtkProsthesisRepresentation::DownClick)
  {
    self->InvokeEvent(vtkCommand::UserEvent + 1, NULL);
  }

  self->GrabFocus(self->EventCallbackCommand);

  // The SetInteractionState has the side effect of highlighting the widget
  reinterpret_cast<vtkProsthesisRepresentation*>(self->WidgetRep)->
    SetInteractionState(interactionState);

  // start the interaction
  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
  self->Render();
}

void vtkProsthesisWidget::EndSelectAction(vtkAbstractWidget* w)
{
  vtkProsthesisWidget* self = reinterpret_cast<vtkProsthesisWidget*>(w);

  // If interaction has started, the widget will be in the active state
  if (self->WidgetState == vtkProsthesisWidget::Start)
  {
    return;
  }

  // Return the widget and representation's states to their defaults
  self->WidgetState = vtkProsthesisWidget::Start;
  reinterpret_cast<vtkProsthesisRepresentation*>(self->WidgetRep)->
    SetInteractionState(vtkProsthesisRepresentation::Outside);
  self->ReleaseFocus();

  self->EventCallbackCommand->SetAbortFlag(1);
  self->EndInteraction();
  self->InvokeEvent(vtkCommand::EndInteractionEvent,NULL);
  self->Render();
}

void vtkProsthesisWidget::MoveAction(vtkAbstractWidget* w)
{
  vtkProsthesisWidget* self = reinterpret_cast<vtkProsthesisWidget*>(w);

  // See whether we're active
  if (self->WidgetState == vtkProsthesisWidget::Start)
  {
    return;
  }

  // compute some info we need for all cases
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, adjust the representation
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->WidgetInteraction(e);

  // moving something
  self->EventCallbackCommand->SetAbortFlag(1);
  self->InvokeEvent(vtkCommand::InteractionEvent,NULL);
  self->Render();
}

void vtkProsthesisWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}