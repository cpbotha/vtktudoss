#include "vtkCustomWidget.h"

#include "vtkCustomRepresentation.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkObjectFactory.h" // Remember to always include this one, otherwise you get a error in vtkStandardNewMacro
#include "vtkWidgetCallbackMapper.h"
#include "vtkEvent.h"
#include "vtkWidgetEvent.h"
#include "vtkRenderer.h"


vtkStandardNewMacro(vtkCustomWidget);

// Constructor
vtkCustomWidget::vtkCustomWidget()
{
  this->WidgetState = vtkCustomWidget::Start;
  // Define widget events
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonPressEvent,
                                          vtkEvent::NoModifier,
                                          0, 0, NULL,
                                          vtkWidgetEvent::Select,
                                          this, vtkCustomWidget::SelectAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonReleaseEvent,
                                          vtkEvent::NoModifier,
                                          0, 0, NULL,
                                          vtkWidgetEvent::EndSelect,
                                          this, vtkCustomWidget::EndSelectAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                          vtkWidgetEvent::Move,
                                          this, vtkCustomWidget::MoveAction);
}

// Destructor
vtkCustomWidget::~vtkCustomWidget()
{
}

void vtkCustomWidget::CreateDefaultRepresentation()
{
  if (!this->WidgetRep)
  {
    this->WidgetRep = vtkCustomRepresentation::New();
  }
}

void vtkCustomWidget::SelectAction(vtkAbstractWidget* w)
{
  std::cout << "↓" << std::endl;

  // We are in a static method, cast to ourself
  vtkCustomWidget* self = reinterpret_cast<vtkCustomWidget*>(w);

  // Get the event position
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];

  // Okay, make sure that the pick is in the current renderer
  if (!self->CurrentRenderer || !self->CurrentRenderer->IsInViewport(X,Y))
  {
    self->WidgetState = vtkCustomWidget::Start;
    return;
  }

  self->WidgetState = vtkCustomWidget::Active;
  // Begin the widget interaction which has the side effect of setting the
  // interaction state.
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  self->WidgetRep->StartWidgetInteraction(e);
  int interactionState = self->WidgetRep->GetInteractionState();
  if (interactionState == vtkCustomRepresentation::Outside)
  {
    return;
  }

  self->GrabFocus(self->EventCallbackCommand);

  // The SetInteractionState has the side effect of highlighting the widget
  reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->
    SetInteractionState(interactionState);

  // start the interaction
  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
  self->Render();
}

void vtkCustomWidget::EndSelectAction(vtkAbstractWidget* w)
{
  std::cout << "↑" << std::endl;
  vtkCustomWidget* self = reinterpret_cast<vtkCustomWidget*>(w);

  // If interaction has started, the widget will be in the active state
  if (self->WidgetState == vtkCustomWidget::Start)
  {
    return;
  }

  // Return the widget and representation's states to their defaults
  self->WidgetState = vtkCustomWidget::Start;
  reinterpret_cast<vtkCustomRepresentation*>(self->WidgetRep)->
    SetInteractionState(vtkCustomRepresentation::Outside);
  self->ReleaseFocus();

  self->EventCallbackCommand->SetAbortFlag(1);
  self->EndInteraction();
  self->InvokeEvent(vtkCommand::EndInteractionEvent,NULL);
  self->Render();
}

void vtkCustomWidget::MoveAction(vtkAbstractWidget* w)
{
  std::cout << "↔";
  vtkCustomWidget* self = reinterpret_cast<vtkCustomWidget*>(w);

  // See whether we're active
  if (self->WidgetState == vtkCustomWidget::Start)
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

void vtkCustomWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}