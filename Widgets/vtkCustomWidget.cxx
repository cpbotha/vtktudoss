#include "vtkCustomWidget.h"

#include "vtkBoxRepresentation.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkObjectFactory.h" // Remember to always include this one, otherwise you get a error in vtkStandardNewMacro
#include "vtkWidgetEventTranslator.h"
#include "vtkWidgetCallbackMapper.h"
#include "vtkEvent.h"
#include "vtkWidgetEvent.h"
#include "vtkRenderWindow.h"


vtkStandardNewMacro(vtkCustomWidget);

// Constructor
vtkCustomWidget::vtkCustomWidget()
{
  this->State = vtkCustomWidget::Start;
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
  std::cout << "Create default representation" << std::endl;
  if (!this->WidgetRep)
  {
    this->WidgetRep = vtkBoxRepresentation::New();
  }
}

void vtkCustomWidget::SelectAction(vtkAbstractWidget* w)
{
  std::cout << "SelectAction" << std::endl;

  // We are in a static method, cast to ourself
  vtkCustomWidget* self = reinterpret_cast<vtkCustomWidget*>(w);

  self->GrabFocus(self->EventCallbackCommand);

  // start the interaction
  self->EventCallbackCommand->SetAbortFlag(1);
  self->StartInteraction();
  self->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
  self->Render();
}

void vtkCustomWidget::EndSelectAction(vtkAbstractWidget* w)
{
  std::cout << "EndSelectAction" << std::endl;
}

void vtkCustomWidget::MoveAction(vtkAbstractWidget* w)
{
  std::cout << "MoveAction" << std::endl;
}

void vtkCustomWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}