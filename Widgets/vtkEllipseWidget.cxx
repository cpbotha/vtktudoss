/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEllipseWidget.cxx,v

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkEllipseWidget.h"
#include "vtkEllipseRepresentation.h"
#include "vtkCommand.h"
#include "vtkCallbackCommand.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkHandleWidget.h"
#include "vtkHandleRepresentation.h"
#include "vtkCoordinate.h"
#include "vtkWidgetCallbackMapper.h"
#include "vtkWidgetEvent.h"
#include "vtkRenderWindow.h"

vtkCxxRevisionMacro(vtkEllipseWidget, "1.10");
vtkStandardNewMacro(vtkEllipseWidget);


// The bidimensional widget observes the handles.
// Here we create the command/observer classes to respond to the 
// slider widgets.
class vtkEllipseWidgetCallback : public vtkCommand
{
public:
  static vtkEllipseWidgetCallback *New() 
    { return new vtkEllipseWidgetCallback; }
  virtual void Execute(vtkObject*, unsigned long eventId, void*)
    {
      switch (eventId)
        {
        case vtkCommand::StartInteractionEvent:
          this->EllipseWidget->StartEllipseInteraction();
          break;
        case vtkCommand::EndInteractionEvent:
          this->EllipseWidget->EndEllipseInteraction();
          break;
        }
    }
  vtkEllipseWidget *EllipseWidget;
};


//----------------------------------------------------------------------
vtkEllipseWidget::vtkEllipseWidget()
{
  this->ManagesCursor = 1;

  // Manage priorities, we want the handles to be lower priority
  if ( this->Priority <= 0.0 )
    {
    this->Priority = 0.01;
    }

  this->FirstClick = 1;

  // The widgets for moving the end points. They observe this widget (i.e.,
  // this widget is the parent to the handles).
  this->Point1Widget = vtkHandleWidget::New();
  this->Point1Widget->SetPriority(this->Priority-0.01);
  this->Point1Widget->SetParent(this);
  this->Point1Widget->ManagesCursorOff();

  this->Point2Widget = vtkHandleWidget::New();
  this->Point2Widget->SetPriority(this->Priority-0.01);
  this->Point2Widget->SetParent(this);
  this->Point2Widget->ManagesCursorOff();

  this->Point3Widget = vtkHandleWidget::New();
  this->Point3Widget->SetPriority(this->Priority-0.01);
  this->Point3Widget->SetParent(this);
  this->Point3Widget->ManagesCursorOff();

  this->Point4Widget = vtkHandleWidget::New();
  this->Point4Widget->SetPriority(this->Priority-0.01);
  this->Point4Widget->SetParent(this);
  this->Point4Widget->ManagesCursorOff();


  // Set up the callbacks on the two handles
  this->EllipseWidgetCallback1 = vtkEllipseWidgetCallback::New();
  this->EllipseWidgetCallback1->EllipseWidget = this;
  this->Point1Widget->AddObserver(vtkCommand::StartInteractionEvent, this->EllipseWidgetCallback1, 
                                  this->Priority);
  this->Point1Widget->AddObserver(vtkCommand::EndInteractionEvent, this->EllipseWidgetCallback1,
                                  this->Priority);

  this->EllipseWidgetCallback2 = vtkEllipseWidgetCallback::New();
  this->EllipseWidgetCallback2->EllipseWidget = this;
  this->Point2Widget->AddObserver(vtkCommand::StartInteractionEvent, this->EllipseWidgetCallback2, 
                                  this->Priority);
  this->Point2Widget->AddObserver(vtkCommand::EndInteractionEvent, this->EllipseWidgetCallback2,
                                  this->Priority);


  this->EllipseWidgetCallback3 = vtkEllipseWidgetCallback::New();
  this->EllipseWidgetCallback3->EllipseWidget = this;
  this->Point3Widget->AddObserver(vtkCommand::StartInteractionEvent, this->EllipseWidgetCallback3, 
                                  this->Priority);
  this->Point3Widget->AddObserver(vtkCommand::EndInteractionEvent, this->EllipseWidgetCallback3,
                                  this->Priority);


  this->EllipseWidgetCallback4 = vtkEllipseWidgetCallback::New();
  this->EllipseWidgetCallback4->EllipseWidget = this;
  this->Point4Widget->AddObserver(vtkCommand::StartInteractionEvent, this->EllipseWidgetCallback4, 
                                  this->Priority);
  this->Point4Widget->AddObserver(vtkCommand::EndInteractionEvent, this->EllipseWidgetCallback4,
                                  this->Priority);


  // These are the event callbacks supported by this widget
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonPressEvent,
                                          vtkWidgetEvent::AddPoint,
                                          this, vtkEllipseWidget::AddPointAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::MouseMoveEvent,
                                          vtkWidgetEvent::Move,
                                          this, vtkEllipseWidget::MoveAction);
  this->CallbackMapper->SetCallbackMethod(vtkCommand::LeftButtonReleaseEvent,
                                          vtkWidgetEvent::EndSelect,
                                          this, vtkEllipseWidget::EndSelectAction);

}

//----------------------------------------------------------------------
vtkEllipseWidget::~vtkEllipseWidget()
{
  this->Point1Widget->RemoveObserver(this->EllipseWidgetCallback1);
  this->Point1Widget->Delete();
  this->EllipseWidgetCallback1->Delete();

  this->Point2Widget->RemoveObserver(this->EllipseWidgetCallback2);
  this->Point2Widget->Delete();
  this->EllipseWidgetCallback2->Delete();

  this->Point3Widget->RemoveObserver(this->EllipseWidgetCallback3);
  this->Point3Widget->Delete();
  this->EllipseWidgetCallback3->Delete();

  this->Point4Widget->RemoveObserver(this->EllipseWidgetCallback4);
  this->Point4Widget->Delete();
  this->EllipseWidgetCallback4->Delete();
}

//----------------------------------------------------------------------
void vtkEllipseWidget::CreateDefaultRepresentation()
{
  if ( ! this->WidgetRep )
    {
    this->WidgetRep = vtkEllipseRepresentation::New();
    }
  vtkEllipseRepresentation::SafeDownCast(this->WidgetRep)->
    InstantiateHandleRepresentation();
}

//----------------------------------------------------------------------
void vtkEllipseWidget::SetEnabled(int enabling)
{
  // The handle widgets are not actually enabled until they are placed.
  // The handle widgets take their representation from the vtkEllipseRepresentation.
  if ( enabling )
    {
    if ( this->FirstClick )
      {
      if (this->WidgetRep)
        {
        vtkEllipseRepresentation::SafeDownCast(this->WidgetRep)->
          Line1VisibilityOff();
        vtkEllipseRepresentation::SafeDownCast(this->WidgetRep)->
          Line2VisibilityOff();
        }
      } // this->FirstClick
    else
      {
      if (this->WidgetRep)
        {
        vtkEllipseRepresentation::SafeDownCast(this->WidgetRep)->
          Line1VisibilityOn();
        vtkEllipseRepresentation::SafeDownCast(this->WidgetRep)->
          Line2VisibilityOn();
        }
      if (this->Point1Widget)
        {
        this->Point1Widget->SetEnabled(1);
        }
      if (this->Point2Widget)
        {
        this->Point2Widget->SetEnabled(1);
        }
      if (this->Point3Widget)
        {
        this->Point3Widget->SetEnabled(1);
        }
      if (this->Point4Widget)
        {
        this->Point4Widget->SetEnabled(1);
        }
      }
    }

  // Done in this wierd order to get everything to work right. This invocation creates the
  // default representation.
  this->Superclass::SetEnabled(enabling);

  if ( enabling )
    {
      if (this->Point1Widget)
        {
        this->Point1Widget->SetRepresentation(
          vtkEllipseRepresentation::SafeDownCast
          (this->WidgetRep)->GetPoint1Representation());
        this->Point1Widget->SetInteractor(this->Interactor);
        this->Point1Widget->GetRepresentation()->SetRenderer(
          this->CurrentRenderer);
        }
      if (this->Point2Widget)
        {
        this->Point2Widget->SetRepresentation(
          vtkEllipseRepresentation::SafeDownCast
          (this->WidgetRep)->GetPoint2Representation());
        this->Point2Widget->SetInteractor(this->Interactor);
        this->Point2Widget->GetRepresentation()->SetRenderer(
          this->CurrentRenderer);
        }
      if (this->Point3Widget)
        {
        this->Point3Widget->SetRepresentation(
          vtkEllipseRepresentation::SafeDownCast
          (this->WidgetRep)->GetPoint3Representation());
        this->Point3Widget->SetInteractor(this->Interactor);
        this->Point3Widget->GetRepresentation()->SetRenderer(
          this->CurrentRenderer);
        }
      if (this->Point4Widget)
        {
        this->Point4Widget->SetRepresentation(
          vtkEllipseRepresentation::SafeDownCast
          (this->WidgetRep)->GetPoint4Representation());
        this->Point4Widget->SetInteractor(this->Interactor);
        this->Point4Widget->GetRepresentation()->SetRenderer(
          this->CurrentRenderer);
        }
    }
  else //disabling widget
    {
    if (this->Point1Widget)
      {
      this->Point1Widget->SetEnabled(0);
      }
    if (this->Point2Widget)
      {
      this->Point2Widget->SetEnabled(0);
      }
    if (this->Point3Widget)
      {
      this->Point3Widget->SetEnabled(0);
      }
    if (this->Point4Widget)
      {
      this->Point4Widget->SetEnabled(0);
      }
    }
}

//----------------------------------------------------------------------
int vtkEllipseWidget::IsMeasureValid()
{
  // we're probably always valid. (durr)
  return 1;    
}

// The following methods are the callbacks that the bidimensional widget responds to. 
//-------------------------------------------------------------------------
void vtkEllipseWidget::AddPointAction(vtkAbstractWidget *w)
{
  vtkEllipseWidget *self = vtkEllipseWidget::SafeDownCast(w);
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);

  if (self->FirstClick)
    {
    self->GrabFocus(self->EventCallbackCommand);
    self->FirstClick = 0;
    self->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);

    self->Point1Widget->SetEnabled(1);
	self->Point2Widget->SetEnabled(1);
	self->Point3Widget->SetEnabled(1);
	self->Point4Widget->SetEnabled(1);

    vtkEllipseRepresentation::SafeDownCast(self->WidgetRep)->StartWidgetDefinition(e);
    //self->InvokeEvent(vtkCommand::PlacePointEvent,(void*)&(self->CurrentHandle));
    vtkEllipseRepresentation::SafeDownCast(self->WidgetRep)->Line1VisibilityOn();
    vtkEllipseRepresentation::SafeDownCast(self->WidgetRep)->Line2VisibilityOn();
    }

  else // !this->FirstClick
    {

	self->HandleLine1Selected = 0;
	self->HandleLine2Selected = 0;
	self->Line1Selected = 0;
	self->Line2Selected = 0;
	self->CenterSelected = 0;
	int modifier = self->Interactor->GetShiftKey() | self->Interactor->GetControlKey();
	int state = self->WidgetRep->ComputeInteractionState(X,Y,modifier);
	if ( state == vtkEllipseRepresentation::Outside )
	  {
	  return;
      }

	self->GrabFocus(self->EventCallbackCommand);
	vtkEllipseRepresentation::SafeDownCast(self->WidgetRep)->StartWidgetManipulation(e);
	if ( state == vtkEllipseRepresentation::NearP1 ||
		 state == vtkEllipseRepresentation::NearP2 )
      {
      self->HandleLine1Selected = 1;
      self->InvokeEvent(vtkCommand::LeftButtonPressEvent,NULL);
      }
    else if ( state == vtkEllipseRepresentation::NearP3 ||
              state == vtkEllipseRepresentation::NearP4 )
      {
      self->HandleLine2Selected = 1;
      self->InvokeEvent(vtkCommand::LeftButtonPressEvent,NULL);
      }
    else if ( state == vtkEllipseRepresentation::OnL1)
      {
      self->WidgetRep->Highlight(1);
      self->Line1Selected = 1;
      self->StartEllipseInteraction();
      }
    else if ( state == vtkEllipseRepresentation::OnL2)
      {
      self->WidgetRep->Highlight(1);
      self->Line2Selected = 1;
      self->StartEllipseInteraction();
      }
    else if ( state == vtkEllipseRepresentation::OnCenter )
      {
      self->WidgetRep->Highlight(1);
      self->CenterSelected = 1;
      self->StartEllipseInteraction();
      }

    } // !this->FirstClick
  
  self->EventCallbackCommand->SetAbortFlag(1);
  self->Render();
}

//-------------------------------------------------------------------------
void vtkEllipseWidget::MoveAction(vtkAbstractWidget *w)
{
  vtkEllipseWidget *self = vtkEllipseWidget::SafeDownCast(w);

  // Do nothing if outside
  if ( self->FirstClick )
    {
    return;
    }

  // Delegate the event consistent with the state
  int X = self->Interactor->GetEventPosition()[0];
  int Y = self->Interactor->GetEventPosition()[1];
  double e[2];
  e[0] = static_cast<double>(X);
  e[1] = static_cast<double>(Y);
  double p1[3], p2[3], slope;

  if ( self->Line1Selected || self->Line2Selected )
    {
    // moving outer portion of line -- rotating
    self->RequestCursorShape(VTK_CURSOR_HAND);
    vtkEllipseRepresentation::SafeDownCast(self->WidgetRep)->
      WidgetInteraction(e);
    self->InvokeEvent(vtkCommand::InteractionEvent,NULL);
    }
  else if ( self->HandleLine1Selected )
    { // moving one of the endpoints of line 1
    reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
      GetPoint1DisplayPosition(p1);
    reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
      GetPoint2DisplayPosition(p2);
    slope = VTK_DOUBLE_MAX;
    if (p1[0] != p2[0])
      {
      slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);
      }
    if ( slope > -1 && slope < 1)
      {
      self->RequestCursorShape(VTK_CURSOR_SIZEWE);
      }
    else
      {
      self->RequestCursorShape(VTK_CURSOR_SIZENS);
      }

    reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
      WidgetInteraction(e);
    self->InvokeEvent(vtkCommand::InteractionEvent,NULL);
    }

  else if ( self->HandleLine2Selected )
    { // moving one of the endpoints of line 2
    reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
      GetPoint3DisplayPosition(p1);
    reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
      GetPoint4DisplayPosition(p2);
    slope = VTK_DOUBLE_MAX;
    if (p1[0] != p2[0])
      {
      slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);
      }
    if ( slope > -1 && slope < 1)
      {
      self->RequestCursorShape(VTK_CURSOR_SIZEWE);
      }
    else
      {
      self->RequestCursorShape(VTK_CURSOR_SIZENS);
      }

    reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
      WidgetInteraction(e);
    self->InvokeEvent(vtkCommand::InteractionEvent,NULL);
    }

  else if ( self->CenterSelected )
    {//grabbing center intersection point
    self->RequestCursorShape(VTK_CURSOR_SIZEALL);
    vtkEllipseRepresentation::SafeDownCast(self->WidgetRep)->
      WidgetInteraction(e);
    self->InvokeEvent(vtkCommand::InteractionEvent,NULL);
    }

  else // just moving around, nothing yet selected
    {
    int state = self->WidgetRep->ComputeInteractionState(X,Y);
    if ( state == vtkEllipseRepresentation::Outside )
      {
      self->RequestCursorShape(VTK_CURSOR_DEFAULT);
      }
    else if ( state == vtkEllipseRepresentation::OnCenter )
      {
      self->RequestCursorShape(VTK_CURSOR_SIZEALL);
      }
    else if ( state == vtkEllipseRepresentation::NearP1 ||
              state == vtkEllipseRepresentation::NearP2 )
      {
      reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
        GetPoint1DisplayPosition(p1);
      reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
        GetPoint2DisplayPosition(p2);
      slope = VTK_DOUBLE_MAX;
      if (p1[0] != p2[0])
        {
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);
        }
      if ( slope > -1 && slope < 1)
        {
        self->RequestCursorShape(VTK_CURSOR_SIZEWE);
        }
      else
        {
        self->RequestCursorShape(VTK_CURSOR_SIZENS);
        }
      }
    else if ( state == vtkEllipseRepresentation::NearP3 ||
              state == vtkEllipseRepresentation::NearP4 )
      {
      reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
        GetPoint3DisplayPosition(p1);
      reinterpret_cast<vtkEllipseRepresentation*>(self->WidgetRep)->
        GetPoint4DisplayPosition(p2);
      slope = VTK_DOUBLE_MAX;
      if (p1[0] != p2[0])
        {
        slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);
        }
      if ( slope > -1 && slope < 1)
        {
        self->RequestCursorShape(VTK_CURSOR_SIZEWE);
        }
      else
        {
        self->RequestCursorShape(VTK_CURSOR_SIZENS);
        }
      }
    else
      {
      self->RequestCursorShape(VTK_CURSOR_HAND);
      }
    }

  self->WidgetRep->BuildRepresentation();
  self->Render();
}

//-------------------------------------------------------------------------
void vtkEllipseWidget::EndSelectAction(vtkAbstractWidget *w)
{
  vtkEllipseWidget *self = vtkEllipseWidget::SafeDownCast(w);

  self->Line1Selected = 0;
  self->Line2Selected = 0;
  self->HandleLine1Selected = 0;
  self->HandleLine2Selected = 0;
  self->CenterSelected = 0;
  self->WidgetRep->Highlight(0);
  self->ReleaseFocus();
  self->WidgetRep->BuildRepresentation();
  int state = self->WidgetRep->GetInteractionState();
  if ( state == vtkEllipseRepresentation::NearP1 ||
       state == vtkEllipseRepresentation::NearP2 ||
       state == vtkEllipseRepresentation::NearP3 ||
       state == vtkEllipseRepresentation::NearP4 )
    {
    self->InvokeEvent(vtkCommand::LeftButtonReleaseEvent,NULL);
    }
  else
    {
    self->EndEllipseInteraction();
    }
  self->EventCallbackCommand->SetAbortFlag(1);
  self->Render();
}

// These are callbacks that are active when the user is manipulating the
// handles of the angle widget.
//----------------------------------------------------------------------
void vtkEllipseWidget::StartEllipseInteraction()
{
  this->Superclass::StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent,NULL);
}

//----------------------------------------------------------------------
void vtkEllipseWidget::EndEllipseInteraction()
{
  this->Superclass::EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent,NULL);
}

//----------------------------------------------------------------------
void vtkEllipseWidget::SetProcessEvents(int pe)
{
  this->Superclass::SetProcessEvents(pe);

  this->Point1Widget->SetProcessEvents(pe);
  this->Point2Widget->SetProcessEvents(pe);
  this->Point3Widget->SetProcessEvents(pe);
  this->Point4Widget->SetProcessEvents(pe);
}

//----------------------------------------------------------------------
void vtkEllipseWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  //Superclass typedef defined in vtkTypeMacro() found in vtkSetGet.h
  this->Superclass::PrintSelf(os,indent);
}
