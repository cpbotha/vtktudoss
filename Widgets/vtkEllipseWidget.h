/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEllipseWidget.h,v

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkEllipseWidget - interact with a 2D ellipse somewhere
// .SECTION Description

// widget to place 2D ellipse.  also instantiate a vtkEllipseRepresentation 
// to make use of this widget.  Oh hell I hate the new widget architecture.
// Complexity has increased ten-fold, very little benefit otherwise.

// .SECTION Event Bindings
// By default, the widget responds to the following VTK events (i.e., it
// watches the vtkRenderWindowInteractor for these events):
// <pre>
//   LeftButtonPressEvent - define a point or manipulate a handle, line,
//                          perform rotation or translate the widget.
//   MouseMoveEvent - position the points, move a line, rotate or translate the widget
//   LeftButtonReleaseEvent - release the selected handle and end interaction
// </pre>
//
// Note that the event bindings described above can be changed using this
// class's vtkWidgetEventTranslator. This class translates VTK events 
// into the vtkEllipseWidget's widget events:
// <pre>
//   vtkWidgetEvent::AddPoint -- (In Define mode:) Add one point; depending on the 
//                               state it may the first, second, third or fourth 
//                               point added. (In Manipulate mode:) If near a handle, 
//                               select the handle. Or if near a line, select the line.
//   vtkWidgetEvent::Move -- (In Define mode:) Position the second, third or fourth 
//                           point. (In Manipulate mode:) Move the handle, line or widget.
//   vtkWidgetEvent::EndSelect -- the manipulation process has completed.
// </pre>
//
// This widget invokes the following VTK events on itself (which observers
// can listen for):
// <pre>
//   vtkCommand::StartInteractionEvent (beginning to interact)
//   vtkCommand::EndInteractionEvent (completing interaction)
//   vtkCommand::InteractionEvent (moving a handle, line or performing rotation)
//   vtkCommand::PlacePointEvent (after a point is positioned; 
//                                call data includes handle id (0,1,2,4))
// </pre>

// .SECTION Thanks
// Charl P. Botha c.p.botha@tudelft.nl for creating and contributing
// this class.
//


// .SECTION See Also
// vtkHandleWidget vtkDistanceWidget


#ifndef __vtkEllipseWidget_h
#define __vtkEllipseWidget_h

#include "vtkAbstractWidget.h"

class vtkEllipseRepresentation;
class vtkHandleWidget;
class vtkEllipseWidgetCallback;


class VTK_EXPORT vtkEllipseWidget : public vtkAbstractWidget
{
public:
  // Description:
  // Instantiate this class.
  static vtkEllipseWidget *New();

  // Description:
  // Standard methods for a VTK class.
  vtkTypeRevisionMacro(vtkEllipseWidget,vtkAbstractWidget);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // The method for activiating and deactiviating this widget. This method
  // must be overridden because it is a composite widget and does more than
  // its superclasses' vtkAbstractWidget::SetEnabled() method.
  virtual void SetEnabled(int);

  // Description:
  // Specify an instance of vtkWidgetRepresentation used to represent this
  // widget in the scene. Note that the representation is a subclass of vtkProp
  // so it can be added to the renderer independent of the widget.
  void SetRepresentation(vtkEllipseRepresentation *r)
    {this->Superclass::SetWidgetRepresentation(reinterpret_cast<vtkWidgetRepresentation*>(r));}
  
  // Description:
  // Create the default widget representation if one is not set. 
  void CreateDefaultRepresentation();

  // Description:
  // A flag indicates whether the bi-dimensional measure is valid. The widget
  // becomes valid after two of the four points are placed.
  int IsMeasureValid();

  // Description:
  // Events. 
  //BTX
  enum
  {
  EndWidgetSelectEvent = 10050
  };
  //ETX

  // Description:
  // Methods to change the whether the widget responds to interaction.
  // Overridden to pass the state to component widgets.
  virtual void SetProcessEvents(int);


protected:
  vtkEllipseWidget();
  ~vtkEllipseWidget();

  int FirstClick;

  int HandleLine1Selected;
  int HandleLine2Selected;
  int Line1Selected;
  int Line2Selected;
  int CenterSelected;

  // Callback interface to capture events when
  // placing the widget.
  static void AddPointAction(vtkAbstractWidget*);
  static void MoveAction(vtkAbstractWidget*);
  static void EndSelectAction(vtkAbstractWidget*);
  
  // The positioning handle widgets
  vtkHandleWidget *Point1Widget;
  vtkHandleWidget *Point2Widget;
  vtkHandleWidget *Point3Widget;
  vtkHandleWidget *Point4Widget;
  vtkEllipseWidgetCallback *EllipseWidgetCallback1;
  vtkEllipseWidgetCallback *EllipseWidgetCallback2;
  vtkEllipseWidgetCallback *EllipseWidgetCallback3;
  vtkEllipseWidgetCallback *EllipseWidgetCallback4;
  
  // Methods invoked when the handles at the
  // end points of the widget are manipulated
  void StartEllipseInteraction();
  virtual void EndEllipseInteraction();
  
//BTX
  friend class vtkEllipseWidgetCallback;
//ETX  

private:
  vtkEllipseWidget(const vtkEllipseWidget&);  //Not implemented
  void operator=(const vtkEllipseWidget&);  //Not implemented
};

#endif
