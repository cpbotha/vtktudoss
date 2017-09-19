/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkEllipseRepresentation.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkEllipseRepresentation - represent the vtkEllipseWidget
// .SECTION Description
// Use this for the vtkEllipseWidget.

// .SECTION Thanks
// Charl P. Botha c.p.botha@tudelft.nl for creating and contributing
// this class.
//

// .SECTION See Also
// vtkAngleWidget vtkHandleRepresentation vtkEllipseRepresentation


#ifndef __vtkEllipseRepresentation_h
#define __vtkEllipseRepresentation_h

#include "vtkWidgetRepresentation.h"

class vtkHandleRepresentation;
class vtkCellArray;
class vtkEllipseSource;
class vtkPoints;
class vtkPolyData;
class vtkPolyDataMapper2D;
class vtkPolyDataMapper;
class vtkTextMapper;
class vtkTransform;
class vtkTransformPolyDataFilter;
class vtkActor2D;
class vtkActor;
class vtkProperty2D;
class vtkTextProperty;


class VTK_EXPORT vtkEllipseRepresentation : public vtkWidgetRepresentation
{
public:
  // Description:
  // Instantiate the class.
  static vtkEllipseRepresentation *New();

  // Description:
  // Standard VTK methods.
  vtkTypeMacro(vtkEllipseRepresentation,vtkWidgetRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Methods to Set/Get the coordinates of the four points defining
  // this representation. Note that methods are available for both
  // display and world coordinates.
  virtual void SetPoint1WorldPosition(double pos[3]);
  virtual void SetPoint2WorldPosition(double pos[3]);
  virtual void SetPoint3WorldPosition(double pos[3]);
  virtual void SetPoint4WorldPosition(double pos[3]);
  virtual void GetPoint1WorldPosition(double pos[3]);
  virtual void GetPoint2WorldPosition(double pos[3]);
  virtual void GetPoint3WorldPosition(double pos[3]);
  virtual void GetPoint4WorldPosition(double pos[3]);
  virtual void SetPoint1DisplayPosition(double pos[3]);
  virtual void SetPoint2DisplayPosition(double pos[3]);
  virtual void SetPoint3DisplayPosition(double pos[3]);
  virtual void SetPoint4DisplayPosition(double pos[3]);
  virtual void GetPoint1DisplayPosition(double pos[3]);
  virtual void GetPoint2DisplayPosition(double pos[3]);
  virtual void GetPoint3DisplayPosition(double pos[3]);
  virtual void GetPoint4DisplayPosition(double pos[3]);

  // Description:
  // Method to get center in world coordinates of the current ellipse
  virtual void GetCenterWorldPosition(double pos[3]);

  // Description:
  // Return length of SemiMajor and Minor axes of the described ellipse
  virtual double GetSemiMajorAxisLength();
  virtual double GetSemiMinorAxisLength();

  // Description:
  // Return vector describing the relevant half-axis
  virtual void GetSemiMajorAxisVector(double v[3]);
  virtual void GetSemiMinorAxisVector(double v[3]);
  
  // Description:
  // Special methods for turning off the lines that define the bi-dimensional
  // measure. Generally these methods are used by the vtkEllipseWidget to
  // control the appearance of the widget. Note: turning off Line1 actually turns
  // off Line1 and Line2.
  vtkSetMacro(Line1Visibility,int);
  vtkGetMacro(Line1Visibility,int);
  vtkBooleanMacro(Line1Visibility,int);
  vtkSetMacro(Line2Visibility,int);
  vtkGetMacro(Line2Visibility,int);
  vtkBooleanMacro(Line2Visibility,int);

  // Description:
  // This method is used to specify the type of handle representation to use
  // for the four internal vtkHandleRepresentations within
  // vtkEllipseRepresentation.  To use this method, create a dummy
  // vtkHandleRepresentation (or subclass), and then invoke this method with
  // this dummy. Then the vtkEllipseRepresentation uses this dummy to
  // clone four vtkHandleRepresentations of the same type. Make sure you set the
  // handle representation before the widget is enabled. (The method
  // InstantiateHandleRepresentation() is invoked by the vtkEllipseWidget
  // for the purposes of cloning.)
  void SetHandleRepresentation(vtkHandleRepresentation *handle);
  void InstantiateHandleRepresentation();

  // Description:
  // Set/Get the handle representations used within the
  // vtkEllipseRepresentation. (Note: properties can be set by
  // grabbing these representations and setting the properties
  // appropriately.)
  vtkGetObjectMacro(Point1Representation,vtkHandleRepresentation);
  vtkGetObjectMacro(Point2Representation,vtkHandleRepresentation);
  vtkGetObjectMacro(Point3Representation,vtkHandleRepresentation);
  vtkGetObjectMacro(Point4Representation,vtkHandleRepresentation);

  // Description:
  // Retrieve the property used to control the appearance of the two
  // orthogonal lines.
  vtkGetObjectMacro(LineProperty,vtkProperty2D);
  vtkGetObjectMacro(SelectedLineProperty,vtkProperty2D);

  // Description:
  // Retrieve the property used to control the appearance of the text
  // labels.
  vtkGetObjectMacro(TextProperty,vtkTextProperty);

  // Description:
  // The tolerance representing the distance to the representation (in
  // pixels) in which the cursor is considered near enough to the
  // representation to be active.
  vtkSetClampMacro(Tolerance,int,1,100);
  vtkGetMacro(Tolerance,int);

  // Description:
  // Return the length of the line defined by (Point1,Point2). This is the
  // distance in the world coordinate system.
  virtual double GetLength1();

  // Description:
  // Return the length of the line defined by (Point3,Point4). This is the
  // distance in the world coordinate system.
  virtual double GetLength2();

  // Description:
  // Specify the format to use for labelling the distance. Note that an empty
  // string results in no label, or a format string without a "%" character
  // will not print the distance value.
  vtkSetStringMacro(LabelFormat);
  vtkGetStringMacro(LabelFormat);

//BTX -- used to communicate about the state of the representation
  enum {Outside=0,NearP1,NearP2,NearP3,NearP4,OnL1,OnL2,OnCenter};
//ETX

  // Description:
  // These are methods that satisfy vtkWidgetRepresentation's API.
  virtual void BuildRepresentation();
  virtual int ComputeInteractionState(int X, int Y, int modify=0);
  virtual void StartWidgetDefinition(double e[2]);
  virtual void Point2WidgetInteraction(double e[2]);
  virtual void Point3WidgetInteraction(double e[2]);
  virtual void StartWidgetManipulation(double e[2]);
  virtual void WidgetInteraction(double e[2]);
  virtual void Highlight(int highlightOn);

  // Description:
  // Methods required by vtkProp superclass.
  virtual void ReleaseGraphicsResources(vtkWindow *w);
  virtual int RenderOverlay(vtkViewport *viewport);

  // Description:
  // Toggle whether to display the label above or below the widget.
  // Defaults to 1
  vtkSetMacro(ShowLabelAboveWidget, int);
  vtkGetMacro(ShowLabelAboveWidget, int);

  // Description:
  // Set/get the id to display in the label.
  void SetID(unsigned long id);
  vtkGetMacro(ID, unsigned long);

  // Description:
  // Get the text shown in the widget's label.
  char* GetLabelText();

  // Description:
  // Get the position of the widget's label in display coordinates.
  double* GetLabelPosition();
  void GetLabelPosition(double pos[3]);

  vtkSetMacro(InitialSemiMajorAxisLength, double);
  vtkGetMacro(InitialSemiMajorAxisLength, double);
  vtkSetMacro(InitialSemiMinorAxisLength, double);
  vtkGetMacro(InitialSemiMinorAxisLength, double);

protected:
  vtkEllipseRepresentation();
  ~vtkEllipseRepresentation();

  // Keep track if modifier is set
  int Modifier;

  // The handle and the rep used to close the handles
  vtkHandleRepresentation *HandleRepresentation;
  vtkHandleRepresentation *Point1Representation;
  vtkHandleRepresentation *Point2Representation;
  vtkHandleRepresentation *Point3Representation;
  vtkHandleRepresentation *Point4Representation;

  // Selection tolerance for the handles
  int Tolerance;

  // Visibility of the lines
  int Line1Visibility;
  int Line2Visibility;

  // Geometry of the lines
  vtkCellArray        *LineCells;
  vtkPoints           *LinePoints;
  vtkPolyData         *LinePolyData;
  vtkPolyDataMapper2D *LineMapper;
  vtkActor2D          *LineActor;
  vtkProperty2D       *LineProperty;
  vtkProperty2D       *SelectedLineProperty;
  vtkEllipseSource    *EllipseSource;
  vtkPolyDataMapper2D *EllipseMapper;
  vtkActor2D          *EllipseActor;
  vtkTransform        *EllipseTrfm;
  vtkTransformPolyDataFilter *EllipseTrfmFilter;

  // The labels for the line lengths
  vtkTextProperty *TextProperty;
  vtkTextMapper   *TextMapper;
  vtkActor2D      *TextActor;

  unsigned long ID;
  int IDInitialized;

  double InitialSemiMajorAxisLength;
  double InitialSemiMinorAxisLength;

  // Internal variables
  double P1World[3];
  double P2World[3];
  double P3World[3];
  double P4World[3];
  double P21World[3];
  double P43World[3];
  double T21;
  double T43;
  double CenterWorld[3];
  double StartEventPositionWorld[4];

  // Format for printing the distance
  char *LabelFormat;

  // toggle to determine whether to place text above or below widget
  int ShowLabelAboveWidget;

  // Helper method
  void ProjectOrthogonalPoint(double x[4], double y[3],
			      double x1[3], double x2[3], double x21[3], 
                              double dir, double xP[3], double yP[3]);

  // given the four points defining the current crosshair thingy, 
  // update the internal ellipse representation
  void UpdateEllipse(double p1[3], double p2[3], 
	                 double p3[3], double p4[3]);


private:
  vtkEllipseRepresentation(const vtkEllipseRepresentation&);  //Not implemented
  void operator=(const vtkEllipseRepresentation&);  //Not implemented
};

#endif
