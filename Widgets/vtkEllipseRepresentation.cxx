/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkEllipseRepresentation.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkEllipseRepresentation.h"
#include "vtkEllipseSource.h"
#include "vtkHandleRepresentation.h"
#include "vtkCoordinate.h"
#include "vtkRenderer.h"
#include "vtkMath.h"
#include "vtkLine.h"
#include "vtkTextProperty.h"
#include "vtkWindow.h"
#include "vtkCellArray.h"
#include "vtkCursor2D.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkActor2D.h"
#include "vtkTextMapper.h"
#include "vtkTextProperty.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkProperty2D.h"
#include "vtkPointHandleRepresentation2D.h"
#include "vtkObjectFactory.h"
#include "vtkInteractorObserver.h"


vtkStandardNewMacro(vtkEllipseRepresentation);


//----------------------------------------------------------------------
vtkEllipseRepresentation::vtkEllipseRepresentation()
{
  // By default, use one of these handles
  this->HandleRepresentation  = vtkPointHandleRepresentation2D::New();

  // why do I have to jump through all these hoops just to switch off
  // the frikking cursor??!

  // these lines taken from the vtkPointHandleRepresentation2D constructor,
  // but now with Axes off.
  vtkCursor2D *cursor2D = vtkCursor2D::New();
  cursor2D->AllOff();
  //cursor2D->AxesOn();
  cursor2D->PointOn();

  ((vtkPointHandleRepresentation2D *)(this->HandleRepresentation))->SetCursorShape(cursor2D->GetOutput());
  cursor2D->GetOutput()->Register(this->HandleRepresentation);

  // end of hoop-jumping ==========================================

  this->Point1Representation = NULL;
  this->Point2Representation = NULL;
  this->Point3Representation = NULL;
  this->Point4Representation = NULL;
  this->InstantiateHandleRepresentation();

  this->Modifier = 0;

  this->Tolerance = 5;
  this->Placed = 0;

  this->Line1Visibility = 1;
  this->Line2Visibility = 1;

  this->SetInitialSemiMajorAxisLength(5);
  this->SetInitialSemiMinorAxisLength(2.5);

  // Create the geometry for the two axes
  this->LineCells = vtkCellArray::New();
  this->LineCells->InsertNextCell(2);
  this->LineCells->InsertCellPoint(0);
  this->LineCells->InsertCellPoint(1);
  this->LineCells->InsertNextCell(2);
  this->LineCells->InsertCellPoint(2);
  this->LineCells->InsertCellPoint(3);
  this->LinePoints = vtkPoints::New();
  this->LinePoints->SetNumberOfPoints(4);
  this->LinePolyData = vtkPolyData::New();
  this->LinePolyData->SetPoints(this->LinePoints);
  this->LinePolyData->SetLines(this->LineCells);
  this->LineMapper = vtkPolyDataMapper2D::New();
  this->LineMapper->SetInputData(this->LinePolyData);
  this->LineProperty = vtkProperty2D::New();
  this->LineActor = vtkActor2D::New();
  this->LineActor->SetProperty(this->LineProperty);
  this->LineActor->SetMapper(this->LineMapper);
  this->SelectedLineProperty = vtkProperty2D::New();
  this->SelectedLineProperty->SetColor(0.0,1.0,0.0);
  this->SelectedLineProperty->SetLineWidth(2.0);

  this->TextProperty = vtkTextProperty::New();
  this->TextProperty->SetBold(1);
  this->TextProperty->SetItalic(1);
  this->TextProperty->SetShadow(1);
  this->TextProperty->SetFontFamilyToArial();
  this->TextMapper = vtkTextMapper::New();
  this->TextMapper->SetTextProperty(this->TextProperty);
  this->TextMapper->SetInput("0.0");
  this->TextActor = vtkActor2D::New();
  this->TextActor->SetMapper(this->TextMapper);
  this->TextActor->VisibilityOff();

  this->LabelFormat = new char[6];
  sprintf(this->LabelFormat,"%s","%0.3g");

  this->EllipseSource = vtkEllipseSource::New();
  this->EllipseSource->SetNumberOfSteps(64);
  this->EllipseTrfmFilter = vtkTransformPolyDataFilter::New();
  this->EllipseTrfm = vtkTransform::New();
  this->EllipseTrfmFilter->SetTransform(this->EllipseTrfm);
  this->EllipseTrfmFilter->SetInputData(this->EllipseSource->GetOutput());
  this->EllipseMapper = vtkPolyDataMapper2D::New();
  this->EllipseMapper->SetInputData(this->EllipseTrfmFilter->GetOutput());
  this->EllipseActor = vtkActor2D::New();
  this->EllipseActor->SetMapper(this->EllipseMapper);

  this->ID = 0;
  this->IDInitialized = 0;

  this->ShowLabelAboveWidget = 1;
}

//----------------------------------------------------------------------
vtkEllipseRepresentation::~vtkEllipseRepresentation()
{
  if ( this->HandleRepresentation )
    {
    this->HandleRepresentation->Delete();
    }
  if ( this->Point1Representation )
    {
    this->Point1Representation->Delete();
    }
  if ( this->Point2Representation )
    {
    this->Point2Representation->Delete();
    }
  if ( this->Point3Representation )
    {
    this->Point3Representation->Delete();
    }
  if ( this->Point4Representation )
    {
    this->Point4Representation->Delete();
    }

  this->LineCells->Delete();
  this->LinePoints->Delete();
  this->LinePolyData->Delete();
  this->LineMapper->Delete();
  this->LineProperty->Delete();
  this->LineActor->Delete();
  this->SelectedLineProperty->Delete();
  this->TextProperty->Delete();
  this->TextMapper->Delete();
  this->TextActor->Delete();
  this->SetLabelFormat(0);

  this->EllipseSource->Delete();
  this->EllipseMapper->Delete();
  this->EllipseActor->Delete();
  this->EllipseTrfm->Delete();
  this->EllipseTrfmFilter->Delete();
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation
::SetHandleRepresentation(vtkHandleRepresentation *handle)
{
  if ( handle == NULL || handle == this->HandleRepresentation )
    {
    return;
    }

  this->Modified();
  this->HandleRepresentation->Delete();
  this->HandleRepresentation = handle;
  this->HandleRepresentation->Register(this);

  this->Point1Representation->Delete();
  this->Point2Representation->Delete();
  this->Point3Representation->Delete();
  this->Point4Representation->Delete();

  this->Point1Representation = NULL;
  this->Point2Representation = NULL;
  this->Point3Representation = NULL;
  this->Point4Representation = NULL;

  this->InstantiateHandleRepresentation();
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint1WorldPosition(double pos[3])
{
  this->Point1Representation->GetWorldPosition(pos);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint2WorldPosition(double pos[3])
{
  this->Point2Representation->GetWorldPosition(pos);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint3WorldPosition(double pos[3])
{
  this->Point3Representation->GetWorldPosition(pos);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint4WorldPosition(double pos[3])
{
  this->Point4Representation->GetWorldPosition(pos);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint1DisplayPosition(double x[3])
{
  this->Point1Representation->SetDisplayPosition(x);
  double p[3];
  this->Point1Representation->GetWorldPosition(p);
  this->Point1Representation->SetWorldPosition(p);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint2DisplayPosition(double x[3])
{
  this->Point2Representation->SetDisplayPosition(x);
  double p[3];
  this->Point2Representation->GetWorldPosition(p);
  this->Point2Representation->SetWorldPosition(p);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint3DisplayPosition(double x[3])
{
  this->Point3Representation->SetDisplayPosition(x);
  double p[3];
  this->Point3Representation->GetWorldPosition(p);
  this->Point3Representation->SetWorldPosition(p);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint4DisplayPosition(double x[3])
{
  this->Point4Representation->SetDisplayPosition(x);
  double p[3];
  this->Point4Representation->GetWorldPosition(p);
  this->Point4Representation->SetWorldPosition(p);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint1WorldPosition(double x[3])
{
  this->Point1Representation->SetWorldPosition(x);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint2WorldPosition(double x[3])
{
  this->Point2Representation->SetWorldPosition(x);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint3WorldPosition(double x[3])
{
  this->Point3Representation->SetWorldPosition(x);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetPoint4WorldPosition(double x[3])
{
  this->Point4Representation->SetWorldPosition(x);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint1DisplayPosition(double pos[3])
{
  this->Point1Representation->GetDisplayPosition(pos);
  pos[2] = 0.0;
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint2DisplayPosition(double pos[3])
{
  this->Point2Representation->GetDisplayPosition(pos);
  pos[2] = 0.0;
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint3DisplayPosition(double pos[3])
{
  this->Point3Representation->GetDisplayPosition(pos);
  pos[2] = 0.0;
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetPoint4DisplayPosition(double pos[3])
{
  this->Point4Representation->GetDisplayPosition(pos);
  pos[2] = 0.0;
}


//----------------------------------------------------------------------
void vtkEllipseRepresentation::InstantiateHandleRepresentation()
{
  if ( ! this->Point1Representation )
    {
    this->Point1Representation = this->HandleRepresentation->NewInstance();
    this->Point1Representation->ShallowCopy(this->HandleRepresentation);
    }

  if ( ! this->Point2Representation )
    {
    this->Point2Representation = this->HandleRepresentation->NewInstance();
    this->Point2Representation->ShallowCopy(this->HandleRepresentation);
    }

  if ( ! this->Point3Representation )
    {
    this->Point3Representation = this->HandleRepresentation->NewInstance();
    this->Point3Representation->ShallowCopy(this->HandleRepresentation);
    }

  if ( ! this->Point4Representation )
    {
    this->Point4Representation = this->HandleRepresentation->NewInstance();
    this->Point4Representation->ShallowCopy(this->HandleRepresentation);
    }
}

//----------------------------------------------------------------------
int vtkEllipseRepresentation::ComputeInteractionState(int X, int Y, int modify)
{
  this->Modifier = modify;

  // See if we are near one of the end points or outside
  double pos1[3], pos2[3], pos3[3], pos4[3];
  this->GetPoint1DisplayPosition(pos1);
  this->GetPoint2DisplayPosition(pos2);
  this->GetPoint3DisplayPosition(pos3);
  this->GetPoint4DisplayPosition(pos4);

  double p1[3], p2[3], p3[3], p4[3], xyz[3];
  double t, closest[3];
  xyz[0] = static_cast<double>(X);
  xyz[1] = static_cast<double>(Y);
  p1[0] = static_cast<double>(pos1[0]);
  p1[1] = static_cast<double>(pos1[1]);
  p2[0] = static_cast<double>(pos2[0]);
  p2[1] = static_cast<double>(pos2[1]);
  p3[0] = static_cast<double>(pos3[0]);
  p3[1] = static_cast<double>(pos3[1]);
  p4[0] = static_cast<double>(pos4[0]);
  p4[1] = static_cast<double>(pos4[1]);
  xyz[2] = p1[2] = p2[2] = p3[2] = p4[2] = 0.0;

  double tol2 = this->Tolerance*this->Tolerance;
  // Check if we are on end points
  if ( vtkMath::Distance2BetweenPoints(xyz,p1) <= tol2 )
    {
    this->InteractionState = vtkEllipseRepresentation::NearP1;
    return this->InteractionState;
    }
  else if ( vtkMath::Distance2BetweenPoints(xyz,p2) <= tol2 )
    {
    this->InteractionState = vtkEllipseRepresentation::NearP2;
    return this->InteractionState;
    }
  else if ( vtkMath::Distance2BetweenPoints(xyz,p3) <= tol2 )
    {
    this->InteractionState = vtkEllipseRepresentation::NearP3;
    return this->InteractionState;
    }
  else if ( vtkMath::Distance2BetweenPoints(xyz,p4) <= tol2 )
    {
    this->InteractionState = vtkEllipseRepresentation::NearP4;
    return this->InteractionState;
    }

  // Compute intersection point.
  double uIntersect, vIntersect;
  vtkLine::Intersection(p1, p2, p3, p4, uIntersect, vIntersect);

  // Check if we are on edges
  int onL1 = (vtkLine::DistanceToLine(xyz,p1,p2,t,closest) <= tol2);
  int onL2 = (vtkLine::DistanceToLine(xyz,p3,p4,t,closest) <= tol2);

  //double xyzParam;

  if ( onL1 && onL2 )
    {
    this->InteractionState = vtkEllipseRepresentation::OnCenter;
    }
  else if ( onL1 )
    {
	this->InteractionState = vtkEllipseRepresentation::OnL1;
    }
  else if ( onL2 )
    {
	this->InteractionState = vtkEllipseRepresentation::OnL2;
    }
  else
    {
    this->InteractionState = vtkEllipseRepresentation::Outside;
    this->Modifier = 0;
    }

  return this->InteractionState;
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::StartWidgetDefinition(double e[2])
{
  // conversion from display to world via SetPoint?DisplayPos
  // ABSOLUTELY does not work (ginormous values for world,
  // although display is reasonable)

  // so now I'm doing an explicit conversion display->world using
  // the renderer.  This seems to work, and I'm too tired to figure
  // out why the BiDimensional way (via this->SetPoint?DisplayPos)
  // does not work.

  double c[3];
  c[0] = e[0];
  c[1] = e[1];
  c[2] = 0.0;

  if (this->Renderer)
    {
    this->Renderer->SetDisplayPoint(c);
    this->Renderer->DisplayToWorld();
    this->Renderer->GetWorldPoint(c);
    }

  //double *c = this->InitialCenter;
  double smaj = this->InitialSemiMajorAxisLength;
  double smin = this->InitialSemiMinorAxisLength;
  double p[3];

  p[0] = c[0] + smaj;
  p[1] = c[1];
  p[2] = c[2];
  this->SetPoint1WorldPosition(p);

  p[0] = c[0] - smaj;
  this->SetPoint2WorldPosition(p);

  p[0] = c[0];
  p[1] = c[1] - smin;
  this->SetPoint3WorldPosition(p);

  p[1] = c[1] + smin;
  this->SetPoint4WorldPosition(p);

  //this->SetPoint1DisplayPosition(pos);
  //this->SetPoint2DisplayPosition(pos);
  //this->SetPoint3DisplayPosition(pos);
  //this->SetPoint4DisplayPosition(pos);

  this->StartEventPosition[0] = c[0];
  this->StartEventPosition[1] = c[1];
  this->StartEventPosition[2] = c[2];
}

//----------------------------------------------------------------------


// cpbotha: both these methods (Point{2,3}WidgetInteraction are only
// called during widget definition.  We want to place this thing
// initially, so we'll throw these ugly things away later.
void vtkEllipseRepresentation::Point2WidgetInteraction(double e[2])
{
  double pos[3],p1[3];
  pos[0] = e[0];
  pos[1] = e[1];
  pos[2] = 0.0;

  // Make sure that the two points are not coincident
  // cpbotha: magic number 2 appearing here... oi.
  this->GetPoint1DisplayPosition(p1);
  if ( ((pos[0]-p1[0])*(pos[0]-p1[0]) + (pos[1]-p1[1])*(pos[1]-p1[1])) < 2 )
    {
    pos[0] += 2;
    }
  this->SetPoint2DisplayPosition(pos);
}

//----------------------------------------------------------------------
// This method is called when Point3 is to be manipulated. Note that Point3
// and Point4 are constrained relative to Line1. As a result, manipulating P3
// results in manipulating P4.
void vtkEllipseRepresentation::Point3WidgetInteraction(double e[2])
{
  double p1[3], p2[3], p3[3], p4[3];
  double slope1[3], slope2[3];

  // Start by getting the coordinates (P1,P2) defining Line1. Also get
  // characterisitics of Line1 including its slope, etc.
  this->GetPoint1WorldPosition(p1);
  this->GetPoint2WorldPosition(p2);
  slope1[0] = p2[0] - p1[0];
  slope1[1] = p2[1] - p1[1];
  slope2[0] = -slope1[1];
  slope2[1] =  slope1[0];
  slope2[2] = 0.0;
  vtkMath::Normalize(slope2);

  // The current position of P3 is constrained to lie along Line1. Also,
  // P4 is placed on the opposite side of Line1.
  double pw[4], t, closest[3];
  if ( this->Renderer )
    {
    this->Renderer->SetDisplayPoint(e[0],e[1],0.0);
    this->Renderer->DisplayToWorld();
    this->Renderer->GetWorldPoint(pw);
    }
  double dist = sqrt(vtkLine::DistanceToLine(pw,p1,p2,t,closest));

  // Set the positions of P3 and P4.
  p3[0] = closest[0] + dist*slope2[0];
  p3[1] = closest[1] + dist*slope2[1];
  p3[2] = pw[2];
  this->SetPoint3WorldPosition(p3);

  p4[0] = closest[0] - dist*slope2[0];
  p4[1] = closest[1] - dist*slope2[1];
  p4[2] = pw[2];
  this->SetPoint4WorldPosition(p4);
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::StartWidgetManipulation(double e[2])
{
  this->StartEventPosition[0] = e[0];
  this->StartEventPosition[1] = e[1];
  this->StartEventPosition[2] = 0.0;

  if ( this->Renderer )
    {
    this->Renderer->SetDisplayPoint(e[0],e[1],0.0);
    this->Renderer->DisplayToWorld();
    this->Renderer->GetWorldPoint(this->StartEventPositionWorld);
    }

  this->GetPoint1WorldPosition(this->P1World);
  this->GetPoint2WorldPosition(this->P2World);
  this->GetPoint3WorldPosition(this->P3World);
  this->GetPoint4WorldPosition(this->P4World);

  int i;
  for (i=0; i<3; i++)
    {
    this->P21World[i] = this->P2World[i] - this->P1World[i];
    this->P43World[i] = this->P4World[i] - this->P3World[i];
    }

  vtkLine::Intersection(this->P1World,this->P2World,
                        this->P3World,this->P4World,
                        this->T21,this->T43);

  // Compute the center point
  for (i=0; i<3; i++)
    {
    this->CenterWorld[i] = ((this->P1World[i] + this->T21*this->P21World[i]) +
                            (this->P3World[i] + this->T43*this->P43World[i]))/2.0;
    }
}

//----------------------------------------------------------------------
// This handles all the nasty special cases when the length of the arms of the
// bidimensional widget become zero. Basically the method prevents the arms
// from getting too short.

// cpbotha: a better description of what exactly this method does is
// desirable.  it's more than just a check on the arm length...

// my modified version does exactly what the original BiDirectional
// version did, but it also moves the opposite handle in the opposite
// direction, so that the axes always remain symmetric
void vtkEllipseRepresentation::ProjectOrthogonalPoint(
  double x[4], double y[3], double x1[3], double x2[3],
  double x21[3], double dir, double xP[3], double yP[3])
{
  double t, closest[3], slope[3], dist;

  // determine the distance from the other (orthogonal) line
  // DistanceToLine computes the squared distance of point x to the
  // line (x1,x2).  closest is the point on the line itself, t is
  // its parametric coordinate
  dist = dir * sqrt(vtkLine::DistanceToLine(x,x1,x2,t,closest));

  // get the closest point on the other line, use its "mate" point to
  // define the projection point, this keeps everything orthogonal.
  vtkLine::DistanceToLine(y,x1,x2,t,closest);

  // Project the point "dist" orthogonal to ray x21.
  // Define an orthogonal line.
  slope[0] = -x21[1];
  slope[1] =  x21[0];
  slope[2] = 0.0;

  // Project out the right distance along the calculated slope
  vtkMath::Normalize(slope);
  xP[0] = closest[0] + dist*slope[0];
  xP[1] = closest[1] + dist*slope[1];
  xP[2] = closest[2] + dist*slope[2];

  // Check to see what side the projection is on, clamp if necessary.
  // Note that closest is modified so that the arms don't end up with
  // zero length.
  if ( ((xP[0]-closest[0])*(x[0]-closest[0]) + (xP[1]-closest[1])*(x[1]-closest[1]) + (xP[2]-closest[2])*(x[2]-closest[2])) < 0.0 )
    {
    // Convert closest point to display coordinates
    double c1[3], c2[3], c21[3], cNew[3], xPNew[4];
    this->Renderer->SetWorldPoint(closest[0],closest[1],closest[2],1.0);
    this->Renderer->WorldToDisplay();
    this->Renderer->GetDisplayPoint(c1);
    // Convert vector in world space to display space
    this->Renderer->SetWorldPoint(closest[0]+dir*slope[0],closest[1]+dir*slope[1],closest[2]+dir*slope[2],1.0);
    this->Renderer->WorldToDisplay();
    this->Renderer->GetDisplayPoint(c2);
    c21[0] = c2[0] - c1[0];
    c21[1] = c2[1] - c1[1];
    c21[2] = c2[2] - c1[2];
    vtkMath::Normalize(c21);

    // Perform vector addition in display space to get new point
    cNew[0] = c1[0] + c21[0];
    cNew[1] = c1[1] + c21[1];
    cNew[2] = c1[2] + c21[2];

    this->Renderer->SetDisplayPoint(cNew[0],cNew[1],cNew[2]);
    this->Renderer->DisplayToWorld();
    this->Renderer->GetWorldPoint(xPNew);

    xP[0] = xPNew[0];
    xP[1] = xPNew[1];
    xP[2] = xPNew[2];
    }

  // cpbotha: just mirror projection on one side to the other side so
  // that ellipse will remain symmetric.
  vtkLine::DistanceToLine(xP,x1,x2,t,closest);
  for (int i = 0; i < 3; i++)
    {
    yP[i] = closest[i] - (xP[i] - closest[i]);
    }


}

//----------------------------------------------------------------------
// This method is tricky because it is constrained by Line1 and Line2.
// This method is invoked after all four points have been placed.
void vtkEllipseRepresentation::WidgetInteraction(double e[2])
{
  // Depending on the state, different motions are allowed.
  if ( this->InteractionState == Outside || ! this->Renderer )
    {
    return;
    }

  // Okay, go to work, convert this event to world coordinates
  double pw[4]; //t, closest[3];
  double p1[3], p2[3], p3[3], p4[3];
  this->Renderer->SetDisplayPoint(e[0],e[1],0.0);
  this->Renderer->DisplayToWorld();
  this->Renderer->GetWorldPoint(pw);

  // depending on the state, perform different operations
  if ( this->InteractionState == OnCenter )
    {
    for (int i=0; i<3; i++)
      {
      p1[i] = this->P1World[i] + (pw[i]-this->StartEventPositionWorld[i]);
      p2[i] = this->P2World[i] + (pw[i]-this->StartEventPositionWorld[i]);
      p3[i] = this->P3World[i] + (pw[i]-this->StartEventPositionWorld[i]);
      p4[i] = this->P4World[i] + (pw[i]-this->StartEventPositionWorld[i]);
      }
    this->SetPoint1WorldPosition(p1);
    this->SetPoint2WorldPosition(p2);
    this->SetPoint3WorldPosition(p3);
    this->SetPoint4WorldPosition(p4);
    }
  else if ( this->InteractionState == OnL1 ||
            this->InteractionState == OnL2) //rotate the representation
    {
    // compute rotation angle and center of rotation
    double sc[3], ec[3], p1c[3], p2c[3], p3c[3], p4c[3];
    for (int i=0; i<3; i++)
      {
      sc[i] = this->StartEventPositionWorld[i] - this->CenterWorld[i];
      ec[i] = pw[i] - this->CenterWorld[i];
      p1c[i] = this->P1World[i] - this->CenterWorld[i];
      p2c[i] = this->P2World[i] - this->CenterWorld[i];
      p3c[i] = this->P3World[i] - this->CenterWorld[i];
      p4c[i] = this->P4World[i] - this->CenterWorld[i];
      }
    double theta = atan2(ec[1],ec[0]) - atan2(sc[1],sc[0]);
    double r1 = vtkMath::Norm(p1c);
    double r2 = vtkMath::Norm(p2c);
    double r3 = vtkMath::Norm(p3c);
    double r4 = vtkMath::Norm(p4c);
    double theta1 = atan2(p1c[1],p1c[0]);
    double theta2 = atan2(p2c[1],p2c[0]);
    double theta3 = atan2(p3c[1],p3c[0]);
    double theta4 = atan2(p4c[1],p4c[0]);

    //rotate the four points
    p1[0] = this->CenterWorld[0] + r1*cos(theta+theta1);
    p1[1] = this->CenterWorld[1] + r1*sin(theta+theta1);
    p2[0] = this->CenterWorld[0] + r2*cos(theta+theta2);
    p2[1] = this->CenterWorld[1] + r2*sin(theta+theta2);
    p3[0] = this->CenterWorld[0] + r3*cos(theta+theta3);
    p3[1] = this->CenterWorld[1] + r3*sin(theta+theta3);
    p4[0] = this->CenterWorld[0] + r4*cos(theta+theta4);
    p4[1] = this->CenterWorld[1] + r4*sin(theta+theta4);
    p1[2] = this->P1World[2];
    p2[2] = this->P2World[2];
    p3[2] = this->P3World[2];
    p4[2] = this->P4World[2];

    this->SetPoint1WorldPosition(p1);
    this->SetPoint2WorldPosition(p2);
    this->SetPoint3WorldPosition(p3);
    this->SetPoint4WorldPosition(p4);
    }

  // pw is e (display coordinate of event) converted to world coords
  // ProjectOrthogonalPoint takes pw, does stuff, and spits out new p
  // we'll change this to spit out two points and not only one...
  else if ( this->InteractionState == NearP1 )
    {
    this->ProjectOrthogonalPoint(pw,this->P2World,this->P3World,this->P4World,this->P43World,-1,p1,p2);
    this->SetPoint1WorldPosition(p1);
    this->SetPoint2WorldPosition(p2);
    }
  else if ( this->InteractionState == NearP2 )
    {
    this->ProjectOrthogonalPoint(pw,this->P1World,this->P3World,this->P4World,this->P43World,1,p2,p1);
    this->SetPoint2WorldPosition(p2);
    this->SetPoint1WorldPosition(p1);
    }
  else if ( this->InteractionState == NearP3 )
     {
     this->ProjectOrthogonalPoint(pw,this->P4World,this->P1World,this->P2World,this->P21World,1,p3,p4);
    this->SetPoint3WorldPosition(p3);
    this->SetPoint4WorldPosition(p4);
    }
  else if ( this->InteractionState == NearP4 )
    {
    this->ProjectOrthogonalPoint(pw,this->P3World,this->P1World,this->P2World,this->P21World,-1,p4,p3);
    this->SetPoint4WorldPosition(p4);
    this->SetPoint3WorldPosition(p3);
    } //near P4
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetCenterWorldPosition(double pos[3])
{
	// we use DistanceToLine to calculate the crossing point of the two lines
	// in world coordinates (pos).  This should be the centre of the ellipse.
	double t;
	double p1[3], p2[3], p3[3];

	this->GetPoint1WorldPosition(p1);
    this->GetPoint2WorldPosition(p2);
	this->GetPoint3WorldPosition(p3);

	vtkLine::DistanceToLine(p3,p1,p2,t,pos);
}

double vtkEllipseRepresentation::GetSemiMajorAxisLength()
{
	double v[3];
	this->GetSemiMajorAxisVector(v);

	double len = 0.0;

	for (int i = 0; i < 3; i++)
	{
		len += v[i] * v[i];
	}

	return sqrt(len);

}

void vtkEllipseRepresentation::GetSemiMajorAxisVector(double v[3])
{
	double p1[3], p2[3];
	this->GetPoint1WorldPosition(p1);
    this->GetPoint2WorldPosition(p2);

	// by convention, p1-p2 is Major
	// semimajor is half of that
	for (int i = 0; i < 3; i++)
	{
		v[i] = (p2[i] - p1[i]) / 2.0;
	}
}

double vtkEllipseRepresentation::GetSemiMinorAxisLength()
{
	double v[3];
	this->GetSemiMinorAxisVector(v);

	double len = 0.0;

	for (int i = 0; i < 3; i++)
	{
		len += v[i] * v[i];
	}

	return sqrt(len);
}

void vtkEllipseRepresentation::GetSemiMinorAxisVector(double v[3])
{
	double p3[3], p4[3];
	this->GetPoint3WorldPosition(p3);
    this->GetPoint4WorldPosition(p4);

	// by convention, p3-p4 is Minor
	// semiminor is half of that
	for (int i = 0; i < 3; i++)
	{
		v[i] = (p4[i] - p3[i]) / 2.0;
	}
}



//----------------------------------------------------------------------
// this is called with DISPLAY coordinates
void vtkEllipseRepresentation::UpdateEllipse(double p1[3], double p2[3], double p3[3], double p4[3])
{
	// update ellipse (factor this out into method)
	// 1. set position (center); intersection of (wp1,wp2) and (wp3,wp4);
	//    also halfway point on each of these lines due to symmetry
	double t, closest[3];
	vtkLine::DistanceToLine(p1,p3,p4,t,closest); // closest is intersection point
	this->EllipseActor->SetPosition(closest);

	// 2. set semimajor and semiminor axes
	// for simplicity's sake, we assume that half of (p1,p2) is the semimajor axis
	double semi_major[3], semi_major_norm, semi_minor[3], semi_minor_norm;
	for (int i = 0; i < 3; i++)
	{
		semi_major[i] = p2[i] - p1[i];
		semi_minor[i] = p4[i] - p3[i];
	}

	semi_major_norm = vtkMath::Norm(semi_major);
	this->EllipseSource->SetSemiMajorAxisLength(semi_major_norm / 2.0);

	semi_minor_norm = vtkMath::Norm(semi_minor);
	this->EllipseSource->SetSemiMinorAxisLength(semi_minor_norm / 2.0);

	// 3. calculate and set angle.
	this->EllipseTrfm->Identity();
	this->EllipseTrfm->RotateZ(atan2(semi_major[1], semi_major[0]) / vtkMath::Pi() * 180);
	this->EllipseTrfmFilter->Update();
}

void vtkEllipseRepresentation::BuildRepresentation()
{
  if ( this->GetMTime() > this->BuildTime ||
       this->Point1Representation->GetMTime() > this->BuildTime ||
       this->Point2Representation->GetMTime() > this->BuildTime ||
       this->Point3Representation->GetMTime() > this->BuildTime ||
       this->Point4Representation->GetMTime() > this->BuildTime ||
       (this->Renderer && this->Renderer->GetVTKWindow() &&
        this->Renderer->GetVTKWindow()->GetMTime() > this->BuildTime) )
    {
    this->Point1Representation->BuildRepresentation();
    this->Point2Representation->BuildRepresentation();
    this->Point3Representation->BuildRepresentation();
    this->Point4Representation->BuildRepresentation();

    // Now bring the lines up to date
    if ( ! this->Line1Visibility )
      {
      return;
      }

    char distStr1[256], distStr2[256];
    double p1[3], p2[3], p3[3], p4[3];
    this->GetPoint1DisplayPosition(p1);
    this->GetPoint2DisplayPosition(p2);
    this->GetPoint3DisplayPosition(p3);
    this->GetPoint4DisplayPosition(p4);

    double wp1[3], wp2[3], wp3[3], wp4[3];
    this->GetPoint1WorldPosition(wp1);
    this->GetPoint2WorldPosition(wp2);
    this->GetPoint3WorldPosition(wp3);
    this->GetPoint4WorldPosition(wp4);

	// i'm not entirely sure why we're getting display coords here
	// (as is all the other geometry)
	this->UpdateEllipse(p1, p2, p3, p4);

    this->LinePoints->SetPoint(0,p1);
    this->LinePoints->SetPoint(1,p2);
    this->LinePoints->SetPoint(2,p3);
    this->LinePoints->SetPoint(3,p4);
    this->LinePoints->Modified();

    this->LineCells->Reset();
    this->LineCells->InsertNextCell(2);
    this->LineCells->InsertCellPoint(0);
    this->LineCells->InsertCellPoint(1);

    if ( this->Line2Visibility )
      {
      this->LineCells->InsertNextCell(2);
      this->LineCells->InsertCellPoint(2);
      this->LineCells->InsertCellPoint(3);
      }

    double line1Dist = sqrt(vtkMath::Distance2BetweenPoints(wp1, wp2));
    double line2Dist = 0;
    if (this->Line2Visibility)
      {
      line2Dist = sqrt(vtkMath::Distance2BetweenPoints(wp3, wp4));
      }
    ostrstream label;
    if (this->IDInitialized)
      {
      label << this->ID << ": ";
      }
    sprintf(distStr1,this->LabelFormat, line1Dist);
    sprintf(distStr2,this->LabelFormat, line2Dist);

    if (line1Dist > line2Dist)
      {
      label << distStr1 << " x " << distStr2 << ends;
      }
    else
      {
      label << distStr2 << " x " << distStr1 << ends;
      }
    this->TextMapper->SetInput(label.str());
    label.rdbuf()->freeze(0);

    // Adjust the font size
    int stringSize[2], *winSize = this->Renderer->GetSize();
    vtkTextMapper::SetRelativeFontSize(this->TextMapper, this->Renderer, winSize,
                                       stringSize, 0.015);

    int maxX = VTK_INT_MIN, maxY = VTK_INT_MIN;
    if (p1[1] > maxY)
      {
      maxX = (int)p1[0];
      maxY = (int)p1[1];
      }
    if (p2[1] > maxY)
      {
      maxX = (int)p2[0];
      maxY = (int)p2[1];
      }
    if (p3[1] > maxY)
      {
      maxX = (int)p3[0];
      maxY = (int)p3[1];
      }
    if (p4[1] > maxY)
      {
      maxX = (int)p4[0];
      maxY = (int)p4[1];
      }
    int minX = VTK_INT_MAX, minY = VTK_INT_MAX;
    if (p1[1] < minY)
      {
      minX = (int)p1[0];
      minY = (int)p1[1];
      }
    if (p2[1] < minY)
      {
      minX = (int)p2[0];
      minY = (int)p2[1];
      }
    if (p3[1] < minY)
      {
      minX = (int)p3[0];
      minY = (int)p3[1];
      }
    if (p4[1] < minY)
      {
      minX = (int)p4[0];
      minY = (int)p4[1];
      }
    int textSize[2];
    this->TextMapper->GetSize(this->Renderer, textSize);
    if (this->ShowLabelAboveWidget)
      {
      this->TextActor->SetPosition(maxX - textSize[0]/2, maxY+9);
      }
    else
      {
      this->TextActor->SetPosition(minX - textSize[0]/2, minY-(textSize[1]+9));
      }

    this->BuildTime.Modified();
    }
}

//----------------------------------------------------------------------
char* vtkEllipseRepresentation::GetLabelText()
{
  return this->TextMapper->GetInput();
}

//----------------------------------------------------------------------
double* vtkEllipseRepresentation::GetLabelPosition()
{
  return this->TextActor->GetPosition();
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::GetLabelPosition(double pos[3])
{
  this->TextActor->GetPositionCoordinate()->GetValue(pos);
}

//----------------------------------------------------------------------
double vtkEllipseRepresentation::GetLength1()
{
  double x1[3], x2[3];

  this->GetPoint1WorldPosition(x1);
  this->GetPoint2WorldPosition(x2);

  return sqrt(vtkMath::Distance2BetweenPoints(x1,x2));
}


//----------------------------------------------------------------------
double vtkEllipseRepresentation::GetLength2()
{
  double x3[3], x4[3];

  this->GetPoint3WorldPosition(x3);
  this->GetPoint4WorldPosition(x4);

  return sqrt(vtkMath::Distance2BetweenPoints(x3,x4));
}


//----------------------------------------------------------------------
void vtkEllipseRepresentation::ReleaseGraphicsResources(vtkWindow *w)
{
  this->LineActor->ReleaseGraphicsResources(w);
  this->TextActor->ReleaseGraphicsResources(w);
  this->EllipseActor->ReleaseGraphicsResources(w);
}


//----------------------------------------------------------------------
int vtkEllipseRepresentation::RenderOverlay(vtkViewport *viewport)
{
  this->BuildRepresentation();

  int count = this->LineActor->RenderOverlay(viewport);
  if ( this->Line1Visibility )
    {
    count += this->TextActor->RenderOverlay(viewport);
	count += this->EllipseActor->RenderOverlay(viewport);
    }

  return count;
}


//----------------------------------------------------------------------
void vtkEllipseRepresentation::Highlight(int highlightOn)
{
  if ( highlightOn )
    {
    this->LineActor->SetProperty(this->SelectedLineProperty);
    }
  else
    {
    this->LineActor->SetProperty(this->LineProperty);
    }
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::SetID(unsigned long id)
{
  this->ID = id;
  this->IDInitialized = 1;
  this->Modified();
}

//----------------------------------------------------------------------
void vtkEllipseRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  //Superclass typedef defined in vtkTypeMacro() found in vtkSetGet.h
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Tolerance: " << this->Tolerance << "\n";

  os << indent << "Length1: " << this->GetLength1() << "\n";
  os << indent << "Length2: " << this->GetLength2() << "\n";

  os << indent << "Line1 Visibility: " << (this->Line1Visibility ? "On\n" : "Off\n");
  os << indent << "Line2 Visibility: " << (this->Line2Visibility ? "On\n" : "Off\n");

  if ( this->TextProperty )
    {
    os << indent << "Text Property:\n";
    this->TextProperty->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Property: (none)\n";
    }

  if ( this->LineProperty )
    {
    os << indent << "Line Property:\n";
    this->LineProperty->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Line Property: (none)\n";
    }

  if ( this->SelectedLineProperty )
    {
    os << indent << "Selected Line Property:\n";
    this->SelectedLineProperty->PrintSelf(os,indent.GetNextIndent());
    }
  else
    {
    os << indent << "Selected Line Property: (none)\n";
    }

  os << indent << "Handle Representation: " << this->HandleRepresentation << "\n";
}
