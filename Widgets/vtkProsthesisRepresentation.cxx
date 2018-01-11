#include "vtkProsthesisRepresentation.h"

#include <cmath>

#include "vtkActor.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataMapper.h"
#include "vtkArrowSource.h"
#include "vtkRegularPolygonSource.h"
#include "vtkAppendPolyData.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkInteractorObserver.h"
#include "vtkCellPicker.h"
#include "vtkCamera.h"
#include "vtkAssemblyPath.h"
#include "vtkWindow.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkTransform.h"
#include "vtkMath.h"

#define PI 3.14159265358979323846

vtkStandardNewMacro(vtkProsthesisRepresentation);

//----------------------------------------------------------------------------
vtkProsthesisRepresentation::vtkProsthesisRepresentation() :
  vtkWidgetRepresentation(),
  ShowOutline(true),
  LeftArrowPolyData(vtkPolyData::New()),
  RightArrowPolyData(vtkPolyData::New()),
  LeftArrowTransformer(vtkTransformPolyDataFilter::New()),
  RightArrowTransformer(vtkTransformPolyDataFilter::New()),
  Transform(vtkTransform::New())
{
  // The initial state
  this->InteractionState = vtkProsthesisRepresentation::Outside;

  // Handle size is in pixels for this widget
  this->HandleSize = 7.5;

  // Set up the initial properties
  this->CreateDefaultProperties();

  this->HandleGeometry = vtkRegularPolygonSource::New();
  this->HandleGeometry->SetNumberOfSides(18);
  this->HandleGeometry->SetRadius(this->HandleSize);
  this->HandleGeometry->SetCenter(0, 0, 0);
  this->HandleGeometry->GeneratePolylineOff();
  this->HandleMapper = vtkPolyDataMapper::New();
  this->HandleMapper->SetInputConnection(this->HandleGeometry->GetOutputPort());
  this->Handle = vtkActor::New();
  this->Handle->SetMapper(this->HandleMapper);
  this->Handle->SetProperty(this->HandleProperty);

  this->RotateHandleGeometry = vtkRegularPolygonSource::New();
  this->RotateHandleGeometry->SetNumberOfSides(18);
  this->RotateHandleGeometry->SetRadius(this->HandleSize);
  this->RotateHandleGeometry->SetCenter(0, 0, 0);
  this->RotateHandleGeometry->GeneratePolylineOff();

  this->ArrowTransform = vtkTransform::New();
  this->ArrowTransform->Identity();

  this->LeftArrowTransformer->SetInputData(this->LeftArrowPolyData);
  this->LeftArrowTransformer->SetTransform(this->ArrowTransform);
  this->RightArrowTransformer->SetInputData(this->RightArrowPolyData);
  this->RightArrowTransformer->SetTransform(this->ArrowTransform);

  // Combine the marker and arrows
  this->RotateGeometryCombiner = vtkAppendPolyData::New();
  this->RotateGeometryCombiner->AddInputConnection(this->RotateHandleGeometry->GetOutputPort());
  this->RotateGeometryCombiner->AddInputConnection(this->LeftArrowTransformer->GetOutputPort());
  this->RotateGeometryCombiner->AddInputConnection(this->RightArrowTransformer->GetOutputPort());

  this->RotateHandleMapper = vtkPolyDataMapper::New();
  this->RotateHandleMapper->SetInputConnection(this->RotateGeometryCombiner->GetOutputPort());
  this->RotateHandle = vtkActor::New();
  this->RotateHandle->SetMapper(this->RotateHandleMapper);
  this->RotateHandle->SetProperty(this->HandleProperty);

  // Construct initial points
  this->Points = vtkPoints::New(VTK_DOUBLE);
  this->Points->SetNumberOfPoints(8); // 8 corners of the bounds

  // Create the outline for the hex
  this->OutlinePolyData = vtkPolyData::New();
  this->OutlinePolyData->SetPoints(this->Points);
  this->OutlineMapper = vtkPolyDataMapper::New();
  this->OutlineMapper->SetInputData(this->OutlinePolyData);
  this->Outline = vtkActor::New();
  this->Outline->SetMapper(this->OutlineMapper);
  this->Outline->SetProperty(this->OutlineProperty);
  {
    vtkCellArray* cells = vtkCellArray::New();
    cells->Allocate(cells->EstimateSize(15,2));
    this->OutlinePolyData->SetLines(cells);
    cells->Delete();
  }

  // Create the outline
  this->GenerateOutline();
  this->GenerateArrow(this->RightArrowPolyData, 0.05, true);
  this->GenerateArrow(this->LeftArrowPolyData, 0.05, false);

  // Define the point coordinates
  double bounds[6] = {-1.0, 1.0, 
                      -1.0, 1.0,
                      -1.0, 1.0};
  this->PlaceWidget(bounds);

  //Manage the picking stuff
  this->HandlePicker = vtkCellPicker::New();
  this->HandlePicker->SetTolerance(0.001);
  this->HandlePicker->AddPickList(this->Handle);
  this->HandlePicker->AddPickList(this->RotateHandle);
  this->HandlePicker->PickFromListOn();
}

//----------------------------------------------------------------------------
vtkProsthesisRepresentation::~vtkProsthesisRepresentation()
{
  this->HandleGeometry->Delete();
  this->HandleMapper->Delete();
  this->Handle->Delete();

  this->RotateHandleGeometry->Delete();
  this->RotateHandleMapper->Delete();
  this->RotateHandle->Delete();

  this->HandlePicker->Delete();
  this->HandleProperty->Delete();
  this->SelectedHandleProperty->Delete();

  this->OutlinePolyData->Delete();
  this->OutlineMapper->Delete();
  this->Outline->Delete();

  this->Points->Delete();
}

//----------------------------------------------------------------------
void vtkProsthesisRepresentation::StartWidgetInteraction(double e[2])
{
  // Store the start position
  this->StartEventPosition[0] = e[0];
  this->StartEventPosition[1] = e[1];
  this->StartEventPosition[2] = 0.0;

  // Store the start position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;

  this->ComputeInteractionState(static_cast<int>(e[0]), static_cast<int>(e[1]), 0);
}

//----------------------------------------------------------------------
void vtkProsthesisRepresentation::WidgetInteraction(double e[2])
{
  // Convert events to appropriate coordinate systems
  vtkCamera* camera = this->Renderer->GetActiveCamera();
  if (!camera)
  {
    return;
  }
  double focalPoint[4], pickPoint[4], prevPickPoint[4];
  double z, vpn[3];
  camera->GetViewPlaneNormal(vpn);

  // Compute the two points defining the motion vector
  double pos[3];
  // NOTE: This position is always the original pick position
  this->HandlePicker->GetPickPosition(pos);
  vtkInteractorObserver::ComputeWorldToDisplay(this->Renderer,
                                               pos[0], pos[1], pos[2],
                                               focalPoint);
  z = focalPoint[2];
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, this->LastEventPosition[0],
                                               this->LastEventPosition[1], z, prevPickPoint);
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, e[0], e[1], z, pickPoint);

  if (this->InteractionState == vtkProsthesisRepresentation::Translating)
  {
    this->Translate(prevPickPoint, pickPoint);
  }
  else if (this->InteractionState == vtkProsthesisRepresentation::Rotating)
  {
    this->Rotate(this->LastEventPosition[0], this->LastEventPosition[1],
                 e[0], e[1], vpn);
  }

  // Store the start position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::Translate(double *p1, double *p2)
{
  double v[3];

  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  v[2] = p2[2] - p1[2];

  // Move the center point
  this->Center[0] += v[0];
  this->Center[1] += v[1];
  this->Center[2] += v[2];

  // Move the corners
  double* pts =
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  for (int i = 0; i < this->Points->GetNumberOfPoints(); i++) {
    *pts++ += v[0];
    *pts++ += v[1];
    *pts++ += v[2];
  }

  this->PositionHandles();
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::Rotate(double previousX, double previousY,
                                     double X, double Y,
                                     double *vpn)
{
  double displayCenter[4];
  vtkInteractorObserver::ComputeWorldToDisplay(this->Renderer,
                                               this->Center[0],
                                               this->Center[1],
                                               this->Center[2],
                                               displayCenter);

  double angle = atan2(Y - displayCenter[1], X - displayCenter[0]);
  double previousAngle = atan2(previousY - displayCenter[1], previousX - displayCenter[0]);

  // The widget will have to be rotated by the change in the angle since the previous event.
  double difference = (angle - previousAngle) * 180.0 / PI;

  //Manipulate the transform to reflect the rotation
  vtkTransform* transform = vtkTransform::New();
  transform->Identity();
  transform->Translate(this->Center[0], this->Center[1], this->Center[2]);
  transform->RotateWXYZ(difference, vpn);
  transform->Translate(-this->Center[0], -this->Center[1], -this->Center[2]);

  //Set the corners
  vtkPoints* newPts = vtkPoints::New(VTK_DOUBLE);
  transform->TransformPoints(this->Points, newPts);

  double* pts = 
    static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  for (int i = 0; i < this->Points->GetNumberOfPoints(); i++, pts += 3)
  {
    this->Points->SetPoint(i, newPts->GetPoint(i));
  }

  newPts->Delete();
  this->PositionHandles();
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::CreateDefaultProperties()
{
  // Handle properties
  this->HandleProperty = vtkProperty::New();
  this->HandleProperty->SetColor(1, 1, 1);

  this->SelectedHandleProperty = vtkProperty::New();
  this->SelectedHandleProperty->SetColor(1, 0, 0);

  // Outline properties
  this->OutlineProperty = vtkProperty::New();
  this->OutlineProperty->SetRepresentationToWireframe();
  this->OutlineProperty->SetAmbient(1.0);
  this->OutlineProperty->SetAmbientColor(1.0, 1.0, 1.0);
  this->OutlineProperty->SetLineWidth(1.0);
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::PlaceWidget(double bds[6])
{
  int i;
  double bounds[6], center[3];

  this->AdjustBounds(bds, bounds, center);

  this->Points->SetPoint(0, bounds[0], bounds[2], bounds[4]);
  this->Points->SetPoint(1, bounds[1], bounds[2], bounds[4]);
  this->Points->SetPoint(2, bounds[1], bounds[3], bounds[4]);
  this->Points->SetPoint(3, bounds[0], bounds[3], bounds[4]);
  this->Points->SetPoint(4, bounds[0], bounds[2], bounds[5]);
  this->Points->SetPoint(5, bounds[1], bounds[2], bounds[5]);
  this->Points->SetPoint(6, bounds[1], bounds[3], bounds[5]);
  this->Points->SetPoint(7, bounds[0], bounds[3], bounds[5]);

  for (i=0; i<6; i++)
  {
    this->InitialBounds[i] = bounds[i];
  }
  this->InitialLength = sqrt((bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
                             (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
                             (bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));

  this->PositionHandles();
  this->ValidPick = 1; //since we have set up widget
  this->SizeHandles();
}

//----------------------------------------------------------------------------
int vtkProsthesisRepresentation::ComputeInteractionState(int X, int Y, int modify)
{
  // Okay, we can process this. Try to pick handles first;
  // if no handles picked, then pick the bounding box.
  if (!this->Renderer || !this->Renderer->IsInViewport(X, Y))
  {
    this->InteractionState = vtkProsthesisRepresentation::Outside;
    return this->InteractionState;
  }

  vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 0., this->HandlePicker);

  if (path != NULL)
  {
    this->ValidPick = 1;
    vtkActor* pickedHandle = reinterpret_cast<vtkActor*>(path->GetFirstNode()->GetViewProp());
    if (pickedHandle == this->Handle)
    {
      this->InteractionState = vtkProsthesisRepresentation::Translating;
    }
    else if (pickedHandle == this->RotateHandle) {
      this->InteractionState = vtkProsthesisRepresentation::Rotating;
    }
  }
  else {
    this->InteractionState = vtkProsthesisRepresentation::Outside;
  }

  return this->InteractionState;
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::BuildRepresentation()
{
  // Rebuild only if necessary
  if ( this->GetMTime() > this->BuildTime ||
       (this->Renderer && this->Renderer->GetVTKWindow() &&
        (this->Renderer->GetVTKWindow()->GetMTime() > this->BuildTime ||
        this->Renderer->GetActiveCamera()->GetMTime() > this->BuildTime)) )
  {
    this->SizeHandles();
    this->BuildTime.Modified();
  }
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::ReleaseGraphicsResources(vtkWindow *w)
{
  // Release the graphics resources associated with the actors of the widget
  this->Handle->ReleaseGraphicsResources(w);
  this->RotateHandle->ReleaseGraphicsResources(w);
  this->Outline->ReleaseGraphicsResources(w);
  // this->Arrow->ReleaseGraphicsResources(w);
}

//----------------------------------------------------------------------------
int vtkProsthesisRepresentation::RenderOpaqueGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();

  // Render all the actors of the widget
  count += this->Outline->RenderOpaqueGeometry(v);
  // count += this->Arrow->RenderOpaqueGeometry(v);
  if (this->Handle->GetVisibility())
  {
    count += this->Handle->RenderOpaqueGeometry(v);
  }
  if (this->RotateHandle->GetVisibility())
  {
    count += this->RotateHandle->RenderOpaqueGeometry(v);
  }

  return count;
}

//----------------------------------------------------------------------------
int vtkProsthesisRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();

  // Render all transparent actors of the widget
  count += this->Outline->RenderTranslucentPolygonalGeometry(v);
  // count += this->Arrow->RenderTranslucentPolygonalGeometry(v);
  if (this->Handle->GetVisibility())
  {
    count += this->Handle->RenderTranslucentPolygonalGeometry(v);
  }
  if (this->RotateHandle->GetVisibility())
  {
    count += this->RotateHandle->RenderTranslucentPolygonalGeometry(v);
  }

  return count;
}

//----------------------------------------------------------------------------
int vtkProsthesisRepresentation::HasTranslucentPolygonalGeometry()
{
  int result=0;
  this->BuildRepresentation();

  result |= this->Outline->HasTranslucentPolygonalGeometry();
  // result |= this->Arrow->HasTranslucentPolygonalGeometry();
  result |= this->Handle->HasTranslucentPolygonalGeometry();
  result |= this->RotateHandle->HasTranslucentPolygonalGeometry();

  return result;
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::PositionHandles()
{
  this->HandleGeometry->SetCenter(this->Center);
  double start[3];
  double up[3];
  this->Points->GetPoint(0, start);
  this->Points->GetPoint(1, up);
  double handlePos[3];
  handlePos[0] = this->Center[0] + (up[0] - start[0]) / 2.0;
  handlePos[1] = this->Center[1] + (up[1] - start[1]) / 2.0;
  handlePos[2] = this->Center[2] + (up[2] - start[2]) / 2.0;
  this->RotateHandleGeometry->SetCenter(handlePos);
  this->GenerateOutline();
  // Required so the handles stay the right size on screen during interaction.
  this->SizeHandles();
  this->UpdateTransform();
  this->GetTransform(this->ArrowTransform);
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::SizeHandles()
{
  this->HandleGeometry->SetRadius(
    this->vtkWidgetRepresentation::SizeHandlesInPixels(1.0,
                                                       this->Center));
  double rotateRadius = this->vtkWidgetRepresentation::SizeHandlesInPixels(1.0, this->RotateHandle->GetCenter());
  this->RotateHandleGeometry->SetRadius(rotateRadius);
  this->GenerateArrow(this->RightArrowPolyData, rotateRadius, true);
  this->GenerateArrow(this->LeftArrowPolyData, rotateRadius, false);
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::HighlightHandle(vtkProp* prop)
{
  // Unhighlight all
  this->Handle->SetProperty(this->HandleProperty);
  this->RotateHandle->SetProperty(this->HandleProperty);

  vtkActor* actorToHighlight = static_cast<vtkActor*>(prop);
  if (actorToHighlight)
  {
    actorToHighlight->SetProperty(this->SelectedHandleProperty);
  }
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::SetInteractionState(int state)
{
  // Clamp to allowable values
  state = (state < vtkProsthesisRepresentation::Outside ? vtkProsthesisRepresentation::Outside :
          (state > vtkProsthesisRepresentation::Scaling ? vtkProsthesisRepresentation::Scaling : state));

  this->InteractionState = state;
  switch (state)
  {
    case vtkProsthesisRepresentation::Translating:
      this->HighlightHandle(this->Handle);
      break;
    case vtkProsthesisRepresentation::Rotating:
      this->HighlightHandle(this->RotateHandle);
      break;
    case vtkProsthesisRepresentation::Scaling:
    default:
      this->HighlightHandle(NULL);
  }
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::SetShowOutline(bool show) {
  if (this->ShowOutline != show) {
    this->ShowOutline = show;
    this->GenerateOutline();
  }
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::GenerateOutline()
{
  // Whatever the case may be, we have to reset the Lines of the
  // OutlinePolyData (i.e. nuke all current line data)
  vtkCellArray* cells = this->OutlinePolyData->GetLines();
  cells->Reset();

  if (this->ShowOutline)
  {
    vtkIdType pts[2];
    pts[0] = 0; pts[1] = 1;
    cells->InsertNextCell(2, pts);
    pts[0] = 1; pts[1] = 2;
    cells->InsertNextCell(2, pts);
    pts[0] = 2; pts[1] = 3;
    cells->InsertNextCell(2, pts);
    pts[0] = 3; pts[1] = 0;
    cells->InsertNextCell(2, pts);
    pts[0] = 4; pts[1] = 5;
    cells->InsertNextCell(2, pts);
    pts[0] = 5; pts[1] = 6;
    cells->InsertNextCell(2, pts);
    pts[0] = 6; pts[1] = 7;
    cells->InsertNextCell(2, pts);
    pts[0] = 7; pts[1] = 4;
    cells->InsertNextCell(2, pts);
    pts[0] = 0; pts[1] = 4;
    cells->InsertNextCell(2, pts);
    pts[0] = 1; pts[1] = 5;
    cells->InsertNextCell(2, pts);
    pts[0] = 2; pts[1] = 6;
    cells->InsertNextCell(2, pts);
    pts[0] = 3; pts[1] = 7;
    cells->InsertNextCell(2, pts);
    this->OutlinePolyData->Modified();
  }
  if (this->OutlineProperty)
  {
    this->OutlineProperty->SetRepresentationToWireframe();
  }
}

void vtkProsthesisRepresentation::GenerateArrow(vtkPolyData* arrowPolyData,
                                                double shaftWidth, 
                                                bool clockwise, bool outline)
{
  // Parameters needed to generate an arrow
  // The radius the arrow is drawn from the center.
  double radius = 0.5;
  // The thickness of the tip.
  double tipSize = 3 * shaftWidth;
  // Angle where the arrow starts
  double startAngle = PI / 2.0;
  // The length of the arrow.
  double arrowLength = 10 * shaftWidth;
  if (clockwise) arrowLength *= -1;
  // The center of the circle the arrow is drawn on.
  double center[3];
  center[0] = center[1] = center[2] = 0.0;
  // Number of segments used to build the shaft section.
  int numSegment = 20;

  // Usefull variables deduced from the parameters.

  // Angle where the arrow stops. If this value is smaller than the start value
  // the arrow will point clockwise, else counter-clockwise.
  // Note: To calculate a angle for a distance on the circumference of a circle
  // the distance must be divided by the circumference and the mutiplied by 2 
  // radians to get the angle; the equation can be simplified to distance 
  // divided by radius, though.
  double endAngle = startAngle + arrowLength / radius;
  // The angle between the start and end of the tip.
  double tipAngle = tipSize / radius;
  // The points include the outside curve's (segments + 1),
  // inside curve's (segments + 1) and 3 points for the tip.
  int numPoints = (numSegment + 1) * 2 + 3;
  // The radius of the outside of the shaft.
  double outsideRadius = radius + shaftWidth / 2.0;
  // The radius of the inside of the shaft.
  double insideRadius = radius - shaftWidth / 2.0;
  // The exact angle the tip starts at.
  double tipStartAngle = endAngle - startAngle;
  // Handle the direction of the arrow correctly...
  tipStartAngle = tipStartAngle < 0 ? tipStartAngle + tipAngle : tipStartAngle - tipAngle;
  // ... and add the start angle to the value.
  tipStartAngle += startAngle;

  vtkPoints* arrowPoints = vtkPoints::New();
  arrowPoints->SetNumberOfPoints(numPoints);
  // Position the shaft points
  double angle = startAngle;
  double angleInc = (tipStartAngle - startAngle) / numSegment;
  for (double i = 0; i <= numSegment; i++, angle += angleInc)
  {
    arrowPoints->SetPoint(i, outsideRadius * sin(angle), 0, outsideRadius * cos(angle));
    arrowPoints->SetPoint(numPoints - i - 1, insideRadius * sin(angle), 0, insideRadius * cos(angle));
  }

  // Position the tip points
  double tipOutsideRadius = radius + tipSize / 2.0;
  double tipInsideRadius = radius - tipSize / 2.0;
  arrowPoints->SetPoint(numSegment + 1, tipOutsideRadius * sin(tipStartAngle), 0, tipOutsideRadius * cos(tipStartAngle));
  arrowPoints->SetPoint(numSegment + 2, radius * sin(endAngle), 0, radius * cos(endAngle));
  arrowPoints->SetPoint(numSegment + 3, tipInsideRadius * sin(tipStartAngle), 0, tipInsideRadius * cos(tipStartAngle));
  arrowPolyData->SetPoints(arrowPoints);
  arrowPoints->Delete();

  // Create the lines of the arrow  
  vtkCellArray* cells = vtkCellArray::New();
  cells->Reset();
  int maxPointId = numPoints - 1;

  if (outline) {
    vtkIdType pts[2];
    for (double i = 0; i < maxPointId; i++)
    {
      pts[0] = i; pts[1] = i+1;
      cells->InsertNextCell(2, pts);
    }
    pts[0] = maxPointId; pts[1] = 0;
    cells->InsertNextCell(2, pts);
    arrowPolyData->SetLines(cells);
  }
  else {
    vtkIdType pts[3];
    for (double i = 0; i < numSegment; i++)
    {
      pts[0] = i; pts[1] = i+1; pts[2] = maxPointId - i - 1;
      cells->InsertNextCell(3, pts);
      pts[0] = maxPointId - i - 1; pts[1] = i+1; pts[2] = maxPointId - i - 2;
      cells->InsertNextCell(3, pts);
    }
    pts[0] = numSegment + 1; pts[1] = numSegment + 2; pts[2] = numSegment + 3;
    cells->InsertNextCell(3, pts);
    arrowPolyData->SetPolys(cells);
  }
  arrowPolyData->Modified();
}

void vtkProsthesisRepresentation::UpdateTransform()
{
  this->Transform->Identity();
  this->Transform->Translate(this->Center);

  // Orientation
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->Identity();
  
  // Compute normals
  double* pts =
     static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double* p0 = pts;
  double* px = pts + 3 * 1;
  double* py = pts + 3 * 3;
  double* pz = pts + 3 * 4;

  double nx[3], ny[3], nz[3]; // The normals vectors
  vtkMath::Subtract(p0, px, nx);
  vtkMath::Subtract(p0, py, ny);
  vtkMath::Subtract(p0, pz, nz);
  vtkMath::Normalize(nx);
  vtkMath::Normalize(ny);
  vtkMath::Normalize(nz);
  for (int i = 0; i < 3; i++)
  {
    matrix->SetElement(i, 0, -nx[i]);
    matrix->SetElement(i, 1, -ny[i]);
    matrix->SetElement(i, 2, -nz[i]);
  }
  this->Transform->Concatenate(matrix);
  matrix->Delete();
}

void vtkProsthesisRepresentation::GetTransform(vtkTransform* t)
{
  t->DeepCopy(this->Transform);
}

double* vtkProsthesisRepresentation::GetBounds()
{
  this->BuildRepresentation();
  return this->Outline->GetBounds();
}
