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
#include "vtkAlgorithmOutput.h"

#define PI 3.14159265358979323846

vtkStandardNewMacro(vtkProsthesisRepresentation);

//----------------------------------------------------------------------------
vtkProsthesisRepresentation::vtkProsthesisRepresentation() :
  vtkWidgetRepresentation(),
  Radius(0.5),
  ShowOutline(true),
  TranslateArrowPolyData(vtkPolyData::New()),
  PlusPolyData(vtkPolyData::New()),
  MinusPolyData(vtkPolyData::New()),
  LeftArrowPolyData(vtkPolyData::New()),
  RightArrowPolyData(vtkPolyData::New()),
  Transform(vtkTransform::New())
{
  this->Center[0] = this->Center[1] = this->Center[2] = 0.0;
  this->HandleColour[0] = this->HandleColour[1] = this->HandleColour[2] = 1.0;
  this->SelectedHandleColour[0] = 0.627;
  this->SelectedHandleColour[1] = 0.909;
  this->SelectedHandleColour[2] = 0.192;

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

  this->TranslateGeometryCombiner = vtkAppendPolyData::New();
  this->TranslateGeometryCombiner->AddInputConnection(this->HandleGeometry->GetOutputPort());
  this->TranslateGeometryCombiner->AddInputData(this->TranslateArrowPolyData);

  this->HandleMapper = vtkPolyDataMapper::New();
  this->TransformPolyData(this->TranslateGeometryCombiner->GetOutputPort(),
                          this->HandleMapper);

  this->Handle = vtkActor::New();
  this->Handle->SetMapper(this->HandleMapper);
  this->Handle->SetProperty(this->HandleProperty);

  this->UpHandleMapper = vtkPolyDataMapper::New();
  this->TransformPolyData(this->PlusPolyData,
                          this->UpHandleMapper);

  this->UpHandle = vtkActor::New();
  this->UpHandle->SetMapper(this->UpHandleMapper);
  this->UpHandle->SetProperty(this->HandleProperty);

  this->DownHandleMapper = vtkPolyDataMapper::New();
  this->TransformPolyData(this->MinusPolyData,
                          this->DownHandleMapper);

  this->DownHandle = vtkActor::New();
  this->DownHandle->SetMapper(this->DownHandleMapper);
  this->DownHandle->SetProperty(this->HandleProperty);

  this->RotateHandleGeometry = vtkRegularPolygonSource::New();
  this->RotateHandleGeometry->SetNumberOfSides(18);
  this->RotateHandleGeometry->SetRadius(this->HandleSize);
  this->RotateHandleGeometry->SetCenter(0, 0, 0);
  this->RotateHandleGeometry->GeneratePolylineOff();

  // Combine the marker and arrows
  this->RotateGeometryCombiner = vtkAppendPolyData::New();
  this->RotateGeometryCombiner->AddInputConnection(this->RotateHandleGeometry->GetOutputPort());
  this->RotateGeometryCombiner->AddInputData(this->LeftArrowPolyData);
  this->RotateGeometryCombiner->AddInputData(this->RightArrowPolyData);

  this->RotateHandleMapper = vtkPolyDataMapper::New();
  this->TransformPolyData(this->RotateGeometryCombiner->GetOutputPort(),
                          this->RotateHandleMapper);

  this->RotateHandle = vtkActor::New();
  this->RotateHandle->SetMapper(this->RotateHandleMapper);
  this->RotateHandle->SetProperty(this->HandleProperty);

  // Construct initial points
  this->Points = vtkPoints::New(VTK_DOUBLE);
  this->Points->SetNumberOfPoints(4); // center, x, y and z axis.

  // Create the outline circle
  this->OutlineGeometry = vtkRegularPolygonSource::New();
  this->OutlineGeometry->SetNumberOfSides(50);
  this->OutlineGeometry->SetRadius(this->Radius);
  this->OutlineGeometry->SetCenter(0, 0, 0);
  this->OutlineGeometry->GeneratePolygonOff();

  this->OutlineMapper = vtkPolyDataMapper::New();
  this->TransformPolyData(this->OutlineGeometry->GetOutputPort(),
                          this->OutlineMapper);

  this->Outline = vtkActor::New();
  this->Outline->SetMapper(this->OutlineMapper);
  this->Outline->SetProperty(this->OutlineProperty);

  this->GenerateTranslateArrow(0.05);
  this->GeneratePlus(1.0);
  this->GenerateMinus(1.0);
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
  this->HandlePicker->AddPickList(this->UpHandle);
  this->HandlePicker->AddPickList(this->DownHandle);
  this->HandlePicker->AddPickList(this->RotateHandle);
  this->HandlePicker->PickFromListOn();
}

void vtkProsthesisRepresentation::TransformPolyData(vtkAlgorithmOutput* algoOutput,
                                                    vtkPolyDataMapper* mapper)
{
  // NOTE: This is a temoporary solution to display the handles in the y-axis
  // direction, since that is how they will be used in HPS desktop.
  // TODO: Implement a feature that always arients the handles to face the
  // camera with a look-at transform matrix.
  vtkTransform* transform = vtkTransform::New();
  transform->Identity();
  transform->RotateX(90);
  vtkTransformPolyDataFilter* transformer = vtkTransformPolyDataFilter::New();
  transformer->SetInputConnection(algoOutput);
  transformer->SetTransform(transform);
  transform->Delete();

  mapper->SetInputConnection(transformer->GetOutputPort());
  transformer->Delete();
}

void vtkProsthesisRepresentation::TransformPolyData(vtkPolyData* polyData, vtkPolyDataMapper* mapper)
{
  // NOTE: This is a temoporary solution to display the handles in the y-axis
  // direction, since that is how they will be used in HPS desktop.
  // TODO: Implement a feature that always arients the handles to face the
  // camera with a look-at transform matrix.
  vtkTransform* transform = vtkTransform::New();
  transform->Identity();
  transform->RotateX(90);
  vtkTransformPolyDataFilter* transformer = vtkTransformPolyDataFilter::New();
  transformer->SetInputData(polyData);
  transformer->SetTransform(transform);
  transform->Delete();

  mapper->SetInputConnection(transformer->GetOutputPort());
  transformer->Delete();
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

  this->OutlineGeometry->Delete();
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
  double z;

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
                 e[0], e[1]);
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
                                         double X, double Y)
{

  // Compute normals
  double* pts =
     static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double* p0 = pts;
  double* py = pts + 3 * 2;
  double ny[3]; // The y-axis is the widget's face normal.
  vtkMath::Subtract(py, p0, ny);
  vtkMath::Normalize(ny);

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
  transform->RotateWXYZ(difference, ny);
  transform->Translate(-this->Center[0], -this->Center[1], -this->Center[2]);

  //Set the corners
  vtkPoints* newPts = vtkPoints::New(VTK_DOUBLE);
  transform->TransformPoints(this->Points, newPts);
  transform->Delete();

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
  this->HandleProperty->SetColor(this->HandleColour);

  this->SelectedHandleProperty = vtkProperty::New();
  this->SelectedHandleProperty->SetColor(this->SelectedHandleColour);

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

  this->Points->SetPoint(0, 0, 0, 0); // center
  this->Points->SetPoint(1, 1, 0, 0); // x
  this->Points->SetPoint(2, 0, 1, 0); // y
  this->Points->SetPoint(3, 0, 0, 1); // z

  for (i=0; i<6; i++)
  {
    this->InitialBounds[i] = bounds[i];
  }
  this->InitialLength = sqrt((bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
                             (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
                             (bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));

  this->PositionHandles();
  this->ValidPick = 1; //since we have set up widget
  this->UpdateHandles();
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
    else if (pickedHandle == this->UpHandle) {
      this->InteractionState = vtkProsthesisRepresentation::UpClick;
    }
    else if (pickedHandle == this->DownHandle) {
      this->InteractionState = vtkProsthesisRepresentation::DownClick;
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
    this->UpdateHandles();
    this->BuildTime.Modified();
  }
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::ReleaseGraphicsResources(vtkWindow *w)
{
  // Release the graphics resources associated with the actors of the widget
  this->Handle->ReleaseGraphicsResources(w);
  this->UpHandle->ReleaseGraphicsResources(w);
  this->DownHandle->ReleaseGraphicsResources(w);
  this->RotateHandle->ReleaseGraphicsResources(w);
  this->Outline->ReleaseGraphicsResources(w);
}

//----------------------------------------------------------------------------
int vtkProsthesisRepresentation::RenderOpaqueGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();

  // Render all the actors of the widget
  if (this->Outline->GetVisibility())
  {
    count += this->Outline->RenderOpaqueGeometry(v);
  }
  if (this->Handle->GetVisibility())
  {
    count += this->Handle->RenderOpaqueGeometry(v);
  }
  if (this->UpHandle->GetVisibility())
  {
    count += this->UpHandle->RenderOpaqueGeometry(v);
  }
  if (this->DownHandle->GetVisibility())
  {
    count += this->DownHandle->RenderOpaqueGeometry(v);
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
  if (this->Outline->GetVisibility())
  {
    count += this->Outline->RenderTranslucentPolygonalGeometry(v);
  }
  if (this->UpHandle->GetVisibility())
  {
    count += this->UpHandle->RenderTranslucentPolygonalGeometry(v);
  }
  if (this->DownHandle->GetVisibility())
  {
    count += this->DownHandle->RenderTranslucentPolygonalGeometry(v);
  }
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
  result |= this->UpHandle->HasTranslucentPolygonalGeometry();
  result |= this->DownHandle->HasTranslucentPolygonalGeometry();
  result |= this->Handle->HasTranslucentPolygonalGeometry();
  result |= this->RotateHandle->HasTranslucentPolygonalGeometry();

  return result;
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::PositionHandles()
{
  this->UpdateTransform();

  vtkTransform* t = vtkTransform::New();
  this->GetTransform(t);
  this->Outline->SetUserTransform(t);

  vtkTransform* tRot = vtkTransform::New();
  tRot->DeepCopy(t);
  tRot->Translate(this->Radius, 0, 0);
  this->RotateHandle->SetUserTransform(tRot);
  tRot->Delete();
  t->Delete();

  this->UpdateHandles();
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::UpdateHandles()
{
  double handleSizeCenter =
    this->vtkWidgetRepresentation::SizeHandlesInPixels(1.0, this->Center);

  // Transform the up and down handles so they're always right side up in 
  // relation to the camera.
  vtkTransform* tUp = vtkTransform::New();
  this->GetAlwaysUpTransform(tUp);

  this->Handle->SetUserTransform(tUp);

  vtkTransform* tPlus = vtkTransform::New();
  tPlus->DeepCopy(tUp);
  tPlus->Translate(-handleSizeCenter * 10, 0, handleSizeCenter * 3);
  this->UpHandle->SetUserTransform(tPlus);
  tPlus->Delete();

  vtkTransform* tMinus = vtkTransform::New();
  tMinus->DeepCopy(tUp);
  tMinus->Translate(-handleSizeCenter * 10, 0,-handleSizeCenter * 3);
  this->DownHandle->SetUserTransform(tMinus);
  tMinus->Delete();
  tUp->Delete();

  // Size the handles based on how far they are from the camera.
  this->HandleGeometry->SetRadius(handleSizeCenter);
  double rotateRadius = this->vtkWidgetRepresentation::SizeHandlesInPixels(1.0, this->RotateHandle->GetCenter());
  this->RotateHandleGeometry->SetRadius(rotateRadius);

  this->GeneratePlus(handleSizeCenter * 6);
  this->GenerateMinus(handleSizeCenter * 6);
  this->GenerateTranslateArrow(handleSizeCenter);
  this->GenerateArrow(this->RightArrowPolyData, rotateRadius, true);
  this->GenerateArrow(this->LeftArrowPolyData, rotateRadius, false);
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::HighlightHandle(vtkProp* prop)
{
  // Unhighlight all
  this->Handle->SetProperty(this->HandleProperty);
  this->RotateHandle->SetProperty(this->HandleProperty);
  this->UpHandle->SetProperty(this->HandleProperty);
  this->DownHandle->SetProperty(this->HandleProperty);

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
          (state > vtkProsthesisRepresentation::DownClick ? vtkProsthesisRepresentation::Outside : state));

  this->InteractionState = state;
  switch (state)
  {
    case vtkProsthesisRepresentation::Translating:
      this->HighlightHandle(this->Handle);
      break;
    case vtkProsthesisRepresentation::Rotating:
      this->HighlightHandle(this->RotateHandle);
      break;
    case vtkProsthesisRepresentation::UpClick:
      this->HighlightHandle(this->UpHandle);
      break;
    case vtkProsthesisRepresentation::DownClick:
      this->HighlightHandle(this->DownHandle);
      break;
    default:
      this->HighlightHandle(NULL);
  }
}

//----------------------------------------------------------------------------
void vtkProsthesisRepresentation::SetShowHandles(bool show)
{
  this->ShowHandles = show;
  this->Handle->SetVisibility(show);
  this->UpHandle->SetVisibility(show);
  this->DownHandle->SetVisibility(show);
  this->RotateHandle->SetVisibility(show);
}

void vtkProsthesisRepresentation::SetShowOutline(bool show) {
  this->ShowOutline = show;
  this->Outline->SetVisibility(show);
}

void vtkProsthesisRepresentation::GeneratePlus(double size, bool outline)
{
  double halfSize = size / 2.0;
  // Half the line thickness
  double lineSize = halfSize / 4.0; // line is a quarter of the size

  int numPoints = 12;
  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(numPoints);
  points->SetPoint( 0,  lineSize,  lineSize, 0);
  points->SetPoint( 1,  halfSize,  lineSize, 0);
  points->SetPoint( 2,  halfSize, -lineSize, 0);
  points->SetPoint( 3,  lineSize, -lineSize, 0);
  points->SetPoint( 4,  lineSize, -halfSize, 0);
  points->SetPoint( 5, -lineSize, -halfSize, 0);
  points->SetPoint( 6, -lineSize, -lineSize, 0);
  points->SetPoint( 7, -halfSize, -lineSize, 0);
  points->SetPoint( 8, -halfSize,  lineSize, 0);
  points->SetPoint( 9, -lineSize,  lineSize, 0);
  points->SetPoint(10, -lineSize,  halfSize, 0);
  points->SetPoint(11,  lineSize,  halfSize, 0);
  this->PlusPolyData->SetPoints(points);
  points->Delete();

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
    this->PlusPolyData->SetLines(cells);
  }
  else {
    vtkIdType pts[3];
    pts[0] = 1; pts[1] = 2; pts[2] = 7;
    cells->InsertNextCell(3, pts);
    pts[0] = 7; pts[1] = 8; pts[2] = 1;
    cells->InsertNextCell(3, pts);
    pts[0] = 4; pts[1] = 5; pts[2] = 10;
    cells->InsertNextCell(3, pts);
    pts[0] = 10; pts[1] = 11; pts[2] = 4;
    cells->InsertNextCell(3, pts);
    this->PlusPolyData->SetPolys(cells);
  }
  this->PlusPolyData->Modified();
}

void vtkProsthesisRepresentation::GenerateMinus(double size, bool outline)
{
  double halfSize = size / 2.0;
  // Half the line thickness
  double lineSize = halfSize / 4.0; // line is a quarter of the size

  int numPoints = 4;
  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(numPoints);
  points->SetPoint(0,  halfSize,  lineSize, 0);
  points->SetPoint(1,  halfSize, -lineSize, 0);
  points->SetPoint(2, -halfSize, -lineSize, 0);
  points->SetPoint(3, -halfSize,  lineSize, 0);
  this->MinusPolyData->SetPoints(points);
  points->Delete();

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
    this->MinusPolyData->SetLines(cells);
  }
  else {
    vtkIdType pts[3];
    pts[0] = 0; pts[1] = 1; pts[2] = 2;
    cells->InsertNextCell(3, pts);
    pts[0] = 2; pts[1] = 3; pts[2] = 0;
    cells->InsertNextCell(3, pts);
    this->MinusPolyData->SetPolys(cells);
  }
  this->MinusPolyData->Modified();
}

void vtkProsthesisRepresentation::GenerateTranslateArrow(double shaftWidth)
{
  double halfShaft = shaftWidth / 2.0;
  // The thickness of the tip.
  double tipSize = 3 * halfShaft;
  // The length of the shaft.
  double shaftLength = 7 * halfShaft;
  // The length of the arrow.
  double arrowLength = tipSize + shaftLength;
  // The center of the circle the arrow is drawn on.
  double center[3];
  center[0] = center[1] = center[2] = 0.0;

  int numPoints = 20; // 8 for the horizontal and vertical bars (shafts) and 
                      // 12 for the 4 triangles (tips).

  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(numPoints);
  points->SetPoint( 0, -shaftLength,  halfShaft, 0);
  points->SetPoint( 1,  shaftLength,  halfShaft, 0);
  points->SetPoint( 2,  shaftLength, -halfShaft, 0);
  points->SetPoint( 3, -shaftLength, -halfShaft, 0);
  points->SetPoint( 4, -halfShaft,  shaftLength, 0);
  points->SetPoint( 5,  halfShaft,  shaftLength, 0);
  points->SetPoint( 6,  halfShaft, -shaftLength, 0);
  points->SetPoint( 7, -halfShaft, -shaftLength, 0);
  points->SetPoint( 8,  shaftLength,  tipSize, 0);
  points->SetPoint( 9,  arrowLength,        0, 0);
  points->SetPoint(10,  shaftLength, -tipSize, 0);
  points->SetPoint(11, -shaftLength,  tipSize, 0);
  points->SetPoint(12, -arrowLength,        0, 0);
  points->SetPoint(13, -shaftLength, -tipSize, 0);
  points->SetPoint(14,  tipSize,  shaftLength, 0);
  points->SetPoint(15,        0,  arrowLength, 0);
  points->SetPoint(16, -tipSize,  shaftLength, 0);
  points->SetPoint(17,  tipSize, -shaftLength, 0);
  points->SetPoint(18,        0, -arrowLength, 0);
  points->SetPoint(19, -tipSize, -shaftLength, 0);
  this->TranslateArrowPolyData->SetPoints(points);
  points->Delete();

  // Create the lines of the arrow  
  vtkCellArray* cells = vtkCellArray::New();
  cells->Reset();
  vtkIdType pts[3];
  // horizontal shaft
  pts[0] = 0; pts[1] = 1; pts[2] = 2;
  cells->InsertNextCell(3, pts);
  pts[0] = 0; pts[1] = 2; pts[2] = 3;
  cells->InsertNextCell(3, pts);
  // vertical shaft
  pts[0] = 4; pts[1] = 5; pts[2] = 6;
  cells->InsertNextCell(3, pts);
  pts[0] = 4; pts[1] = 6; pts[2] = 7;
  // right arrow
  cells->InsertNextCell(3, pts);
  pts[0] = 8; pts[1] = 9; pts[2] = 10;
  cells->InsertNextCell(3, pts);
  // left arrow
  pts[0] = 11; pts[1] = 12; pts[2] = 13;
  cells->InsertNextCell(3, pts);
  // top arrow
  pts[0] = 14; pts[1] = 15; pts[2] = 16;
  cells->InsertNextCell(3, pts);
  // bottom arrow
  pts[0] = 17; pts[1] = 18; pts[2] = 19;
  cells->InsertNextCell(3, pts);
  this->TranslateArrowPolyData->SetPolys(cells);
  this->TranslateArrowPolyData->Modified();
}

void vtkProsthesisRepresentation::GenerateArrow(vtkPolyData* arrowPolyData,
                                                double shaftWidth, 
                                                bool clockwise, bool outline)
{
  // Parameters needed to generate an arrow
  // The radius the arrow is drawn from the center.
  double radius = this->Radius;
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
  double xOff = radius * sin(angle);
  double yOff = radius * cos(angle);

  for (double i = 0; i <= numSegment; i++, angle += angleInc)
  {
    arrowPoints->SetPoint(i, outsideRadius * sin(angle) - xOff, outsideRadius * cos(angle) - yOff, 0);
    arrowPoints->SetPoint(numPoints - i - 1, insideRadius * sin(angle) - xOff, insideRadius * cos(angle) - yOff, 0);
  }

  // Position the tip points
  double tipOutsideRadius = radius + tipSize / 2.0;
  double tipInsideRadius = radius - tipSize / 2.0;
  arrowPoints->SetPoint(numSegment + 1, tipOutsideRadius * sin(tipStartAngle) - xOff, tipOutsideRadius * cos(tipStartAngle) - yOff, 0);
  arrowPoints->SetPoint(numSegment + 2, radius * sin(endAngle) - xOff, radius * cos(endAngle) - yOff, 0);
  arrowPoints->SetPoint(numSegment + 3, tipInsideRadius * sin(tipStartAngle) - xOff, tipInsideRadius * cos(tipStartAngle) - yOff, 0);
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
  double* py = pts + 3 * 2;
  double* pz = pts + 3 * 3;

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

void vtkProsthesisRepresentation::GetAlwaysUpTransform(vtkTransform* t)
{
  t->Identity();
  t->Translate(this->Center);
  if (!this->Renderer) {
    return;
  }

  vtkCamera* camera = this->Renderer->GetActiveCamera();
  if (!camera)
  {
    return;
  }

  // Orientation
  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->Identity();

  // Compute normals
  double* pts =
     static_cast<vtkDoubleArray*>(this->Points->GetData())->GetPointer(0);
  double* p0 = pts;
  double* px = pts + 3 * 1;
  double* py = pts + 3 * 2;
  double* pz = pts + 3 * 3;

  double nx[3], ny[3], nz[3]; // The normals vectors
  vtkMath::Subtract(p0, px, nx);
  vtkMath::Subtract(p0, py, ny);
  vtkMath::Subtract(p0, pz, nz);
  vtkMath::Normalize(ny);

  vtkMath::Cross(camera->GetViewUp(), ny, nx);
  vtkMath::Normalize(nx);
  // TODO: Make this robust if the widget is changed to handle any up vector.
  if (nx[0] == 0 && nx[1] == 0 && nx[2] == 0) {
    std::cout << "Error: Up same as y-axis." << std::endl;
    return;
  }

  vtkMath::Cross(ny, nx, nz);

  for (int i = 0; i < 3; i++)
  {
    matrix->SetElement(i, 0, nx[i]);
    matrix->SetElement(i, 1, ny[i]);
    matrix->SetElement(i, 2, nz[i]);
  }
  t->Concatenate(matrix);
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

void vtkProsthesisRepresentation::SetHandleColour(double* rgb)
{
  this->HandleProperty->SetColor(rgb);
}

void vtkProsthesisRepresentation::SetHandleColour(double r, double g, double b)
{
  this->HandleProperty->SetColor(r, g, b);

}

void vtkProsthesisRepresentation::SetSelectedHandleColour(double* rgb)
{
  this->SelectedHandleProperty->SetColor(rgb);
}

void vtkProsthesisRepresentation::SetSelectedHandleColour(double r, double g, double b)
{
  this->SelectedHandleProperty->SetColor(r, g, b);

}

void vtkProsthesisRepresentation::SetRadius(double value)
{
  if (value <= 0.0) value = 0.01;
  this->Radius = value;
  this->OutlineGeometry->SetRadius(value);
  this->PositionHandles();
}
