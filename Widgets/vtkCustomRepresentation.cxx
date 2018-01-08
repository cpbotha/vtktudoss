#include "vtkCustomRepresentation.h"
#include "vtkActor.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
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


vtkStandardNewMacro(vtkCustomRepresentation);

//----------------------------------------------------------------------------
vtkCustomRepresentation::vtkCustomRepresentation() :
  ShowOutline(true)
{
  // The initial state
  this->InteractionState = vtkCustomRepresentation::Outside;

  // Handle size is in pixels for this widget
  this->HandleSize = 5.0;

  // Set up the initial properties
  this->CreateDefaultProperties();

  // Create the handles
  this->HandleGeometry = vtkSphereSource::New();
  this->HandleGeometry->SetThetaResolution(16);
  this->HandleGeometry->SetPhiResolution(8);
  this->HandleMapper = vtkPolyDataMapper::New();
  this->HandleMapper->SetInputConnection(this->HandleGeometry->GetOutputPort());
  this->Handle = vtkActor::New();
  this->Handle->SetMapper(this->HandleMapper);
  this->Handle->SetProperty(this->HandleProperty);

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
  vtkCellArray* cells = vtkCellArray::New();
  cells->Allocate(cells->EstimateSize(15,2));
  this->OutlinePolyData->SetLines(cells);
  cells->Delete();

  // Create the outline
  this->GenerateOutline();

  // Define the point coordinates
  double bounds[6] = {-1.0, 1.0, 
                      -1.0, 1.0,
                      -1.0, 1.0};
  this->PlaceWidget(bounds);

  //Manage the picking stuff
  this->HandlePicker = vtkCellPicker::New();
  this->HandlePicker->SetTolerance(0.001);
  this->HandlePicker->AddPickList(this->Handle);
  this->HandlePicker->PickFromListOn();
}

//----------------------------------------------------------------------------
vtkCustomRepresentation::~vtkCustomRepresentation()
{
  this->HandleGeometry->Delete();
  this->HandleMapper->Delete();
  this->Handle->Delete();
  this->HandlePicker->Delete();
  this->HandleProperty->Delete();
  this->SelectedHandleProperty->Delete();
  this->Points = vtkPoints::New(VTK_DOUBLE);
  this->Points->SetNumberOfPoints(8); // 8 corners of the bounds

  this->OutlinePolyData->Delete();
  this->OutlineMapper->Delete();
  this->Outline->Delete();

  this->Points->Delete();
}

//----------------------------------------------------------------------
void vtkCustomRepresentation::StartWidgetInteraction(double e[2])
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
void vtkCustomRepresentation::WidgetInteraction(double e[2])
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
  this->HandlePicker->GetPickPosition(pos);
  vtkInteractorObserver::ComputeWorldToDisplay(this->Renderer,
                                               pos[0], pos[1], pos[2],
                                               focalPoint);
  z = focalPoint[2];
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, this->LastEventPosition[0],
                                               this->LastEventPosition[1], z, prevPickPoint);
  vtkInteractorObserver::ComputeDisplayToWorld(this->Renderer, e[0], e[1], z, pickPoint);

  if (this->InteractionState == vtkCustomRepresentation::Translating)
  {
    this->Translate(prevPickPoint, pickPoint);
  }

  // Store the start position
  this->LastEventPosition[0] = e[0];
  this->LastEventPosition[1] = e[1];
  this->LastEventPosition[2] = 0.0;
}

//----------------------------------------------------------------------------
void vtkCustomRepresentation::Translate(double *p1, double *p2)
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
  for (int i = 0; i < 8; i++) {
    *pts++ += v[0];
    *pts++ += v[1];
    *pts++ += v[2];
  }

  this->PositionHandles();
}

//----------------------------------------------------------------------------
void vtkCustomRepresentation::CreateDefaultProperties()
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
void vtkCustomRepresentation::PlaceWidget(double bds[6])
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
  this->InitialLength = sqrt((bounds[1]-bounds[0])*(bounds[1]-bounds[0]) +
                             (bounds[3]-bounds[2])*(bounds[3]-bounds[2]) +
                             (bounds[5]-bounds[4])*(bounds[5]-bounds[4]));

  this->PositionHandles();
  this->ValidPick = 1; //since we have set up widget
  this->SizeHandles();
}

//----------------------------------------------------------------------------
int vtkCustomRepresentation::ComputeInteractionState(int X, int Y, int modify)
{
  // Okay, we can process this. Try to pick handles first;
  // if no handles picked, then pick the bounding box.
  if (!this->Renderer || !this->Renderer->IsInViewport(X, Y))
  {
    this->InteractionState = vtkCustomRepresentation::Outside;
    return this->InteractionState;
  }

  vtkAssemblyPath* path = this->GetAssemblyPath(X, Y, 0., this->HandlePicker);

  if ( path != NULL )
  {
    this->ValidPick = 1;
    this->Handle = reinterpret_cast<vtkActor*>(path->GetFirstNode()->GetViewProp());
    this->InteractionState = vtkCustomRepresentation::Translating;
  }
  else {
    this->InteractionState = vtkCustomRepresentation::Outside;
  }

  return this->InteractionState;
}

//----------------------------------------------------------------------------
void vtkCustomRepresentation::BuildRepresentation()
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
void vtkCustomRepresentation::ReleaseGraphicsResources(vtkWindow *w)
{
  // render the handles
  this->Handle->ReleaseGraphicsResources(w);
  this->Outline->ReleaseGraphicsResources(w);
}

//----------------------------------------------------------------------------
int vtkCustomRepresentation::RenderOpaqueGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();
  count += this->Outline->RenderOpaqueGeometry(v);
  // render the handles
  if(this->Handle->GetVisibility())
  {
    count += this->Handle->RenderOpaqueGeometry(v);
  }

  return count;
}

//----------------------------------------------------------------------------
int vtkCustomRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport *v)
{
  int count=0;
  this->BuildRepresentation();

  count += this->Outline->RenderTranslucentPolygonalGeometry(v);
  // render the handles
  if(this->Handle->GetVisibility())
  {
    count += this->Handle->RenderTranslucentPolygonalGeometry(v);
  }

  return count;
}

//----------------------------------------------------------------------------
int vtkCustomRepresentation::HasTranslucentPolygonalGeometry()
{
  int result=0;
  this->BuildRepresentation();

  result |= this->Outline->HasTranslucentPolygonalGeometry();
  result |= this->Handle->HasTranslucentPolygonalGeometry();

  return result;
}

//----------------------------------------------------------------------------
void vtkCustomRepresentation::PositionHandles()
{
  this->HandleGeometry->SetCenter(this->Center);
  this->GenerateOutline();
}

//----------------------------------------------------------------------------
void vtkCustomRepresentation::SizeHandles()
{
  double radius =
      this->vtkWidgetRepresentation::SizeHandlesInPixels(1.5, this->Center);
  this->HandleGeometry->SetRadius(radius);
}

void vtkCustomRepresentation::HighlightHandle(vtkProp* prop)
{
  if (prop == NULL)
  {
    this->Handle->SetProperty(this->HandleProperty);
    return;
  }

  this->Handle = static_cast<vtkActor*>(prop);
  if (this->Handle)
  {
    this->Handle->SetProperty(this->SelectedHandleProperty);
  }
}

//----------------------------------------------------------------------------
void vtkCustomRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void vtkCustomRepresentation::SetInteractionState(int state)
{
  // Clamp to allowable values
  state = (state < vtkCustomRepresentation::Outside ? vtkCustomRepresentation::Outside :
          (state > vtkCustomRepresentation::Scaling ? vtkCustomRepresentation::Scaling : state));

  this->InteractionState = state;
  switch (state)
  {
    case vtkCustomRepresentation::Translating:
      this->HighlightHandle(this->Handle);
      break;
    case vtkCustomRepresentation::Rotating:
    case vtkCustomRepresentation::Scaling:
    default:
      this->HighlightHandle(NULL);
  }
}

void vtkCustomRepresentation::SetShowOutline(bool show) {
  if (this->ShowOutline != show) {
    this->ShowOutline = show;
    this->GenerateOutline();
  }
}

void vtkCustomRepresentation::GenerateOutline()
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