#ifndef vtkCustomRepresentation_h
#define vtkCustomRepresentation_h

#include "vtkInteractionWidgetsModule.h" // For export macro
#include "vtkWidgetRepresentation.h"

class vtkActor;
class vtkPolyDataMapper;
class vtkSphereSource;
class vtkCellPicker;
class vtkProperty;
class vtkPoints;
class vtkPolyData;


class VTKINTERACTIONWIDGETS_EXPORT vtkCustomRepresentation : public vtkWidgetRepresentation
{
public:
  // Instantiate the class.
  static vtkCustomRepresentation *New();

  // Standard methods for the class.
  vtkTypeMacro(vtkCustomRepresentation,vtkWidgetRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;

  // These are methods that satisfy vtkWidgetRepresentation's API.
  void PlaceWidget(double bounds[6]) VTK_OVERRIDE;
  void BuildRepresentation() VTK_OVERRIDE;
  int  ComputeInteractionState(int X, int Y, int modify=0) VTK_OVERRIDE;
  void StartWidgetInteraction(double e[2]) VTK_OVERRIDE;
  void WidgetInteraction(double e[2]) VTK_OVERRIDE;

  // Methods supporting, and required by, the rendering process.
  void ReleaseGraphicsResources(vtkWindow*) VTK_OVERRIDE;
  int  RenderOpaqueGeometry(vtkViewport*) VTK_OVERRIDE;
  int  RenderTranslucentPolygonalGeometry(vtkViewport*) VTK_OVERRIDE;
  int  HasTranslucentPolygonalGeometry() VTK_OVERRIDE;

  double Center[3];
  vtkSetVector3Macro(Center, double);
  vtkGetVector3Macro(Center, double);

  bool ShowOutline;
  void SetShowOutline(bool);
  vtkGetMacro(ShowOutline, bool);
  vtkBooleanMacro(ShowOutline, bool);

  // The interaction state may be set from a widget (e.g., vtkCustomWidget) or
  // other object. This controls how the interaction with the widget
  // proceeds. Normally this method is used as part of a handshaking
  // process with the widget: First ComputeInteractionState() is invoked that
  // returns a state based on geometric considerations (i.e., cursor near a
  // widget feature), then based on events, the widget may modify this further.
  void SetInteractionState(int state);

  // Used to manage the state of the widget.
  enum {Outside=0, Translating, Rotating, Scaling};

protected:
  // Constructor and destructor.
  vtkCustomRepresentation();
  ~vtkCustomRepresentation() VTK_OVERRIDE;

  // Manage how the representation appears.
  double LastEventPosition[3];

  // Important points on the widget.
  vtkPoints* Points;

  // The translation handle.
  vtkActor* Handle;
  vtkPolyDataMapper* HandleMapper;
  vtkSphereSource* HandleGeometry;
  void HighlightHandle(vtkProp *prop);
  virtual void PositionHandles();
  virtual void SizeHandles();

  // Do the picking
  vtkCellPicker *HandlePicker;

  // Properties used to control the appearance of selected objects and
  // the manipulator in general.
  vtkProperty* HandleProperty;
  vtkProperty* SelectedHandleProperty;
  vtkProperty* OutlineProperty;
  virtual void CreateDefaultProperties();

  // wireframe outline
  vtkActor* Outline;
  vtkPolyDataMapper* OutlineMapper;
  vtkPolyData* OutlinePolyData;
  // Generate an outline on the bounds of the widget
  void GenerateOutline();

  // Methods to update the widget
  virtual void Translate(double *p1, double *p2);

private:
  vtkCustomRepresentation(const vtkCustomRepresentation&) VTK_DELETE_FUNCTION;
  void operator=(const vtkCustomRepresentation&) VTK_DELETE_FUNCTION;
};

#endif
