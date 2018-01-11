#ifndef vtkProsthesisRepresentation_h
#define vtkProsthesisRepresentation_h

#include "vtkInteractionWidgetsModule.h" // For export macro
#include "vtkWidgetRepresentation.h"

class vtkActor;
class vtkPolyDataMapper;
class vtkArrowSource;
class vtkRegularPolygonSource;
class vtkCellPicker;
class vtkProperty;
class vtkPoints;
class vtkPolyData;
class vtkTransform;
class vtkAppendPolyData;
class vtkTransformPolyDataFilter;
class vtkAlgorithmOutput;


class VTKINTERACTIONWIDGETS_EXPORT vtkProsthesisRepresentation : public vtkWidgetRepresentation
{
public:
  // Instantiate the class.
  static vtkProsthesisRepresentation *New();

  // Standard methods for the class.
  vtkTypeMacro(vtkProsthesisRepresentation,vtkWidgetRepresentation);
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

  virtual void GetTransform(vtkTransform *t);

  double Center[3];
  vtkSetVectorMacro(Center, double, 3);
  vtkGetVectorMacro(Center, double, 3);

  double Radius;
  void SetRadius(double value);
  vtkGetMacro(Radius, double);

  double HandleColour[3];
  void SetHandleColour(double*);
  void SetHandleColour(double, double, double);
  vtkGetVector3Macro(HandleColour, double);

  double SelectedHandleColour[3];
  void SetSelectedHandleColour(double*);
  void SetSelectedHandleColour(double, double, double);
  vtkGetVector3Macro(SelectedHandleColour, double);

  // Flag used to show/hide the handles of the widget
  bool ShowHandles;
  void SetShowHandles(bool);
  vtkGetMacro(ShowHandles, bool);
  vtkBooleanMacro(ShowHandles, bool);

  // Flag used to show/hide the outline of the widget
  bool ShowOutline;
  void SetShowOutline(bool);
  vtkGetMacro(ShowOutline, bool);
  vtkBooleanMacro(ShowOutline, bool);

  double* GetBounds();

  // The interaction state may be set from a widget (e.g., vtkCustomWidget) or
  // other object. This controls how the interaction with the widget
  // proceeds. Normally this method is used as part of a handshaking
  // process with the widget: First ComputeInteractionState() is invoked that
  // returns a state based on geometric considerations (i.e., cursor near a
  // widget feature), then based on events, the widget may modify this further.
  void SetInteractionState(int state);

  // Used to manage the state of the widget.
  enum {Outside=0, Translating, Rotating, UpClick, DownClick};

protected:
  // Constructor and destructor.
  vtkProsthesisRepresentation();
  ~vtkProsthesisRepresentation() VTK_OVERRIDE;

  // Manage how the representation appears.
  double LastEventPosition[3];

  // Important points on the widget.
  vtkPoints* Points;

  // The translation handle.
  vtkActor* Handle;
  vtkPolyDataMapper* HandleMapper;
  vtkAppendPolyData* TranslateGeometryCombiner;
  vtkRegularPolygonSource* HandleGeometry;
  vtkPolyData* TranslateArrowPolyData;
  void GenerateTranslateArrow(double shaftWidth = 0.05);

  // The up handles.
  vtkActor* UpHandle;
  vtkPolyDataMapper* UpHandleMapper;
  vtkPolyData* PlusPolyData;
  void GeneratePlus(double size, bool outline = false);

  // The down handle.
  vtkActor* DownHandle;
  vtkPolyDataMapper* DownHandleMapper;
  vtkPolyData* MinusPolyData;
  void GenerateMinus(double size, bool outline = false);

  // The rotation handle.
  vtkActor* RotateHandle;
  vtkPolyDataMapper* RotateHandleMapper;
  vtkAppendPolyData* RotateGeometryCombiner;
  vtkRegularPolygonSource* RotateHandleGeometry;
  vtkPolyData* LeftArrowPolyData;
  vtkPolyData* RightArrowPolyData;

  // GenerateArrow parameters:
  // shaftWidth : The thickness of the shaft.
  // clockwise : Determines which way the arrow points.
  // outline : If true, only draws the arrow's outline.
  void GenerateArrow(vtkPolyData* arrowPolyData,
                     double shaftWidth = 0.05, bool clockwise = true, 
                     bool outline = false);

  void HighlightHandle(vtkProp *prop);
  virtual void PositionHandles();
  virtual void UpdateHandles();

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
  vtkRegularPolygonSource* OutlineGeometry;
  vtkPolyDataMapper* OutlineMapper;

  // Methods to update the widget
  virtual void Translate(double *p1, double *p2);
  virtual void Rotate(double previousX, double previousY, double X, double Y, double *vpn);

  vtkTransform* Transform;
  void UpdateTransform();
  void GetAlwaysUpTransform(vtkTransform* t);

private:
  vtkProsthesisRepresentation(const vtkProsthesisRepresentation&) VTK_DELETE_FUNCTION;
  void operator=(const vtkProsthesisRepresentation&) VTK_DELETE_FUNCTION;

  void TransformPolyData(vtkAlgorithmOutput* algoOutput, vtkPolyDataMapper* mapper);
  void TransformPolyData(vtkPolyData* polyData, vtkPolyDataMapper* mapper);
};

#endif
