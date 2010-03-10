/*=========================================================================
vtkAdvancedAngleWidget
=========================================================================*/

#include "vtk3DWidget.h"

class vtkActor;
class vtkArrowSource;
class vtkAssembly;
class vtkCellPicker;
//class vtkPointPicker;
class vtkDataSetMapper;
class vtkFloatArray;
class vtkPlanes;
class vtkPoints;
class vtkPolyData;
class vtkRenderer;
class vtkPolyDataMapper;
class vtkProp;
class vtkProperty;
class vtkSphereSource;
class vtkSTLReader;
class vtkTextActor3D;
class vtkTransform;
class vtkTriangleStrip;

#ifndef __vtkAdvancedAngleWidget_h
#define __vtkAdvancedAngleWidget_h

class VTK_EXPORT vtkAdvancedAngleWidget : public vtk3DWidget
{
public:
	void PrintSelf(ostream& os, vtkIndent indent);
	vtkTypeRevisionMacro(vtkAdvancedAngleWidget,vtk3DWidget);
	static vtkAdvancedAngleWidget *New();

	double widgetcolor[3];
	double pos[3];
	double posoffset[3];
	double center[3];	
	double initialaxis[3];
	double widgetaxis[3];
	double bounds[6];
	double minArcPoint[3];
	double maxArcPoint[3];
	double totalangle;
	double Scale;
	double mirrorfonts;
	int upperlowerbool;
	vtkTransform *orientation;
	vtkTransform *widgetorientation;
	double WidgetScale;
	double theta; //rotation angle
	double rotation_axis[3]; //axis of rotation
	bool showContext;
	bool piechart;

	double filterRangeMin; // 0 < x < 180
	double filterRangeMax; // 0 < x < 180
	double absoluteRangeMin1;
	double absoluteRangeMax1;
	double absoluteRangeMin2;
	double absoluteRangeMax2;
	double axis_length;
	double initial_axis_length;

	int arcparts;

	// Description:
	// Methods that satisfy the superclass' API.
	virtual void SetEnabled(int);
	virtual void PlaceWidget(double bds[6]);
	virtual void SetInitialPosition(double position[3]);
	virtual void SetInitialAxis(double widgetorientation[3]);
	virtual void Rotate(double angle, double axis[3]);
	virtual void SetColor(double color[3]);
	double GetXOffset();
	double GetYOffset();
	double GetZOffset();
	void GetOriginalPosition(double point[3]);
	void GetRotationTransform(vtkTransform *t);	
	void PlaceWidget()
		{this->Superclass::PlaceWidget();}
	void PlaceWidget(double xmin, double xmax, double ymin, double ymax, 
					double zmin, double zmax)
		{this->Superclass::PlaceWidget(xmin,xmax,ymin,ymax,zmin,zmax);}

	virtual void Realign();
	void CreateDefaultProperties();
	void SetWidgetRange1(double min, double max);
	void SetWidgetRange2(double min, double max);
	void UpdateRangeIndicators();
	void ToggleSize();

	vtkSetMacro(TranslationEnabled,int);
	vtkGetMacro(TranslationEnabled,int);
	vtkBooleanMacro(TranslationEnabled,int);
	vtkSetMacro(showContext, bool);
	vtkGetMacro(showContext, bool);
	vtkSetMacro(piechart, bool);
	vtkGetMacro(piechart, bool);
	vtkSetMacro(ScalingEnabled,int);
	vtkGetMacro(ScalingEnabled,int);
	vtkSetMacro(LabelsEnabled, int);
	vtkGetMacro(LabelsEnabled, int);
	vtkBooleanMacro(LabelsEnabled, int);
	vtkBooleanMacro(ScalingEnabled,int);
	vtkSetMacro(RotationEnabled,int);
	vtkGetMacro(RotationEnabled,int);
	vtkBooleanMacro(RotationEnabled,int);
	vtkSetVectorMacro(center,double,3);
	vtkGetVectorMacro(center,double,3);
	vtkSetVectorMacro(pos,double,3);
	vtkGetVectorMacro(pos,double,3);
	vtkSetVectorMacro(posoffset,double,3);
	vtkGetVectorMacro(posoffset,double,3);
	vtkSetMacro(filterRangeMin, double);
	vtkSetMacro(filterRangeMax, double);
	vtkGetMacro(filterRangeMin, double);
	vtkGetMacro(filterRangeMax, double);
	vtkGetMacro(mirrorfonts, double);
	vtkSetMacro(mirrorfonts, double);
	vtkGetMacro(totalangle,double);
	vtkSetMacro(totalangle,double);
	vtkGetVectorMacro(initialaxis, double,3);
	vtkSetMacro(Scale,double);
	vtkGetMacro(Scale,double);
	vtkSetMacro(axis_length,double);
	vtkGetMacro(axis_length,double);	
	vtkSetVectorMacro(widgetcolor,double,3);
	vtkGetVectorMacro(widgetcolor,double,3);

	vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

	virtual void GetTransform(vtkTransform *t);
	virtual void GetTransformForActor(vtkTransform *t);
	virtual void GetInitialOrientationTransform(vtkTransform *t);

protected:
    vtkAdvancedAngleWidget(void);
	~vtkAdvancedAngleWidget(void);
	//BTX - manage the state of the widget
	char *FileName;

	int State;
	enum WidgetState
	{
		Start=0,
		Moving,
		Scaling,
		Outside
	};
	//ETX
	    
	// Handles the events
	static void ProcessEvents(vtkObject* object, 
								unsigned long event,
								void* clientdata, 
								void* calldata);
	virtual bool IsFlat();

	// ProcessEvents() dispatches to these methods.
	virtual void OnMouseMove();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();

	virtual void SubtractVectors(double vector1[3], double vector2[3], double *result);
	virtual void AddVectors(double vector1[3], double vector2[3], double *result);

	//virtual void OnMiddleButtonDown();
	//virtual void OnMiddleButtonUp();
	//virtual void OnRightButtonDown();
	//virtual void OnRightButtonUp();

	vtkRenderer		  *Renderer;

	vtkActor		  *clippingDummy1;
	vtkPolyDataMapper *clippingDummyMapper1;
	vtkSphereSource   *clippingDummySphere1;

	vtkAssembly		  *axisAssembly;
	vtkPolyData		  *ArcZ1;

	vtkPoints            *innerArcPoints;
	vtkActor		     *innerArcA;
	vtkDataSetMapper     *innerArcMapA;	
	vtkPolyData			 *innerboogje1;
	vtkActor		     *innerArcB;
	vtkDataSetMapper     *innerArcMapB;	
	vtkPolyData			 *innerboogje2;

	vtkPoints            *outerArcPoints;
	vtkActor		     *outerArcA;
	vtkDataSetMapper     *outerArcMapA;
	vtkPolyData			 *outerboogje1;
	vtkActor		     *outerArcB;
	vtkDataSetMapper     *outerArcMapB;	
	vtkPolyData			 *outerboogje2;

	vtkPoints            *outer2ArcPoints;
	vtkActor		     *outer2ArcA;
	vtkDataSetMapper     *outer2ArcMapA;
	vtkPolyData			 *outer2boogje1;
	vtkActor		     *outer2ArcB;
	vtkDataSetMapper     *outer2ArcMapB;	
	vtkPolyData			 *outer2boogje2;

	vtkSTLReader		 *stlreader;
	vtkDataSetMapper     *nudgemapper;	
	vtkActor		     *nudgeactor;



	vtkTextActor3D *textActorMin;
	vtkTextActor3D *textActorMax;

	// glyphs representing hot spots (e.g., handles)
	vtkActor          **Handle;
	vtkPolyDataMapper **HandleMapper;
	vtkSphereSource	  **HandleGeometry;
	//virtual void PositionHandles();
	int HighlightHandle(vtkProp *prop); //returns cell id

	// Do the picking
	vtkCellPicker *HandlePicker;
	//vtkPointPicker *HandlePicker;
	vtkActor *CurrentHandle;
	int      CurrentAxis;

	//virtual void Translate(double *p1, double *p2);
	// Transform the hexahedral points (used for rotations)
	vtkTransform *Transform;

	vtkProperty *HandlePropertyDarkBlue;
	vtkProperty *HandlePropertyBlue;
	vtkProperty *HandlePropertyLightBlue;

	vtkProperty *HandlePropertyDarkYellow;
	vtkProperty *HandlePropertyYellow;
	vtkProperty *HandlePropertyLightYellow;

	vtkProperty *HandlePropertyOrange;

	vtkProperty *HandlePropertyNudge;

	vtkProperty *HandlePropertyHighlight;
	vtkProperty *HandlePropertyInvalid;
	vtkProperty *HandlePropertyValid;
	
	void CreatePrimaryArc(float min, float max);
	void CreateArcRange1(float min, float max);
	void CreateArcRange2(float min, float max);
	void CreateNudge();
	void DestroyPrimaryArc();
	void DestroyArcRange1();
	void DestroyArcRange2();

	void AxisTranslate(int axis, double *p1, double *p2);
	void GetDirection(const double Nx[3],const double Ny[3], const double Nz[3], double dir[3]);
	double AxisRotate(int axisconstraint, double *p1, double *p2);
	bool IllegalPosition();
	//void HighlightOutline(int highlight);

	// Control whether scaling, rotation, and translation are supported
	int TranslationEnabled;
	int ScalingEnabled;
	int RotationEnabled;
	int LabelsEnabled;
};


#endif




