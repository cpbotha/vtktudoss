

/*=========================================================================
AxesWidget : test class
=========================================================================*/
#include "vtkAdvancedAngleWidget.h"

#include "math.h"

#include "vtkActor.h"
#include "vtkArrowSource.h"
#include "vtkAssemblyNode.h"
#include "vtkAssemblyPath.h"
#include "vtkAssembly.h"
#include "vtkCallbackCommand.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellPicker.h"
#include "vtkPointPicker.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlanes.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPointData.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRendererCollection.h"
#include "vtkSphereSource.h"
#include "vtkSTLReader.h"
#include "vtkTextActor3D.h"
#include "vtkTextProperty.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkTriangleStrip.h"
#include "vtkViewport.h"
#include "vtkDataSetMapper.h"

// we're using snprintf to convert from integer to string representation
#ifdef _MSC_VER  //all MS compilers define this (version)
     #define snprintf _snprintf
#endif

vtkCxxRevisionMacro(vtkAdvancedAngleWidget, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkAdvancedAngleWidget);

// Constructor
vtkAdvancedAngleWidget::vtkAdvancedAngleWidget()
{
	int i, j;
	this->filterRangeMin = 90.0;
	this->filterRangeMax = 90.0;
	this->absoluteRangeMin1 = 90.0;
	this->absoluteRangeMax1 = 90.0;
	this->widgetcolor[0] = 0.0;
	this->widgetcolor[1] = 0.0;
	this->widgetcolor[2] = 1.0;
	this->axis_length = 0.8;
	this->mirrorfonts = 1.0;
	this->piechart = true;
	
	this->innerArcA = vtkActor::New();
	this->innerArcB = vtkActor::New();
	this->outerArcA = vtkActor::New();
	this->outerArcB = vtkActor::New();
	this->outer2ArcA = vtkActor::New();
	this->outer2ArcB = vtkActor::New();
	
	this->nudgeactor = vtkActor::New();


    // we have to init these as well so Peter's destruction in the
    // the dtor doesn't break everything...
    this->innerArcMapA = NULL;
	this->innerArcMapB = NULL;
	this->outerArcMapA = NULL;
	this->outerArcMapB = NULL;
	this->outer2ArcMapA = NULL;
	this->outer2ArcMapB = NULL;


	this->HandlePropertyDarkBlue = vtkProperty::New();
	this->HandlePropertyBlue = vtkProperty::New();
	this->HandlePropertyLightBlue = vtkProperty::New();
	
	this->HandlePropertyDarkYellow = vtkProperty::New();
	this->HandlePropertyYellow = vtkProperty::New();
	this->HandlePropertyLightYellow = vtkProperty::New();

	this->HandlePropertyOrange = vtkProperty::New();
	this->HandlePropertyNudge = vtkProperty::New();

    this->HandlePropertyHighlight = vtkProperty::New();
	this->CreateDefaultProperties();

	this->State = vtkAdvancedAngleWidget::Start;
	this->EventCallbackCommand->SetCallback(vtkAdvancedAngleWidget::ProcessEvents);
	// Enable/disable the translation, rotation, and scaling of the widget
	this->TranslationEnabled = 1;
	this->RotationEnabled = 0;
	this->ScalingEnabled = 0;
	this->LabelsEnabled = 0;

	/*this->Renderer = vtkRenderer::New();
	this->Renderer->SetLayer( 1 );
	this->Renderer->InteractiveOff();*/

	this->Scale = 1.0;

	this->bounds[0] = -0.5;
	this->bounds[1] = 0.5;
	this->bounds[2] = -0.5;
	this->bounds[3] = 0.5;
	this->bounds[4] = -0.5;
	this->bounds[5] = 0.5;

	this->pos[0] = 0.0;
	this->pos[1] = 0.0;
	this->pos[2] = 0.0;
	this->posoffset[0] = 0.0;
	this->posoffset[1] = 0.0;
	this->posoffset[2] = 0.0;
	this->orientation = vtkTransform::New();
	this->widgetorientation = vtkTransform::New();
	this->orientation->Identity();
	this->widgetorientation->Identity();
	this->WidgetScale = 1.0;

	vtkTextProperty *tp1 = vtkTextProperty::New();
	this->textActorMin = vtkTextActor3D::New();
	this->textActorMax = vtkTextActor3D::New();
	tp1 = this->textActorMin->GetTextProperty();
    tp1->SetFontFamilyToArial();
	tp1->SetFontSize(48);
    tp1->ShadowOn();
    tp1->BoldOn();
    tp1->SetJustificationToRight();
    tp1->SetVerticalJustificationToBottom();
	this->textActorMax->SetTextProperty(tp1);
    //textActor->GetTextProperty()->SetColor(scene->fontColor);

    //this->CurrentRenderer->AddActor2D(textActor);

	//this->CurrentRenderer->AddActor2D(minText);


	//Manage the picking stuff
	this->HandlePicker = vtkCellPicker::New();
	//this->HandlePicker = vtkPointPicker::New();
	this->HandlePicker->SetTolerance(0.001);
	this->HandlePicker->AddPickList(this->innerArcA);	
	this->HandlePicker->AddPickList(this->nudgeactor);
	this->HandlePicker->PickFromListOn();

	this->CurrentHandle = NULL;

	this->Transform = vtkTransform::New();

	this->FileName = NULL;
}



// Destructor
vtkAdvancedAngleWidget::~vtkAdvancedAngleWidget()
{	
    // thank you cpbotha, segfault preventor!
	
    // instantiated in CreatePrimaryArc
    if (this->innerArcMapA)
        this->innerArcMapA->Delete();
    // instantiated in ctor
	this->innerArcA->Delete();
    // instantiated in CreatePrimaryArc
    if (this->innerArcMapB)
        this->innerArcMapB->Delete();
	// instantiated in ctor
    this->innerArcB->Delete();
    // instantiated in CreateArcRange1
    if (this->outerArcMapA)
        this->outerArcMapA->Delete();
    // instantiated in ctor
	this->outerArcA->Delete();
    // instantiated in CreateArgRange1
    if (this->outerArcMapB)
        this->outerArcMapB->Delete();
    // instantiated in ctor
	this->outerArcB->Delete();
    // instantiated in CreateArcRange2
    if (this->outer2ArcMapA)
        this->outer2ArcMapA->Delete();
    // instantiated in ctor
	this->outer2ArcA->Delete();
    // instantiated in CreateArcRange2
    if (this->outer2ArcMapB)
        this->outer2ArcMapB->Delete();
    // instantiated in ctor
	this->outer2ArcB->Delete();
    // instantiated in ctor
	this->nudgeactor->Delete();

	this->Transform->Delete();
    this->HandlePicker->Delete();
	this->SetFileName(0);
}



//	Range Indicator 1
void vtkAdvancedAngleWidget::CreateArcRange1(float min, float max) 
{
	//arcZ
	float Pi;
    float a1[3], a2[3];
	double arcparts1, arcparts2, arcparts3;
	int i, j;

	if (min < 0){
		arcparts1 = 0;
	}
	else
		arcparts1 = (int)((min)/2);
	arcparts2 = (int)((max - min)/2);
	if (max>180)
		arcparts3 = 0;
	else
		arcparts3 = (int)(90-(max)/2);

	Pi = 3.141593;
	this->rotation_axis[0]=1;this->rotation_axis[1]=0;this->rotation_axis[2]=0;
	this->theta=0;

	a1[0]=0.0;a1[1]=0.0;a1[2]=1.0;
	a2[0]=1.0;a2[1]=0.0;a2[2]=0.0;
	this->outerArcPoints = vtkPoints::New();
	this->outerArcPoints->SetNumberOfPoints((arcparts1+arcparts2 + arcparts3 + 3) * 2 +2);

	float filterRangeMin_dec, filterRangeMax_dec;
	filterRangeMin_dec = min / 180.0;
	filterRangeMax_dec = max / 180.0;

	this->outerArcPoints->InsertPoint(0, 0.0, 0.1, 0.0);
	this->outerArcPoints->InsertPoint(1, 0.0, -0.1, 0.0);
	for(i=0; i < arcparts1+1; i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.07)*sin( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v2 = (this->axis_length-0.07)*cos( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v3 = (this->axis_length-0.03)*sin( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v4 = (this->axis_length-0.03)*cos( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));

		this->outerArcPoints->InsertPoint(2+(i*2), v2, 0, v1);
		this->outerArcPoints->InsertPoint(2+(i*2)+1, v4, 0, v3);
    }

	for(i=arcparts1+1; i < (arcparts1+1 + arcparts2+1); i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.1)*sin( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v2 = (this->axis_length-0.1)*cos( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v3 = this->axis_length*sin( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v4 = this->axis_length*cos( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));

		this->outerArcPoints->InsertPoint(2+(i*2), v2, 0, v1);
		this->outerArcPoints->InsertPoint(2+(i*2)+1, v4, 0, v3);
    }

	for(i=arcparts1+1+arcparts2+1; i < (arcparts1+1 + arcparts2+1 + arcparts3+1); i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.07)*sin( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v2 = (this->axis_length-0.07)*cos( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v3 = (this->axis_length-0.03)*sin( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v4 = (this->axis_length-0.03)*cos( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));

		this->outerArcPoints->InsertPoint(2+(i*2), v2, 0, v1);
		this->outerArcPoints->InsertPoint(2+(i*2)+1, v4, 0, v3);
    }



	////printf("%f %f %f\n", this->minArcPoint[0], this->minArcPoint[1], this->minArcPoint[2]);
	////printf("%f %f %f\n", this->maxArcPoint[0], this->maxArcPoint[1], this->maxArcPoint[2]);

	vtkCellArray *cells2A = vtkCellArray::New();	
	vtkCellArray *cells2B = vtkCellArray::New();	
	j=1;
	for(i=0; i<arcparts1; i++)
	{	
		cells2A->InsertNextCell(4);		
		cells2A->InsertCellPoint(j*2);			
		cells2A->InsertCellPoint(j*2+1);
		cells2A->InsertCellPoint(j*2+3);
		cells2A->InsertCellPoint(j*2+2);
		j=j+1;
	}
	j=j+1;

	for(i=0; i<arcparts2; i++)
	{	
		cells2B->InsertNextCell(4);		
		cells2B->InsertCellPoint(j*2);			
		cells2B->InsertCellPoint(j*2+1);
		cells2B->InsertCellPoint(j*2+3);
		cells2B->InsertCellPoint(j*2+2);
		j=j+1;
	}
	j=j+1;

	for(i=0; i<arcparts3; i++)
	{	
		cells2A->InsertNextCell(4);		
		cells2A->InsertCellPoint(j*2);			
		cells2A->InsertCellPoint(j*2+1);
		cells2A->InsertCellPoint(j*2+3);
		cells2A->InsertCellPoint(j*2+2);
		j=j+1;
	}

	//printf("%i", j-1);


	this->outerboogje1 = vtkPolyData::New();
	this->outerboogje1->SetPolys(cells2A);	
	this->outerboogje1->SetPoints(this->outerArcPoints);
	cells2A->Delete();

	this->outerboogje2 = vtkPolyData::New();
	this->outerboogje2->SetPolys(cells2B);	
	this->outerboogje2->SetPoints(this->outerArcPoints);
	cells2B->Delete();
	
	this->outerArcMapA = vtkDataSetMapper::New();
	this->outerArcMapA->SetInput(this->outerboogje1);
	this->outerArcA->SetMapper(this->outerArcMapA);
	this->outerArcA->SetProperty(this->HandlePropertyDarkBlue);
	this->outerArcMapB = vtkDataSetMapper::New();
	this->outerArcMapB->SetInput(this->outerboogje2);
	this->outerArcB->SetMapper(this->outerArcMapB);
	this->outerArcB->SetProperty(this->HandlePropertyLightBlue);
}

void vtkAdvancedAngleWidget::DestroyArcRange1()
{
	this->outerArcPoints->Delete();
    this->outerboogje1->Delete();
	this->outerboogje2->Delete();
	this->outerArcMapA->Delete();
	this->outerArcMapB->Delete();
}

//	Range Indicator 1
void vtkAdvancedAngleWidget::CreateArcRange2(float min, float max) 
{
	//arcZ
	float Pi;
    float a1[3], a2[3];
	double arcparts1, arcparts2, arcparts3;
	int i, j;

	if (min < 0){
		arcparts1 = 0;
	}
	else
		arcparts1 = (int)((min)/2);
	arcparts2 = (int)((max - min)/2);
	if (max>180)
		arcparts3 = 0;
	else
		arcparts3 = (int)(90-(max)/2);

	Pi = 3.141593;
	this->rotation_axis[0]=1;this->rotation_axis[1]=0;this->rotation_axis[2]=0;
	this->theta=0;

	a1[0]=0.0;a1[1]=0.0;a1[2]=1.0;
	a2[0]=1.0;a2[1]=0.0;a2[2]=0.0;
	this->outer2ArcPoints = vtkPoints::New();
	this->outer2ArcPoints->SetNumberOfPoints((arcparts1+arcparts2 + arcparts3 + 3) * 2 +2);

	float filterRangeMin_dec, filterRangeMax_dec;
	filterRangeMin_dec = min / 180.0;
	filterRangeMax_dec = max / 180.0;

	float indicatoroffset;
	indicatoroffset = 0.1;
	this->outer2ArcPoints->InsertPoint(0, 0.0, 0.1, 0.0);
	this->outer2ArcPoints->InsertPoint(1, 0.0, -0.1, 0.0);
	for(i=0; i < arcparts1+1; i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.07 + indicatoroffset)*sin( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v2 = (this->axis_length-0.07 + indicatoroffset)*cos( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v3 = (this->axis_length-0.03 + indicatoroffset)*sin( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v4 = (this->axis_length-0.03 + indicatoroffset)*cos( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));

		this->outer2ArcPoints->InsertPoint(2+(i*2), v2, 0, v1);
		this->outer2ArcPoints->InsertPoint(2+(i*2)+1, v4, 0, v3);
    }

	for(i=arcparts1+1; i < (arcparts1+1 + arcparts2+1); i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.1 + indicatoroffset)*sin( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v2 = (this->axis_length-0.1 + indicatoroffset)*cos( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v3 = (this->axis_length + indicatoroffset)*sin( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v4 = (this->axis_length + indicatoroffset)*cos( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));

		this->outer2ArcPoints->InsertPoint(2+(i*2), v2, 0, v1);
		this->outer2ArcPoints->InsertPoint(2+(i*2)+1, v4, 0, v3);
    }

	for(i=arcparts1+1+arcparts2+1; i < (arcparts1+1 + arcparts2+1 + arcparts3+1); i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.07 + indicatoroffset)*sin( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v2 = (this->axis_length-0.07 + indicatoroffset)*cos( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v3 = (this->axis_length-0.03 + indicatoroffset)*sin( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v4 = (this->axis_length-0.03 + indicatoroffset)*cos( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));

		this->outer2ArcPoints->InsertPoint(2+(i*2), v2, 0, v1);
		this->outer2ArcPoints->InsertPoint(2+(i*2)+1, v4, 0, v3);
    }



	////printf("%f %f %f\n", this->minArcPoint[0], this->minArcPoint[1], this->minArcPoint[2]);
	////printf("%f %f %f\n", this->maxArcPoint[0], this->maxArcPoint[1], this->maxArcPoint[2]);

	vtkCellArray *cells2A = vtkCellArray::New();	
	vtkCellArray *cells2B = vtkCellArray::New();	
	j=1;
	for(i=0; i<arcparts1; i++)
	{	
		cells2A->InsertNextCell(4);		
		cells2A->InsertCellPoint(j*2);			
		cells2A->InsertCellPoint(j*2+1);
		cells2A->InsertCellPoint(j*2+3);
		cells2A->InsertCellPoint(j*2+2);
		j=j+1;
	}
	j=j+1;

	for(i=0; i<arcparts2; i++)
	{	
		cells2B->InsertNextCell(4);		
		cells2B->InsertCellPoint(j*2);			
		cells2B->InsertCellPoint(j*2+1);
		cells2B->InsertCellPoint(j*2+3);
		cells2B->InsertCellPoint(j*2+2);
		j=j+1;
	}
	j=j+1;

	for(i=0; i<arcparts3; i++)
	{	
		cells2A->InsertNextCell(4);		
		cells2A->InsertCellPoint(j*2);			
		cells2A->InsertCellPoint(j*2+1);
		cells2A->InsertCellPoint(j*2+3);
		cells2A->InsertCellPoint(j*2+2);
		j=j+1;
	}

	//printf("%i", j-1);


	this->outer2boogje1 = vtkPolyData::New();
	this->outer2boogje1->SetPolys(cells2A);	
	this->outer2boogje1->SetPoints(this->outer2ArcPoints);
	cells2A->Delete();

	this->outer2boogje2 = vtkPolyData::New();
	this->outer2boogje2->SetPolys(cells2B);	
	this->outer2boogje2->SetPoints(this->outer2ArcPoints);
	cells2B->Delete();
	
	this->outer2ArcMapA = vtkDataSetMapper::New();
	this->outer2ArcMapA->SetInput(this->outer2boogje1);
	this->outer2ArcA->SetMapper(this->outer2ArcMapA);
	this->outer2ArcA->SetProperty(this->HandlePropertyDarkYellow);
	this->outer2ArcMapB = vtkDataSetMapper::New();
	this->outer2ArcMapB->SetInput(this->outer2boogje2);
	this->outer2ArcB->SetMapper(this->outer2ArcMapB);
	this->outer2ArcB->SetProperty(this->HandlePropertyYellow);
}


void vtkAdvancedAngleWidget::DestroyArcRange2()
{
	this->outer2ArcPoints->Delete();
    this->outer2boogje1->Delete();
	this->outer2boogje2->Delete();
	this->outer2ArcMapA->Delete();
	this->outer2ArcMapB->Delete();
}





void vtkAdvancedAngleWidget::SetWidgetRange1(double min, double max)
{
	if (this->filterRangeMin == this->absoluteRangeMin1)
		this->filterRangeMin = min;
	if (this->filterRangeMax == this->absoluteRangeMax1)
		this->filterRangeMax = max;
	this->absoluteRangeMin1 = min;
	this->absoluteRangeMax1 = max;
}

void vtkAdvancedAngleWidget::SetWidgetRange2(double min, double max)
{
	if (this->filterRangeMin == this->absoluteRangeMin2)
		this->filterRangeMin = min;
	if (this->filterRangeMax == this->absoluteRangeMax2)
		this->filterRangeMax = max;
	this->absoluteRangeMin2 = min;
	this->absoluteRangeMax2 = max;
}

void vtkAdvancedAngleWidget::UpdateRangeIndicators()
{
	this->DestroyArcRange1();
	this->CreateArcRange1(this->absoluteRangeMin1, this->absoluteRangeMax1);
	this->DestroyArcRange2();
	this->CreateArcRange2(this->absoluteRangeMin2, this->absoluteRangeMax2);
	this->DestroyPrimaryArc();
	this->CreatePrimaryArc(this->filterRangeMin, this->filterRangeMax);
	this->Interactor->Render();
}


void vtkAdvancedAngleWidget::CreatePrimaryArc(float min, float max)
{
	float Pi;
    float a1[3], a2[3];
	int i, j;
	double arcparts1, arcparts2, arcparts3;

	if (min < 0){
		arcparts1 = 0;
	}
	else
		arcparts1 = (int)((min)/2);
	arcparts2 = (int)((max - min)/2);
	if (max>180)
		arcparts3 = 0;
	else
		arcparts3 = (int)(90-(max)/2);
	Pi = 3.141593;

	a1[0]=0.0;a1[1]=0.0;a1[2]=1.0;
	a2[0]=1.0;a2[1]=0.0;a2[2]=0.0;
	this->innerArcPoints = vtkPoints::New();
	this->innerArcPoints->SetNumberOfPoints((arcparts1+arcparts2 + arcparts3 + 3)*2+1);

	float filterRangeMin_dec, filterRangeMax_dec;
	filterRangeMin_dec = min / 180.0;
	filterRangeMax_dec = max / 180.0;

	this->innerArcPoints->InsertPoint(0, 0.0, 0.0, 0.0);

	for(i=0; i < arcparts1+1; i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.12)*sin( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v2 = (this->axis_length-0.12)*cos( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v3 = (this->axis_length-0.3)*sin( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		v4 = (this->axis_length-0.3)*cos( 1.5*Pi - (i/((float)arcparts1)) * (Pi*(filterRangeMin_dec)));
		this->innerArcPoints->InsertPoint(1+(i*2), v2, 0, v1);
		this->innerArcPoints->InsertPoint(2+(i*2), v4, 0, v3);
    }



	for(i=arcparts1+1; i < (arcparts1+1 + arcparts2+1); i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.12)*sin( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v2 = (this->axis_length-0.12)*cos( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v3 = (this->axis_length-0.3)*sin( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		v4 = (this->axis_length-0.3)*cos( 1.5*Pi - ((i-(arcparts1+1))/((float)arcparts2)) * (Pi*(filterRangeMax_dec - filterRangeMin_dec)) - (Pi * filterRangeMin_dec));
		this->innerArcPoints->InsertPoint(1+(i*2), v2, 0, v1);
		this->innerArcPoints->InsertPoint(2+(i*2), v4, 0, v3);
    }

	for(i=arcparts1+1+arcparts2+1; i < (arcparts1+1 + arcparts2+1 + arcparts3+1); i++)
    {
		float v1,v2,v3,v4;
		v1 = (this->axis_length-0.12)*sin( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v2 = (this->axis_length-0.12)*cos( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v3 = (this->axis_length-0.3)*sin( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		v4 = (this->axis_length-0.3)*cos( 1.5*Pi - ((i-(arcparts1+1+arcparts2+1))/((float)arcparts3)) * (Pi*(1.0-filterRangeMax_dec)) - (Pi * (filterRangeMax_dec)));
		this->innerArcPoints->InsertPoint(1+(i*2), v2, 0, v1);
		this->innerArcPoints->InsertPoint(2+(i*2), v4, 0, v3);
    }

	vtkCellArray *cells2A = vtkCellArray::New();	
	vtkCellArray *cells2B = vtkCellArray::New();	

	if (this->piechart){
		j=0;
		for(i=0; i<arcparts1; i++)
		{	
			cells2B->InsertNextCell(3);		
			cells2B->InsertCellPoint(0);
			cells2B->InsertCellPoint(j*2+1);
			cells2B->InsertCellPoint(j*2+3);
			j=j+1;
		}
		j=j+1;

		for(i=0; i<arcparts2; i++)
		{	
			cells2A->InsertNextCell(3);		
			cells2A->InsertCellPoint(0);			
			cells2A->InsertCellPoint(j*2+1);
			cells2A->InsertCellPoint(j*2+3);
			j=j+1;
		}
		j=j+1;

		for(i=0; i<arcparts3; i++)
		{	
			cells2B->InsertNextCell(3);		
			cells2B->InsertCellPoint(0);			
			cells2B->InsertCellPoint(j*2+1);
			cells2B->InsertCellPoint(j*2+3);
			j=j+1;
		}
	}else{
		j=0;
		for(i=0; i<arcparts1; i++)
		{	
			cells2B->InsertNextCell(4);		
			cells2B->InsertCellPoint(j*2+1);
			cells2B->InsertCellPoint(j*2+3);
			cells2B->InsertCellPoint(j*2+4);
			cells2B->InsertCellPoint(j*2+2);
			j=j+1;
		}
		j=j+1;

		for(i=0; i<arcparts2; i++)
		{	
			cells2A->InsertNextCell(4);		
			cells2A->InsertCellPoint(j*2+1);
			cells2A->InsertCellPoint(j*2+3);
			cells2A->InsertCellPoint(j*2+4);
			cells2A->InsertCellPoint(j*2+2);
			j=j+1;
		}
		j=j+1;

		for(i=0; i<arcparts3; i++)
		{	
			cells2B->InsertNextCell(4);		
			cells2B->InsertCellPoint(j*2+1);
			cells2B->InsertCellPoint(j*2+3);
			cells2B->InsertCellPoint(j*2+4);
			cells2B->InsertCellPoint(j*2+2);
			j=j+1;
		}
	}

	//printf("%i", j-1);


	this->innerboogje1 = vtkPolyData::New();
	this->innerboogje1->SetPolys(cells2A);	
	this->innerboogje1->SetPoints(this->innerArcPoints);
	cells2A->Delete();

	this->innerboogje2 = vtkPolyData::New();
	this->innerboogje2->SetPolys(cells2B);	
	this->innerboogje2->SetPoints(this->innerArcPoints);
	cells2B->Delete();
	
	this->innerArcMapA = vtkDataSetMapper::New();
	this->innerArcMapA->SetInput(this->innerboogje1);
	this->innerArcA->SetMapper(this->innerArcMapA);
	this->innerArcA->SetProperty(this->HandlePropertyOrange);
	this->innerArcMapB = vtkDataSetMapper::New();
	this->innerArcMapB->SetInput(this->innerboogje2);
	this->innerArcB->SetMapper(this->innerArcMapB);
	this->innerArcB->SetProperty(this->HandlePropertyDarkBlue);




	if (this->LabelsEnabled)
	{
		if (this->State == vtkAdvancedAngleWidget::Moving)
		{
			double oldpos[3];
			double newpos[3];
			char buf[10];
			double scalevalue;
			scalevalue = this->Scale/240.0;
			vtkTransform *genTransf = vtkTransform::New();
			if (this->upperlowerbool==0 || this->upperlowerbool==2){	
				oldpos[0] = (this->axis_length +0.2)*cos( 1.5*Pi - (Pi * filterRangeMin_dec));
				oldpos[2] = (this->axis_length +0.2)*sin( 1.5*Pi - (Pi * filterRangeMin_dec)) - 0.05;
				oldpos[1] = 0.0;
				genTransf->SetMatrix(this->innerArcA->GetUserMatrix());
				genTransf->TransformPoint(oldpos, newpos);
				
				this->textActorMin->SetPosition(0,0,0);
				this->textActorMin->SetOrientation(0,0,0);
				vtkTransform *temptransform = vtkTransform::New();
				temptransform->SetMatrix(this->innerArcA->GetUserMatrix());
				temptransform->PreMultiply();
				temptransform->RotateX(90);								
				this->textActorMin->SetScale(-scalevalue * mirrorfonts, scalevalue, scalevalue);
				this->textActorMin->SetOrientation(temptransform->GetOrientation());
				this->textActorMin->SetPosition(newpos);

				// itoa(this->filterRangeMin, buf, 10);
                // standard alternative:
                snprintf(buf, 10, "%d", int(this->filterRangeMin));
				this->textActorMin->SetInput(buf);
			}
			if (this->upperlowerbool==1 || this->upperlowerbool==2){	
				oldpos[0] = (this->axis_length +0.2)*cos( 1.5*Pi - (Pi * filterRangeMax_dec));
				oldpos[2] = (this->axis_length +0.2)*sin( 1.5*Pi - (Pi * filterRangeMax_dec)) - 0.05;
				//oldpos[0] = this->maxArcPoint[0];// - 0.2*this->Scale;
				oldpos[1] = 0.0;
				//oldpos[2] = this->maxArcPoint[2];// - 0.05*this->Scale;				
				genTransf->SetMatrix(this->innerArcA->GetUserMatrix());
				genTransf->TransformPoint(oldpos, newpos);
				
				this->textActorMax->SetPosition(0,0,0);
				this->textActorMax->SetOrientation(0,0,0);
				vtkTransform *temptransform = vtkTransform::New();
				temptransform->SetMatrix(this->innerArcA->GetUserMatrix());
				temptransform->PreMultiply();
				temptransform->RotateX(90);								
				this->textActorMax->SetScale(-scalevalue * mirrorfonts, scalevalue, scalevalue);
				this->textActorMax->SetOrientation(temptransform->GetOrientation());
				this->textActorMax->SetPosition(newpos);

				//itoa(this->filterRangeMax, buf, 10);
                // alternative:
                snprintf(buf, 10, "%d", int(this->filterRangeMax));
				this->textActorMax->SetInput(buf);
			}
		}			
	}
	
}

void vtkAdvancedAngleWidget::DestroyPrimaryArc()
{
	this->innerArcPoints->Delete();
    this->innerboogje1->Delete();
	this->innerboogje2->Delete();
	this->innerArcMapA->Delete();
	this->innerArcMapB->Delete();
}


void vtkAdvancedAngleWidget::CreateDefaultProperties()
{  
  //this->widgetcolor[0], this->widgetcolor[1], this->widgetcolor[2]);

  this->HandlePropertyBlue->SetColor(0.0, 0.0, 1.0);
  this->HandlePropertyBlue->SetAmbient(1);
  this->HandlePropertyBlue->SetDiffuse(0);
  this->HandlePropertyDarkBlue->SetColor(0.0, 0.0, 0.5);
  this->HandlePropertyDarkBlue->SetAmbient(1);
  this->HandlePropertyDarkBlue->SetDiffuse(0);
  this->HandlePropertyLightBlue->SetColor(0.5, 0.5, 1.0);
  this->HandlePropertyLightBlue->SetAmbient(1);
  this->HandlePropertyLightBlue->SetDiffuse(0);

  this->HandlePropertyYellow->SetColor(1.0, 1.0, 0.0);  
  this->HandlePropertyYellow->SetAmbient(1);
  this->HandlePropertyYellow->SetDiffuse(0);
  this->HandlePropertyDarkYellow->SetColor(0.5, 0.5, 0.0);
  this->HandlePropertyDarkYellow->SetAmbient(1);
  this->HandlePropertyDarkYellow->SetDiffuse(0);
  this->HandlePropertyLightYellow->SetColor(1.0, 0.2, 0.2);
  this->HandlePropertyLightYellow->SetAmbient(1);
  this->HandlePropertyLightYellow->SetDiffuse(0);

  this->HandlePropertyOrange->SetColor(1.0, 0.5, 0.0);  
  this->HandlePropertyOrange->SetAmbient(1);
  this->HandlePropertyOrange->SetDiffuse(0);

  this->HandlePropertyNudge->SetColor(1.0, 1.0, 1.0);  

  double widgetcolorbright[3];
  widgetcolorbright[0] = this->widgetcolor[0]+0.5;
  widgetcolorbright[1] = this->widgetcolor[1]+0.5;
  widgetcolorbright[2] = this->widgetcolor[2]+0.5;
  
  for (int i=0; i<3; i++)
	  if (widgetcolorbright[i]>1.0)
		  widgetcolorbright[i] = 1.0;
  

  this->HandlePropertyHighlight->SetColor(1,1,0);
  this->HandlePropertyHighlight->SetAmbient(1);
  this->HandlePropertyHighlight->SetDiffuse(0);
}

void vtkAdvancedAngleWidget::SetColor(double color[3])
{
	this->HandlePropertyNudge->SetColor(color);  
}


void vtkAdvancedAngleWidget::SetEnabled(int enabling)
{
	int k;

  if ( ! this->Interactor )
    {
    vtkErrorMacro(<<"The interactor must be set prior to enabling/disabling widget");
    return;
    }

  if (!this->FileName || (this->FileName && (0==strlen(this->FileName))))
    {
    vtkErrorMacro(<<"A FileName must be specified.");
    return;
    }


  if ( enabling ) //------------------------------------------------------------
    {
		vtkDebugMacro(<<"Enabling widget");

		if ( this->Enabled ) //already enabled, just return
		  {
		  return;
		  }
	    
		if ( ! this->CurrentRenderer )
		  {this->SetCurrentRenderer((vtkRenderer *)(this->Interactor->GetRenderWindow()->GetRenderers()->GetItemAsObject(1)));
		  /*this->SetCurrentRenderer(this->Interactor->FindPokedRenderer(
			this->Interactor->GetLastEventPosition()[0],
			this->Interactor->GetLastEventPosition()[1]));*/

		  //this->SetCurrentRenderer((vtkRenderer *)(this->CurrentRenderer->GetRenderWindow()->GetRenderers()->GetItemAsObject(1)));
		  if (this->CurrentRenderer == NULL)
			{
			return;
			}
		  }
		this->Enabled = 1;

		vtkCamera *cam = this->CurrentRenderer->GetActiveCamera();
		double *vp = this->CurrentRenderer->GetViewport();

		vtkRenderWindow* renwin = this->CurrentRenderer->GetRenderWindow();
		//renwin->AddRenderer( this->Renderer );


		/*if (renwin->GetNumberOfLayers() < 2)
		{
			renwin->SetNumberOfLayers( 2 );
		}*/

		// listen to the following events
		vtkRenderWindowInteractor *i = this->Interactor;
		i->AddObserver(vtkCommand::MouseMoveEvent, this->EventCallbackCommand, 
					   this->Priority);
		i->AddObserver(vtkCommand::LeftButtonPressEvent, 
					   this->EventCallbackCommand, this->Priority);
		i->AddObserver(vtkCommand::LeftButtonReleaseEvent, 
					   this->EventCallbackCommand, this->Priority);

		// Add the various actors
		// Add the outline

		this->CurrentRenderer->AddActor(this->outerArcA);
		this->CurrentRenderer->AddActor(this->outerArcB);
		this->CurrentRenderer->AddActor(this->outer2ArcA);
		this->CurrentRenderer->AddActor(this->outer2ArcB);
		this->CurrentRenderer->AddActor(this->innerArcA);
		if (this->showContext){
			this->CurrentRenderer->AddActor(this->innerArcB);
		}		
		this->CurrentRenderer->AddActor(this->nudgeactor);
		this->CurrentRenderer->AddActor(this->textActorMin);
		this->CurrentRenderer->AddActor(this->textActorMax);


		//this->Renderer->SetActiveCamera(cam);

		this->InvokeEvent(vtkCommand::EnableEvent,NULL);
    }

  else //disabling-------------------------------------------------------------
    {
		vtkDebugMacro(<<"Disabling widget");
		if ( ! this->Enabled ) //already disabled, just return
		  {
		  return;
		  }
		this->Enabled = 0;
		// don't listen for events any more
		this->Interactor->RemoveObserver(this->EventCallbackCommand);
		// turn off the handles  //other layer: Renderer instead of CurrentRenderer
		this->CurrentRenderer->RemoveActor(this->innerArcA);
		this->CurrentRenderer->RemoveActor(this->innerArcB);
		this->CurrentRenderer->RemoveActor(this->outerArcA);
		this->CurrentRenderer->RemoveActor(this->outerArcB);
		this->CurrentRenderer->RemoveActor(this->outer2ArcA);
		this->CurrentRenderer->RemoveActor(this->outer2ArcB);
		this->CurrentRenderer->RemoveActor(this->nudgeactor);
		this->CurrentHandle = NULL;
		this->InvokeEvent(vtkCommand::DisableEvent,NULL);
		this->SetCurrentRenderer(NULL);
    }

  this->Interactor->Render();
}

void vtkAdvancedAngleWidget::ProcessEvents(vtkObject* vtkNotUsed(object), 
                                 unsigned long event,
                                 void* clientdata, 
                                 void* vtkNotUsed(calldata))
{
  vtkAdvancedAngleWidget* self = reinterpret_cast<vtkAdvancedAngleWidget *>( clientdata );

  //okay, let's do the right thing
  switch(event)
    {
    case vtkCommand::LeftButtonPressEvent:
      self->OnLeftButtonDown();
      break;
    case vtkCommand::LeftButtonReleaseEvent:
      self->OnLeftButtonUp();
      break;
    case vtkCommand::MouseMoveEvent:
      self->OnMouseMove();
      break;
    }
}

#define VTK_AVERAGE(a,b,c) \
  c = (a + b)/2.0;

#define VTK_AVSIZE(a,b,c) \
  c = (a - b);\
  if (c<0) c=-c;

void vtkAdvancedAngleWidget::PlaceWidget(double bds[6])
{
	double bounds[6], centertemp[3];  
	double s1,s2,s3, p;
	this->AdjustBounds(bds,bounds,centertemp);
	VTK_AVERAGE(bds[0],bds[1],this->center[0])
	VTK_AVERAGE(bds[2],bds[3],this->center[1])
	VTK_AVERAGE(bds[4],bds[5],this->center[2])

	VTK_AVSIZE(bds[0],bds[1],s1);
	VTK_AVSIZE(bds[2],bds[3],s2);
	VTK_AVSIZE(bds[4],bds[5],s3);
	p = ((s1 + s2 + s3)/3);   
	this->WidgetScale = p;
	//reset the translation and rotation data
	this->pos[0]=0;
	this->pos[1]=0;
	this->pos[2]=0;
	this->orientation->Identity();
	this->totalangle=0;

	this->CreateNudge();
	this->CreatePrimaryArc(this->filterRangeMin, this->filterRangeMax);
	this->CreateArcRange1(this->absoluteRangeMin1, this->absoluteRangeMax1);	
	this->CreateArcRange2(this->absoluteRangeMin2, this->absoluteRangeMax2);	

	this->Realign();
}

void vtkAdvancedAngleWidget::CreateNudge()
{
	this->stlreader = vtkSTLReader::New();
	this->stlreader->SetFileName(this->FileName);
	this->stlreader->Update();
	
	this->nudgemapper = vtkDataSetMapper::New();
	this->nudgemapper->SetInput(this->stlreader->GetOutput());
	this->nudgeactor->SetMapper(this->nudgemapper);
	this->nudgeactor->SetProperty(this->HandlePropertyNudge);
}

void vtkAdvancedAngleWidget::SetInitialPosition(double position[3])
{
	this->center[0] = position[0] + this->posoffset[0];
	this->center[1] = position[1] + this->posoffset[1];
	this->center[2] = position[2] + this->posoffset[2];
	this->Realign();
}

void vtkAdvancedAngleWidget::SetInitialAxis(double primaryAxis[3])
{
	double beginas[3];
	double doelas[3];
	double N1[3];
	double dottie;
	double a;

	beginas[0] = 0.0;
	beginas[1] = 1.0;
	beginas[2] = 0.0;
	this->widgetaxis[0] = primaryAxis[0];
	this->widgetaxis[1] = primaryAxis[1];
	this->widgetaxis[2] = primaryAxis[2];

	a = vtkMath::Normalize(this->widgetaxis);
	vtkMath::Cross(beginas, this->widgetaxis, N1);     
	dottie = acos(vtkMath::Dot(beginas,this->widgetaxis) / (sqrt(pow(beginas[0],2) + pow(beginas[1],2) + pow(beginas[2],2))*sqrt(pow(this->widgetaxis[0],2) + pow(this->widgetaxis[1],2) + pow(this->widgetaxis[2],2))));
	dottie = vtkMath::DegreesFromRadians(dottie);
	//printf("orientatie-as verschil RotationArcWidget %f %f %f", N1[0], N1[1], N1[2]);
	//printf("orientatie verschil RotationArcWidget %f", dottie);

	this->initialaxis[0] = N1[1]; //raar, maar das iets met hoe de widget geplaatst wordt t.o.v. de assen :S
	this->initialaxis[1] = N1[2]; //Het gaat overigens over het verschil met de initiele as.
	this->initialaxis[2] = N1[0];

	this->widgetorientation->Identity();	
	this->widgetorientation->RotateWXYZ(dottie,N1[0],N1[1],N1[2]);
	this->Realign();
}



#undef VTK_AVSIZE

void vtkAdvancedAngleWidget::Realign()
{
  vtkTransform *t = vtkTransform::New();
  this->GetTransform(t);
  vtkMatrix4x4 *temp = t->GetMatrix();
  this->innerArcA->SetUserMatrix(temp);
  this->innerArcB->SetUserMatrix(temp);
  this->outerArcA->SetUserMatrix(temp);
  this->outerArcB->SetUserMatrix(temp);
  this->outer2ArcA->SetUserMatrix(temp);
  this->outer2ArcB->SetUserMatrix(temp);
  this->nudgeactor->SetUserMatrix(temp);
}

#undef VTK_AVERAGE

void vtkAdvancedAngleWidget::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os,indent);
}


void vtkAdvancedAngleWidget::ToggleSize()
{
	if (this->axis_length == 0.3){
		//this->axis_length = this->initial_axis_length;
		//this->nudgeactor->SetScale(1.0, 1.0, 1.0);

		//this->initial_axis_length = this->axis_length;
		float increment = (this->initial_axis_length - this->axis_length)/5.0;
		for (int i=0; i<5; i++){
			this->center[0] = this->center[0] + this->posoffset[0]/7.0;
			this->center[1] = this->center[1] + this->posoffset[1]/7.0;
			this->center[2] = this->center[2] + this->posoffset[2]/7.0;
			/*this->posoffset[0] = this->posoffset[0] * 1.1764705882352941176470588235294;
			this->posoffset[1] = this->posoffset[1] * 1.1764705882352941176470588235294;
			this->posoffset[2] = this->posoffset[2] * 1.1764705882352941176470588235294;*/
			this->axis_length = this->axis_length + increment;	
			this->UpdateRangeIndicators();
			this->GetCurrentRenderer()->GetRenderWindow()->Render();
			this->nudgeactor->SetScale(0.5+(i+1)*0.1, 1.0, 0.5+(i+1)*0.1);
			this->Realign();
		}
	}
	else
	{
		this->initial_axis_length = this->axis_length;
		float increment = (this->initial_axis_length - 0.3)/5.0;
		for (int i=0; i<5; i++){
			this->center[0] = this->center[0] - this->posoffset[0]/7.0;
			this->center[1] = this->center[1] - this->posoffset[1]/7.0;
			this->center[2] = this->center[2] - this->posoffset[2]/7.0;
			this->axis_length = this->axis_length - increment;	
			this->UpdateRangeIndicators();
			this->GetCurrentRenderer()->GetRenderWindow()->Render();
			this->nudgeactor->SetScale(1.0-(i+1)*0.1, 1.0, 1.0-(i+1)*0.1);
			this->Realign();
		}
		this->axis_length = 0.3;		
	}
	//this->UpdateRangeIndicators();
	
}


void vtkAdvancedAngleWidget::OnLeftButtonDown()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  // Okay, we can process this. Try to pick handles first;
  // if no handles picked, then pick the bounding box.
  if (!this->CurrentRenderer || !this->CurrentRenderer->IsInViewport(X, Y))
    {
    this->State = vtkAdvancedAngleWidget::Outside;
    return;
    }
  
  vtkAssemblyPath *path;
  this->HandlePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->HandlePicker->GetPath();
  if ( path != NULL )
		{
			this->State = vtkAdvancedAngleWidget::Moving;
			this->HighlightHandle(path->GetFirstNode()->GetViewProp());
			this->HandlePicker->GetPickPosition(this->LastPickPosition);

			double transformed_minArcPoint[3];
			double transformed_maxArcPoint[3];
			vtkTransform *genTransf = vtkTransform::New();
			genTransf->SetMatrix(this->innerArcA->GetUserMatrix());
			genTransf->TransformPoint(this->minArcPoint, transformed_minArcPoint);
			genTransf->TransformPoint(this->maxArcPoint, transformed_maxArcPoint);
			
			//get pick position
			//determine absolute angle between center-pick and center-midvector
			//see if this is larger than crossangle1
			//see if this is smaller than crossangle2
			double dottie, dottieangle, controlehoek, pickedangle, minangle, maxangle, minor, major;
			double v1[3], v2[3], tmpvec1[3], transformed_tmpvec1[3], zaxis[3], yaxis[3];
			double tmpvec2[3], transformed_tmpvec2[3], y_transformed[3], controlevector[3];

			zaxis[0] = 0.0;
			zaxis[1] = 0.0;
			zaxis[2] = -1.0;
			yaxis[0] = 0.0;
			yaxis[1] = 1.0;
			yaxis[2] = 0.0;
			

			//printf("v1: %f %f %f\n", *v1, *(v1+1), *(v1+2));

			
			//v1 = this->LastPickPosition - this->center;
			this->SubtractVectors(this->LastPickPosition, this->center, v1);
			vtkMath::Normalize(v1);
			//printf("v1: %f %f %f\n", *v1, *(v1+1), *(v1+2));

			//v2 = transform(this->center+[1,0,0]) - this->center;
			//this->AddVectors(this->center, xaxis, tmpvec1);
			genTransf->TransformPoint(zaxis, transformed_tmpvec1);
			this->SubtractVectors(transformed_tmpvec1, this->center, v2);
			vtkMath::Normalize(v2);
			//printf("v2: %f %f %f\n", *v2, *(v2+1), *(v2+2));

			dottie = vtkMath::Dot(v1, v2); 
			//printf("dottie: %f\n", dottie);
			//controlehoek = vtkMath::Round( vtkMath::Dot(v1, transform(this->center+[0,0,1]) - this->center) );
			//this->AddVectors(this->center, yaxis, tmpvec2);
			//genTransf->TransformPoint(tmpvec2, transformed_tmpvec2);
			
			genTransf->TransformPoint(yaxis, transformed_tmpvec2);
			this->SubtractVectors(transformed_tmpvec2, this->center, y_transformed);
			vtkMath::Normalize(y_transformed);
			vtkMath::Cross(v1, v2, controlevector);
			vtkMath::Normalize(controlevector);
			//printf("controlevector: %f %f %f\n", *controlevector, *(controlevector+1), *(controlevector+2));
			//printf("control: %f\n", vtkMath::Dot(controlevector, y_transformed));

			dottieangle = acos(dottie / (sqrt(pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2))*sqrt(pow(v2[0],2) + pow(v2[1],2) + pow(v2[2],2))));
			pickedangle = -vtkMath::DegreesFromRadians(dottieangle) * vtkMath::Dot(controlevector, y_transformed);

			if (pickedangle < this->absoluteRangeMin1 && pickedangle < this->absoluteRangeMin2)
				pickedangle = 180.0 + (180.0+pickedangle);
			if (pickedangle > this->absoluteRangeMax1 && pickedangle > this->absoluteRangeMax2)
				pickedangle = -180.0 - (180.0 -pickedangle);
			
			/*printf("angle: %f\n", pickedangle);
			printf("this->filterRangeMin %f\n", this->filterRangeMin);
			printf("this->filterRangeMax: %f\n", this->filterRangeMax);*/

			
			minor = this->filterRangeMin + (this->filterRangeMax - this->filterRangeMin)/3.0;
			major = this->filterRangeMin + 2.0*(this->filterRangeMax - this->filterRangeMin)/3.0;

			if (pickedangle < minor)
				this->upperlowerbool = 0;
			else if (pickedangle < major)
				this->upperlowerbool = 2;
			else
				this->upperlowerbool = 1;


			this->ValidPick = 1;
		}
  else
		{
			this->HighlightHandle(NULL);
			this->State = vtkAdvancedAngleWidget::Outside;
			return;
		}

  this->EventCallbackCommand->SetAbortFlag(1);
  this->StartInteraction();
  this->InvokeEvent(vtkCommand::StartInteractionEvent, NULL);
  this->Interactor->Render();
}

void vtkAdvancedAngleWidget::SubtractVectors(double vector1[3], double vector2[3], double *resultvector)
{
	for (int i=0; i<3; i++)
		resultvector[i] = vector1[i] - vector2[i];
}

void vtkAdvancedAngleWidget::AddVectors(double vector1[3], double vector2[3], double *resultvector)
{
	for (int i=0; i<3; i++)
		resultvector[i] = vector1[i] + vector2[i];
}

void vtkAdvancedAngleWidget::OnLeftButtonUp()
{
  int X = this->Interactor->GetEventPosition()[0];
  int Y = this->Interactor->GetEventPosition()[1];

  if ( this->State == vtkAdvancedAngleWidget::Outside ||
       this->State == vtkAdvancedAngleWidget::Start )
    {
    return;
    }

  this->State = vtkAdvancedAngleWidget::Start;
  this->HighlightHandle(NULL);
  this->textActorMin->SetInput("");
  this->textActorMax->SetInput("");

  vtkAssemblyPath *path;
  this->HandlePicker->Pick(X,Y,0.0,this->CurrentRenderer);
  path = this->HandlePicker->GetPath();
  if ( path != NULL )
  {
		if (path->GetFirstNode()->GetViewProp() == this->nudgeactor )
		{
			this->ToggleSize();
		}
  }

  this->EventCallbackCommand->SetAbortFlag(1);
  this->EndInteraction();
  this->InvokeEvent(vtkCommand::EndInteractionEvent, NULL);
  this->Interactor->Render(); 
}

void vtkAdvancedAngleWidget::OnMouseMove()
{
  int i,j=0;

	// See whether we're active
  if ( this->State == vtkAdvancedAngleWidget::Outside || 
       this->State == vtkAdvancedAngleWidget::Start )
    {
    return;
    }

  double X = (double)this->Interactor->GetEventPosition()[0];
  double Y = (double)this->Interactor->GetEventPosition()[1];

  // Do different things depending on state
  // Calculations everybody does
  double focalPoint[4], pickPointOld[4], prevPickPointOld[4], prev[4], pick[4];
  double u, *as, beginas[4], temp[4], temp2[4], temppoint[4], temppoint2[4], temppoint3[3], temppoint4[3], temppoint5, temppoint6, point_on_pickline1[3], point_on_pickline2[3], point_on_pickingplane1[4], point_on_pickingplane2[4];
  double *pickPoint,*prevPickPoint;
  double z; 
  double vpn[3];
  
  // Process the motion
  if ( this->State == vtkAdvancedAngleWidget::Moving )
    {
    // Okay to process
    if ( this->CurrentHandle )
      {
  		if (this->CurrentHandle == this->innerArcA)
			{beginas[0]=0.0; beginas[1]=1.0; beginas[2]=0.0;}
      }
	}

  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  if ( !camera )
    {
    return;
    }

  camera->GetViewPlaneNormal(vpn);

  int *windowsize;
  double aspect;
  double *clippingplanes;
  windowsize = this->CurrentRenderer->GetRenderWindow()->GetSize();
  aspect = windowsize[0]/windowsize[1];
  clippingplanes = camera->GetClippingRange();
	  
  // Compute the two points defining the motion vector
  this->ComputeWorldToDisplay(this->LastPickPosition[0], this->LastPickPosition[1],
                              this->LastPickPosition[2], focalPoint);  
  z = focalPoint[2];

  if (this->CurrentHandle == this->innerArcA)
  {
	this->ComputeDisplayToWorld(double(this->Interactor->GetLastEventPosition()[0]),
							double(this->Interactor->GetLastEventPosition()[1]),
							z, temppoint);

	vtkTransform *t = vtkTransform::New();
	this->widgetorientation->GetOrientationWXYZ(temp);
	t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
	this->orientation->GetOrientationWXYZ(temp2);
	t->RotateWXYZ(temp2[0],temp2[1],temp2[2],temp2[3]);
	as = t->TransformDoublePoint(beginas);

	//printf("AAAAAAAAAAAAAS %f %f %f\n", *as, *(as+1), *(as+2));
	//printf("temp %f %f %f %f\n", *temp, *(temp+1), *(temp+2), *(temp+3));
	//printf("temp2 %f %f %f %f\n", *temp2, *(temp2+1), *(temp2+2), *(temp2+3));

	point_on_pickline1[0] = temppoint[0];
	point_on_pickline1[1] = temppoint[1];
	point_on_pickline1[2] = temppoint[2];
	point_on_pickline2[0] = camera->GetPosition()[0];
	point_on_pickline2[1] = camera->GetPosition()[1];
	point_on_pickline2[2] = camera->GetPosition()[2];  
	temppoint3[0] = this->pos[0] + this->center[0] - point_on_pickline2[0];
	temppoint3[1] = this->pos[1] + this->center[1] - point_on_pickline2[1];
	temppoint3[2] = this->pos[2] + this->center[2] - point_on_pickline2[2];
	temppoint4[0] = point_on_pickline1[0] - point_on_pickline2[0];
	temppoint4[1] = point_on_pickline1[1] - point_on_pickline2[1];
	temppoint4[2] = point_on_pickline1[2] - point_on_pickline2[2];
	temppoint5 = vtkMath::Dot(as, temppoint3);
	temppoint6 = vtkMath::Dot(as, temppoint4);
	u = temppoint5 / temppoint6;
	temppoint2[0] = u*(point_on_pickline2[0] - point_on_pickline1[0]);
	temppoint2[1] = u*(point_on_pickline2[1] - point_on_pickline1[1]);
	temppoint2[2] = u*(point_on_pickline2[2] - point_on_pickline1[2]);  
	point_on_pickingplane1[0] = point_on_pickline2[0] - temppoint2[0];
	point_on_pickingplane1[1] = point_on_pickline2[1] - temppoint2[1];
	point_on_pickingplane1[2] = point_on_pickline2[2] - temppoint2[2];

	this->ComputeDisplayToWorld(double(X), double(Y), z, temppoint);
	point_on_pickline1[0] = temppoint[0];
	point_on_pickline1[1] = temppoint[1];
	point_on_pickline1[2] = temppoint[2];
	point_on_pickline2[0] = camera->GetPosition()[0];
	point_on_pickline2[1] = camera->GetPosition()[1];
	point_on_pickline2[2] = camera->GetPosition()[2];  
	temppoint3[0] = this->pos[0] + this->center[0] - point_on_pickline2[0];
	temppoint3[1] = this->pos[1] + this->center[1] - point_on_pickline2[1];
	temppoint3[2] = this->pos[2] + this->center[2] - point_on_pickline2[2];
	temppoint4[0] = point_on_pickline1[0] - point_on_pickline2[0];
	temppoint4[1] = point_on_pickline1[1] - point_on_pickline2[1];
	temppoint4[2] = point_on_pickline1[2] - point_on_pickline2[2];
	temppoint5 = vtkMath::Dot(as, temppoint3);
	temppoint6 = vtkMath::Dot(as, temppoint4);
	u = temppoint5 / temppoint6;
	temppoint2[0] = u*(point_on_pickline2[0] - point_on_pickline1[0]);
	temppoint2[1] = u*(point_on_pickline2[1] - point_on_pickline1[1]);
	temppoint2[2] = u*(point_on_pickline2[2] - point_on_pickline1[2]);   
	point_on_pickingplane2[0] = point_on_pickline2[0] - temppoint2[0];
	point_on_pickingplane2[1] = point_on_pickline2[1] - temppoint2[1];
	point_on_pickingplane2[2] = point_on_pickline2[2] - temppoint2[2];
    
	


	/*vtkSphereSource *spheertje = vtkSphereSource::New();
	spheertje->SetRadius(0.05);
	vtkPolyDataMapper *spheertjemap = vtkPolyDataMapper::New();
	spheertjemap->SetInput(spheertje->GetOutput());
	vtkActor *spheertjeact = vtkActor::New();
	spheertjeact->SetMapper(spheertjemap);
	spheertjeact->SetPosition(point_on_pickingplane1[0], point_on_pickingplane1[1], point_on_pickingplane1[2]);
	this->CurrentRenderer->AddActor(spheertjeact);*/
  }

  double added_angle;
  

  // Process the motion
  if ( this->State == vtkAdvancedAngleWidget::Moving )
    {
    // Okay to process
    if ( this->CurrentHandle )
      {
        if (this->CurrentHandle == this->innerArcA)
		  {
		    added_angle = this->AxisRotate(2, point_on_pickingplane2, point_on_pickingplane1);
			if (this->upperlowerbool==1){
				this->filterRangeMax = this->filterRangeMax + added_angle;
				if (this->filterRangeMax > this->absoluteRangeMax1 && this->filterRangeMax > this->absoluteRangeMax2)
					if (this->absoluteRangeMax1 > this->absoluteRangeMax2)
						this->filterRangeMax = this->absoluteRangeMax1;
					else
						this->filterRangeMax = this->absoluteRangeMax2;
				if (this->filterRangeMax < this->filterRangeMin+2)
					this->filterRangeMax = this->filterRangeMin+2;
				}
			else if (this->upperlowerbool==0){
				this->filterRangeMin = this->filterRangeMin + added_angle;
				if (this->filterRangeMin > this->filterRangeMax-2)
					this->filterRangeMin = this->filterRangeMax-2;			
				if (this->filterRangeMin < this->absoluteRangeMin1 && this->filterRangeMin < this->absoluteRangeMin2)
					if (this->absoluteRangeMin1 < this->absoluteRangeMin2)
						this->filterRangeMin = this->absoluteRangeMin1;
					else
						this->filterRangeMin = this->absoluteRangeMin2;
				}
			else
			{
				if (this->absoluteRangeMax1 > this->absoluteRangeMax2){
					if (this->filterRangeMax + added_angle > this->absoluteRangeMax1)
						added_angle = added_angle-((this->filterRangeMax + added_angle) - this->absoluteRangeMax1);
				}
				else{
					if (this->filterRangeMax + added_angle > this->absoluteRangeMax2)
						added_angle = added_angle-((this->filterRangeMax + added_angle) - this->absoluteRangeMax2);
				}
				if (this->absoluteRangeMin1 < this->absoluteRangeMin2){
					if (this->filterRangeMin + added_angle < this->absoluteRangeMin1)
						added_angle = added_angle-((this->filterRangeMin + added_angle) - this->absoluteRangeMin1);					
				}
				else{
					if (this->filterRangeMin + added_angle < this->absoluteRangeMin2)
						added_angle = added_angle-((this->filterRangeMin + added_angle) - this->absoluteRangeMin2);					
				}
				this->filterRangeMin = this->filterRangeMin + added_angle;
				this->filterRangeMax = this->filterRangeMax + added_angle;
			}
			
		  	
  			this->DestroyPrimaryArc();
			this->CreatePrimaryArc(this->filterRangeMin, this->filterRangeMax);
			this->CurrentRenderer->AddActor(this->innerArcA);
			if (this->showContext){
				this->CurrentRenderer->AddActor(this->innerArcB);
			}

			this->HighlightHandle(this->innerArcA);//path->GetFirstNode()->GetViewProp());
			//this->Interactor->Render();
		  }
		else
		  {
		  this->AxisRotate(2, point_on_pickingplane2, point_on_pickingplane1);
          }
      }
	}

  // Interact, if desired
  this->EventCallbackCommand->SetAbortFlag(1);
  //this->StartInteraction();
  this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
  this->Interactor->Render();
}

int vtkAdvancedAngleWidget::HighlightHandle(vtkProp *prop)
{
  // first unhighlight anything pickable
  if ( this->CurrentHandle )
    {
		this->innerArcA->SetProperty(this->HandlePropertyOrange);
		this->nudgeactor->SetProperty(this->HandlePropertyNudge);
		//this->innerArcB->SetProperty(this->HandlePropertyDarkBlue);
	}

  this->CurrentHandle = (vtkActor *)prop;

  if ( this->CurrentHandle )
    {
		if ( this->CurrentHandle == this->innerArcA ){
			this->CurrentHandle->SetProperty(this->HandlePropertyHighlight);
		}
		else if ( this->CurrentHandle == this->nudgeactor )
		{
			this->CurrentHandle->SetProperty(this->HandlePropertyHighlight);
		}
    
	

	/*if(this->IsFlat())
	{
		this->CurrentHandle->GetProperty()->SetAmbient(1);
		this->CurrentHandle->GetProperty()->SetDiffuse(0);
	}
    */
	if ( this->CurrentHandle == this->innerArcA )
		return 1;
    }  
  return -1;
}

bool vtkAdvancedAngleWidget::IsFlat()
{
	if ( this->CurrentHandle == this->innerArcA )
		return true;
	else
		return false;
}

void vtkAdvancedAngleWidget::GetDirection(const double Nx[3],const double Ny[3], const double Nz[3], double dir[3])
{
  double dotNy, dotNz;
  double y[3];

  if(vtkMath::Dot(Nx,Nx)!=0)
    {
    dir[0] = Nx[0];
    dir[1] = Nx[1];
    dir[2] = Nx[2];
    }
  else 
    {
    dotNy = vtkMath::Dot(Ny,Ny);
    dotNz = vtkMath::Dot(Nz,Nz);
    if(dotNy != 0 && dotNz != 0)
      {
      vtkMath::Cross(Ny,Nz,dir);
      }
    else if(dotNy != 0)
      {
      //dir must have been initialized to the 
      //corresponding coordinate direction before calling
      //this method
      vtkMath::Cross(Ny,dir,y);
      vtkMath::Cross(y,Ny,dir);
      }
    else if(dotNz != 0)
      {
      //dir must have been initialized to the 
      //corresponding coordinate direction before calling
      //this method
      vtkMath::Cross(Nz,dir,y);
      vtkMath::Cross(y,Nz,dir);
      }
    }
}

void vtkAdvancedAngleWidget::AxisTranslate(int axis, double *p1, double *p2)
{
	int i;
	double v[3];
	double temp[4];
	double dir[3] = { 0,0,0};	
	double axisdirection[3];

	if (axis==0){
		axisdirection[0]=1;axisdirection[1]=0;axisdirection[2]=0;}
	else if (axis==1){
		axisdirection[0]=0;axisdirection[1]=0;axisdirection[2]=-1;}
	else {
		axisdirection[0]=0;axisdirection[1]=1;axisdirection[2]=0;}

	vtkTransform *t = vtkTransform::New();
	vtkTransform *t2 = vtkTransform::New();
	this->GetTransform(t);
	t->GetOrientationWXYZ(temp);
	t2->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);

	double *rotated_axisdirection;
	rotated_axisdirection = t2->TransformDoublePoint(axisdirection);

	double *ViewDirection;
	double *ViewPosition;
	ViewDirection = this->CurrentRenderer->GetActiveCamera()->GetDirectionOfProjection();
	ViewPosition = this->CurrentRenderer->GetActiveCamera()->GetPosition();

	double N1[3], N2[3];
	vtkMath::Cross(rotated_axisdirection, ViewDirection, N1);
	vtkMath::Normalize(N1);
	vtkMath::Cross(ViewDirection, N1, N2);
	vtkMath::Normalize(N2);

	double T, T2;
	double pointInPlane[3];
	double pointInPlane2[3];
	double temp7[3], temp8[3];
	double temp1, temp2, temp3,temp4, temp5, temp6;
	temp7[0] = p1[0]-pos[0];
	temp7[1] = p1[1]-pos[1];
	temp7[2] = p1[2]-pos[2];
	T = vtkMath::Dot(N2,temp7)/vtkMath::Dot(N2, rotated_axisdirection);
	temp1 = rotated_axisdirection[0]*(T);
	temp2 = rotated_axisdirection[1]*(T);
	temp3 = rotated_axisdirection[2]*(T);
	pointInPlane[0] = pos[0] + temp1;
	pointInPlane[1] = pos[1] + temp2;
	pointInPlane[2] = pos[2] + temp3;

	temp8[0] = p2[0]-pos[0];
	temp8[1] = p2[1]-pos[1];
	temp8[2] = p2[2]-pos[2];
	T2 = vtkMath::Dot(N2,temp8)/vtkMath::Dot(N2, rotated_axisdirection);
	temp4 = rotated_axisdirection[0]*(T2);
	temp5 = rotated_axisdirection[1]*(T2);
	temp6 = rotated_axisdirection[2]*(T2);
	pointInPlane2[0] = pos[0] + temp4;
	pointInPlane2[1] = pos[1] + temp5;
	pointInPlane2[2] = pos[2] + temp6;

	for (i=0; i<3; i++)
	{
		v[i] = pointInPlane2[i] - pointInPlane[i];
		this->pos[i]+=v[i];
	}

	this->Realign();
}

bool vtkAdvancedAngleWidget::IllegalPosition()
{
	if(this->pos[1]<-5 || this->pos[1]>5)
		return true;
	else if (sqrt(pow((this->pos[0]),2)+pow((this->pos[2]),2))>2.5)
		return true;
	else 
		return false;
}

double vtkAdvancedAngleWidget::AxisRotate(int axisconstraint, double *p1, double *p2)
{     
	double x_as[3], y_as[3], z_as[3], beginas[3], temp[4], *new_x_as, *new_y_as, *new_z_as, normal1[3], normal2[3], crossie[3], dottie;

	if (axisconstraint == 0)
		{beginas[0]=0;beginas[1]=0;beginas[2]=1;}
	else if (axisconstraint == 1)
		{beginas[0]=1;beginas[1]=0;beginas[2]=0;}
	else
		{beginas[0]=0;beginas[1]=1;beginas[2]=0;}
	
	vtkTransform *t = vtkTransform::New();
	this->widgetorientation->GetOrientationWXYZ(temp);
	t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
	this->orientation->GetOrientationWXYZ(temp);
	t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
	x_as[0]=1.0;x_as[1]=0.0;x_as[2]=0.0;
	y_as[0]=0.0;y_as[1]=1.0;y_as[2]=0.0;
	z_as[0]=0.0;z_as[1]=0.0;z_as[2]=1.0;

	new_x_as = t->TransformDoublePoint(x_as);
	x_as[0]=*new_x_as;x_as[1]=*(new_x_as+1);x_as[2]=*(new_x_as+2);
	new_y_as = t->TransformDoublePoint(y_as);
	y_as[0]=*new_y_as;y_as[1]=*(new_y_as+1);y_as[2]=*(new_y_as+2);
	new_z_as = t->TransformDoublePoint(z_as);
	z_as[0]=*new_z_as;z_as[1]=*(new_z_as+1);z_as[2]=*(new_z_as+2);

	normal1[0] = p1[0]-center[0]-pos[0];
	normal1[1] = p1[1]-center[1]-pos[1];
	normal1[2] = p1[2]-center[2]-pos[2];
	normal2[0] = p2[0]-center[0]-pos[0];
	normal2[1] = p2[1]-center[1]-pos[1];
	normal2[2] = p2[2]-center[2]-pos[2];
	vtkMath::Normalize(normal1);
	vtkMath::Normalize(normal2);

	//vtkSphereSource *spheertje = vtkSphereSource::New();
	//vtkPolyDataMapper *spheertjemap = vtkPolyDataMapper::New();
	//spheertjemap->SetInput(spheertje->GetOutput());
	//this->spheertjeact = vtkActor::New();
	//this->spheertjeact->SetMapper(spheertjemap);
	//this->spheertjeact->SetPosition(center[0],center[1],center[2]);
	//this->Renderer->AddActor(this->spheertjeact);

	vtkMath::Cross(normal1,normal2,crossie);
	dottie = acos(vtkMath::Dot(normal1,normal2) / (sqrt(pow(normal1[0],2) + pow(normal1[1],2) + pow(normal1[2],2))*sqrt(pow(normal2[0],2) + pow(normal2[1],2) + pow(normal2[2],2))));
	dottie = vtkMath::DegreesFromRadians(dottie);

	if (axisconstraint == 0)
		if ((z_as[0]*crossie[0] + z_as[1]*crossie[1] + z_as[2]*crossie[2])<0)
			dottie = abs(dottie)*-1;
		else
			dottie = abs(dottie);
	if (axisconstraint == 1)
		if ((x_as[0]*crossie[0] + x_as[1]*crossie[1] + x_as[2]*crossie[2])<0)
			dottie = abs(dottie)*-1;
		else
			dottie = abs(dottie);
	if (axisconstraint == 2)
		if ((y_as[0]*crossie[0] + y_as[1]*crossie[1] + y_as[2]*crossie[2])<0)
			dottie = abs(dottie);
		else
			dottie = abs(dottie)*-1;
			
	return dottie;
	
}

double vtkAdvancedAngleWidget::GetXOffset()
{
  return this->pos[0];
}
double vtkAdvancedAngleWidget::GetYOffset()
{
  return this->pos[2];
}
double vtkAdvancedAngleWidget::GetZOffset()
{
  return this->pos[1];
}

void vtkAdvancedAngleWidget::GetOriginalPosition(double point[3])
{
	point[0] = this->center[0];
	point[1] = this->center[1];
	point[2] = this->center[2];
	//printf("cplusplus zegt: %f, %f, %f\n", center[0], center[1], center[2]);
}

void vtkAdvancedAngleWidget::GetRotationTransform(vtkTransform *t)
{
	t->Identity();
	double temp[4];
	this->orientation->GetOrientationWXYZ(temp);
	t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);	
}

void vtkAdvancedAngleWidget::GetInitialOrientationTransform(vtkTransform *t)
{
	t->Identity();
	double temp[4];
	this->widgetorientation->GetOrientationWXYZ(temp);
	t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);	
}

void vtkAdvancedAngleWidget::Rotate(double angle, double axis[3])
{
	this->theta = angle;
	this->rotation_axis[0]=axis[0];
	this->rotation_axis[1]=axis[1];
	this->rotation_axis[2]=axis[2];
	this->Realign();
}

void vtkAdvancedAngleWidget::GetTransform(vtkTransform *t)
{
  double translate[3];
  double temp[4];

  this->widgetorientation->RotateWXYZ(this->theta,this->rotation_axis);
  this->totalangle = this->totalangle + this->theta;
  //printf("theta: %f\n", this->theta);
  //printf("totalangle: %f\n", this->totalangle);
  this->theta = 0;

  t->Identity();

    t->Translate(this->center[0],this->center[1],this->center[2]);
    t->Scale(this->Scale,this->Scale,this->Scale);
    t->Translate(-this->center[0],-this->center[1],-this->center[2]);
 
  this->widgetorientation->GetOrientationWXYZ(temp);
  //t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
  // Translation
  translate[0] = (1/this->WidgetScale)*this->pos[0];
  translate[1] = (1/this->WidgetScale)*this->pos[1];
  translate[2] = (1/this->WidgetScale)*this->pos[2];
  t->Translate(translate[0], translate[1], translate[2]);
  //t->RotateWXYZ(-temp[0],temp[1],temp[2],temp[3]);

  t->Translate(this->center[0],this->center[1],this->center[2]);
  this->widgetorientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
  this->orientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
}

void vtkAdvancedAngleWidget::GetTransformForActor(vtkTransform *t)
{
  //double translate[3];
  double temp[4];

  this->GetTransform(t);
  this->widgetorientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(-temp[0],temp[1],temp[2],temp[3]);
    /*
  // Translation
  translate[0] = this->pos[0];
  translate[1] = this->pos[1];
  translate[2] = this->pos[2];

  this->widgetorientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
  t->Translate(translate[0], translate[1], translate[2]);
  t->RotateWXYZ(-temp[0],temp[1],temp[2],temp[3]);

  t->Translate(this->center[0],this->center[1],this->center[2]);
  this->widgetorientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
  this->orientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(temp[0],temp[1],temp[2],temp[3]);
  this->widgetorientation->GetOrientationWXYZ(temp);
  t->RotateWXYZ(-temp[0],temp[1],temp[2],temp[3]);

  t->Translate(-this->center[0],-this->center[1],-this->center[2]);*/
}
