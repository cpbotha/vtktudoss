/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkOpenGLPolyDataMapper.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtktudExtendedOpenGLPolyDataMapper.h"

#include <math.h>
#include <stdio.h>

#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCommand.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkMatrix4x4.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkOpenGLRenderer.h"
#include "vtkPlane.h"
#include "vtkPlaneCollection.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkProperty.h"
#include "vtkTimerLog.h"
#include "vtkTriangle.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkOpenGLTexture.h"
#include "vtkImageData.h"

#include "vtkOpenGLExtensionManager.h"
#include "vtkgl.h"


# include "vtkOpenGL.h"

#include <math.h>

vtkCxxRevisionMacro(vtktudExtendedOpenGLPolyDataMapper, "$Revision: 1.00 $");
vtkStandardNewMacro(vtktudExtendedOpenGLPolyDataMapper);

// Construct empty object.
vtktudExtendedOpenGLPolyDataMapper::vtktudExtendedOpenGLPolyDataMapper()
{
}

// Destructor (don't call ReleaseGraphicsResources() since it is virtual
vtktudExtendedOpenGLPolyDataMapper::~vtktudExtendedOpenGLPolyDataMapper()
{
}

double distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
   return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}

void vtktudExtendedOpenGLPolyDataMapper::RenderPiece(vtkRenderer *ren, vtkActor *act)
{
  vtkPolyData *input= this->GetInput();
  vtkPlaneCollection *clipPlanes;
  vtkPlane *plane;
  int i, numClipPlanes;
  double planeEquation[4];

  glDepthMask(false);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE); //blended, soft, additive
  glLineWidth(2);

	/*double *pos;
	double origin[3];
	ren->GetRenderWindow()->Render();
	pos = ren->GetActiveCamera()->GetPosition();
	float d = sqrt(vtkMath::Distance2BetweenPoints(pos, origin));
	float fogColor[4]={0.3,0.3,0.3,0.0};
	glFogf( GL_FOG_START, d-30.0 );
	glFogf( GL_FOG_END, d+30.0 );
	glFogfv(GL_FOG_COLOR, fogColor);
	glFogf( GL_FOG_DENSITY, 0.5);
	glFogi( GL_FOG_MODE, GL_LINEAR );
	
	glEnable(GL_FOG );*/



  //glEnable(GL_FOG);
  //glFogi(GL_FOG_MODE, GL_EXP2);//  - rate of fade mode = GL_EXP, GL_EXP2 or GL_LINEAR 
  //GLfloat d = 0.5;
  //glFogf(GL_FOG_DENSITY,d);


		/*vtkOpenGLExtensionManager *extensions=vtkOpenGLExtensionManager::New();
		extensions->SetRenderWindow(ren->GetRenderWindow());
		int supports_GL_ARB_shader_objects=extensions->ExtensionSupported("GL_ARB_shader_objects");
		int supports_GL_ARB_vertex_shader=extensions->ExtensionSupported("GL_ARB_vertex_shader");
		int supports_GL_ARB_fragment_shader=extensions->ExtensionSupported("GL_ARB_fragment_shader");
		if(supports_GL_ARB_shader_objects && supports_GL_ARB_vertex_shader && supports_GL_ARB_fragment_shader)
		{
		extensions->LoadExtension("GL_ARB_vertex_shader");
		extensions->LoadExtension("GL_ARB_shader_objects");
		extensions->LoadExtension("GL_ARB_fragment_shader");
		}
		extensions->Delete();
		GLuint shader=vtk::CreateShaderObjectARB(vtkgl::FRAGMENT_SHADER_ARB);
*/

  //
  // make sure that we've been properly initialized
  //
  if (ren->GetRenderWindow()->CheckAbortStatus())
    {
    return;
    }
  
  if ( input == NULL ) 
    {
    vtkErrorMacro(<< "No input!");
    return;
    }
  else
    {
    this->InvokeEvent(vtkCommand::StartEvent,NULL);
    if (!this->Static)
      {
      input->Update();
      }
    this->InvokeEvent(vtkCommand::EndEvent,NULL);

    vtkIdType numPts = input->GetNumberOfPoints();
    if (numPts == 0)
      {
      vtkDebugMacro(<< "No points!");
      return;
      }
    } 

  if ( this->LookupTable == NULL )
    {
    this->CreateDefaultLookupTable();
    }

  // make sure our window is current
  ren->GetRenderWindow()->MakeCurrent();

  clipPlanes = this->ClippingPlanes;

  if (clipPlanes == NULL)
    {
    numClipPlanes = 0;
    }
  else
    {
    numClipPlanes = clipPlanes->GetNumberOfItems();
    if (numClipPlanes > 6)
      {
      vtkErrorMacro(<< "OpenGL guarantees at most 6 additional clipping planes");
      }
    }

  for (i = 0; i < numClipPlanes; i++)
    {
     glEnable(static_cast<GLenum>(GL_CLIP_PLANE0+i));
    }

  if ( clipPlanes )
    {
    vtkMatrix4x4 *actorMatrix = vtkMatrix4x4::New();
    act->GetMatrix( actorMatrix );
    actorMatrix->Invert();
    
    double origin[4], normal[3], point[4];
    
    for (i = 0; i < numClipPlanes; i++)
      {    
      plane = static_cast<vtkPlane *>(clipPlanes->GetItemAsObject(i));
      
      plane->GetOrigin(origin);
      plane->GetNormal(normal);
      
      point[0] = origin[0] + normal[0];
      point[1] = origin[1] + normal[1];
      point[2] = origin[2] + normal[2];
      
      origin[3] = point[3] = 1.0;
      
      actorMatrix->MultiplyPoint( origin, origin );
      actorMatrix->MultiplyPoint( point, point );
      
      if ( origin[3] != 1.0 )
        {
        origin[0] /= origin[3];
        origin[1] /= origin[3];
        origin[2] /= origin[3];
        }
      
      if ( point[3] != 1.0 )
        {
        point[0] /= point[3];
        point[1] /= point[3];
        point[2] /= point[3];
        }
      
      normal[0] = point[0] - origin[0];
      normal[1] = point[1] - origin[1];
      normal[2] = point[2] - origin[2];
      
      planeEquation[0] = normal[0];
      planeEquation[1] = normal[1];
      planeEquation[2] = normal[2];
      planeEquation[3] = -(planeEquation[0]*origin[0]+
                           planeEquation[1]*origin[1]+
                           planeEquation[2]*origin[2]);
      glClipPlane(static_cast<GLenum>(GL_CLIP_PLANE0+i),planeEquation);
      }
    
    actorMatrix->Delete();  
    }

  // For vertex coloring, this sets this->Colors as side effect.
  // For texture map coloring, this sets ColorCoordinates
  // and ColorTextureMap as a side effect.
  // I moved this out of the conditional because it is fast.
  // Color arrays are cached. If nothing has changed, 
  // then the scalars do not have to be regenerted. 
  this->MapScalars(act->GetProperty()->GetOpacity());
  // If we are coloring by texture, then load the texture map.
  if (this->ColorTextureMap)
    {
    if (this->InternalColorTexture == 0)
      {
      this->InternalColorTexture = vtkOpenGLTexture::New();
      this->InternalColorTexture->RepeatOff();
      }
    this->InternalColorTexture->SetInput(this->ColorTextureMap);
    // Keep color from interacting with texture.
    float info[4];
    info[0] = info[1] = info[2] = info[3] = 1.0;
    glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, info );
    }

  //
  // if something has changed regenerate colors and display lists
  // if required
  //
  int noAbort=1;
  if ( this->GetMTime() > this->BuildTime || 
       input->GetMTime() > this->BuildTime ||
       act->GetProperty()->GetMTime() > this->BuildTime ||
       ren->GetRenderWindow() != this->LastWindow)
    {
    if (!this->ImmediateModeRendering && 
        !this->GetGlobalImmediateModeRendering())
      {
      this->ReleaseGraphicsResources(ren->GetRenderWindow());
      this->LastWindow = ren->GetRenderWindow();
      
      // If we are coloring by texture, then load the texture map.
      // Use Map as indicator, because texture hangs around.
      if (this->ColorTextureMap)
        {
        this->InternalColorTexture->Load(ren);
        }
      
      // get a unique display list id
      this->ListId = glGenLists(1);
      glNewList(this->ListId,GL_COMPILE);

      noAbort = this->Draw(ren,act);
      glEndList();

      // Time the actual drawing
      this->Timer->StartTimer();
      glCallList(this->ListId);
      this->Timer->StopTimer();      
      }
    else
      {
      this->ReleaseGraphicsResources(ren->GetRenderWindow());
      this->LastWindow = ren->GetRenderWindow();
      }
    if (noAbort)
      {
      this->BuildTime.Modified();
      }
    }
  // if nothing changed but we are using display lists, draw it
  else
    {
    if (!this->ImmediateModeRendering && 
        !this->GetGlobalImmediateModeRendering())
      {
      // If we are coloring by texture, then load the texture map.
      // Use Map as indicator, because texture hangs around.
      if (this->ColorTextureMap)
        {
        this->InternalColorTexture->Load(ren);
        }

      // Time the actual drawing
      this->Timer->StartTimer();
      glCallList(this->ListId);
      this->Timer->StopTimer();      
      }
    }
   
  // if we are in immediate mode rendering we always
  // want to draw the primitives here
  if (this->ImmediateModeRendering ||
      this->GetGlobalImmediateModeRendering())
    {
    // If we are coloring by texture, then load the texture map.
    // Use Map as indicator, because texture hangs around.
    if (this->ColorTextureMap)
      {
      this->InternalColorTexture->Load(ren);
      }
    // Time the actual drawing
    this->Timer->StartTimer();
    this->Draw(ren,act);
    this->Timer->StopTimer();      
    }

  this->TimeToDraw = this->Timer->GetElapsedTime();

  // If the timer is not accurate enough, set it to a small
  // time so that it is not zero
  if ( this->TimeToDraw == 0.0 )
    {
    this->TimeToDraw = 0.0001;
    }

  for (i = 0; i < numClipPlanes; i++)
    {
    glDisable(static_cast<GLenum>(GL_CLIP_PLANE0+i));
    }

  glDepthMask(true);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_BLEND);
  glDisable(GL_FOG);
  //glBlendFunc(GL_ONE, GL_ONE);//GL_SRC_ALPHA, GL_DST_COLOR);
  glLineWidth(1);
}

