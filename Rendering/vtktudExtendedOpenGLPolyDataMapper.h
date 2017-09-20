#ifndef __vtktudExtendedOpenGLPolyDataMapper_h
#define __vtktudExtendedOpenGLPolyDataMapper_h

#include "vtkOpenGLPolyDataMapper.h"
#include "vtkOpenGL.h"

class vtkCellArray;
class vtkPoints;
class vtkProperty;
class vtkRenderWindow;
class vtkOpenGLRenderer;
class vtkOpenGLTexture;


class VTK_EXPORT vtktudExtendedOpenGLPolyDataMapper : public vtkOpenGLPolyDataMapper
{
public:
  static vtktudExtendedOpenGLPolyDataMapper *New();
  vtkTypeMacro(vtktudExtendedOpenGLPolyDataMapper,vtkOpenGLPolyDataMapper);


  virtual void RenderPiece(vtkRenderer *ren, vtkActor *a);

protected:
  vtktudExtendedOpenGLPolyDataMapper();
  ~vtktudExtendedOpenGLPolyDataMapper();

private:
	double distance(double x1, double y1, double z1, double x2, double y2, double z2);
};

#endif
