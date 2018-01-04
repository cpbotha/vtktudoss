#include "vtkCustomWidget.h"

#include "vtkSphereSource.h"
#include "vtkTransform.h"

// we're using snprintf to convert from integer to string representation
#ifdef _MSC_VER  //all MS compilers define this (version)
     #define snprintf _snprintf
#endif


vtkStandardNewMacro(vtkCustomWidget);

// Constructor
vtkCustomWidget::vtkCustomWidget()
{
}

// Destructor
vtkCustomWidget::~vtkCustomWidget()
{
}

void vtkCustomWidget::PlaceWidget(double bds[6])
{
}

void vtkCustomWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
