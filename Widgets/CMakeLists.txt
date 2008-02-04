#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Widgets_SRCS
   vtkEllipseRepresentation.cxx
   vtkEllipseWidget.cxx
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (vtktudossWidgets ${Widgets_SRCS})
TARGET_LINK_LIBRARIES(vtktudossWidgets vtktudossGraphics vtkGraphics vtkRendering vtkFiltering vtkHybrid vtkWidgets vtkCommon)

# Create the vtktudossWidgets Tcl library
IF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudossWidgetsTCL WidgetsTCL_SRCS
  "${Widgets_SRCS}" "")
  ADD_LIBRARY (vtktudossWidgetsTCL ${WidgetsTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossWidgetsTCL vtkFilteringTCL vtkWidgetsTCL vtktudossWidgets)
ENDIF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)

# Create the vtktudossWidgets Python library.
IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  VTK_WRAP_PYTHON3 (vtktudossWidgetsPython WidgetsPython_SRCS "${Widgets_SRCS}")
  ADD_LIBRARY (vtktudossWidgetsPythonD ${WidgetsPython_SRCS})
  ADD_LIBRARY (vtktudossWidgetsPython MODULE
  vtktudossWidgetsPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudossWidgetsPython vtktudossWidgetsPythonD)
  TARGET_LINK_LIBRARIES (vtktudossWidgetsPythonD
			             vtkGraphicsPythonD
			 vtkRenderingPythonD
			 vtkFilteringPythonD
                         vtkHybridPythonD
                         vtkWidgetsPythonD 
                         vtktudossWidgets)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)

# Create the vtktudossWidgets Java library.
IF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudossWidgetsJava ImagingJava_SRCS "${Widgets_SRCS}")
  ADD_LIBRARY (vtktudossWidgetsJava SHARED ${Widgetsava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossWidgetsJava vtktudossWidgets)
ENDIF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
