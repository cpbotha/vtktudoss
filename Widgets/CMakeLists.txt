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

# Create the vtktudWidgets C++ library
ADD_LIBRARY (vtktudWidgets ${Widgets_SRCS})
TARGET_LINK_LIBRARIES(vtktudWidgets vtktudGraphics vtktudCommon vtkGraphics vtkRendering vtkFiltering vtkHybrid vtkWidgets vtkCommon)

# Create the vtktudWidgets Tcl library
IF (VTK_WRAP_TCL AND VTKTUD_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudWidgetsTCL WidgetsTCL_SRCS
  "${Widgets_SRCS}" "")
  ADD_LIBRARY (vtktudWidgetsTCL ${WidgetsTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudWidgetsTCL vtktudCommonTCL vtkFilteringTCL vtkWidgetsTCL vtktudWidgets)
ENDIF (VTK_WRAP_TCL AND VTKTUD_WRAP_TCL)

# Create the vtktudWidgets Python library.
IF (VTK_WRAP_PYTHON AND VTKTUD_WRAP_PYTHON)
  VTK_WRAP_PYTHON3 (vtktudWidgetsPython WidgetsPython_SRCS "${Widgets_SRCS}")
  ADD_LIBRARY (vtktudWidgetsPythonD ${WidgetsPython_SRCS})
  ADD_LIBRARY (vtktudWidgetsPython MODULE
  vtktudWidgetsPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudWidgetsPython vtktudWidgetsPythonD)
  TARGET_LINK_LIBRARIES (vtktudWidgetsPythonD
                         vtktudCommonPythonD
			 vtkGraphicsPythonD
			 vtkRenderingPythonD
			 vtkFilteringPythonD
                         vtkHybridPythonD
                         vtkWidgetsPythonD 
                         vtktudWidgets)
ENDIF (VTK_WRAP_PYTHON AND VTKTUD_WRAP_PYTHON)

# Create the vtktudWidgets Java library.
IF (VTK_WRAP_JAVA AND VTKTUD_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudWidgetsJava ImagingJava_SRCS "${Widgets_SRCS}")
  ADD_LIBRARY (vtktudWidgetsJava SHARED ${Widgetsava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudWidgetsJava vtktudCommonJava vtktudWidgets)
ENDIF (VTK_WRAP_JAVA AND VTKTUD_WRAP_JAVA)
