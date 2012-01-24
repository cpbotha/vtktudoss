#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Graphics_SRCS
   vtkEllipseSource.cxx
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (vtktudossGraphics ${Graphics_SRCS})
TARGET_LINK_LIBRARIES(vtktudossGraphics vtkRendering vtkGraphics)

# Create the vtktudossGraphics Tcl library
IF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudossGraphicsTCL GraphicsTCL_SRCS
  "${Graphics_SRCS}" "")
  ADD_LIBRARY (vtktudossGraphicsTCL ${GraphicsTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossGraphicsTCL vtkFilteringTCL vtkGraphicsTCL vtktudossGraphics)
ENDIF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)

# Create the vtktudossGraphics Python library.
IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  VTK_WRAP_PYTHON3 (vtktudossGraphicsPython GraphicsPython_SRCS "${Graphics_SRCS}")
  ADD_LIBRARY (vtktudossGraphicsPythonD ${GraphicsPython_SRCS})
  ADD_LIBRARY (vtktudossGraphicsPython MODULE
  vtktudossGraphicsPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudossGraphicsPython vtktudossGraphicsPythonD)
  TARGET_LINK_LIBRARIES (vtktudossGraphicsPythonD
			             vtkRenderingPythonD
                         vtkGraphicsPythonD 
                         vtktudossGraphics)
  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(vtktudossGraphicsPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)

# Create the vtktudossGraphics Java library.
IF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudossGraphicsJava ImagingJava_SRCS "${Graphics_SRCS}")
  ADD_LIBRARY (vtktudossGraphicsJava SHARED ${Graphicsava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossGraphicsJava vtktudossGraphics)
ENDIF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
