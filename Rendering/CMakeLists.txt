#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Rendering_SRCS
   vtktudExtendedOpenGLPolyDataMapper.cxx
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (vtktudossRendering ${Rendering_SRCS})
TARGET_LINK_LIBRARIES(vtktudossRendering vtkRendering vtkRendering)

# cpbotha: I have absolutely NO idea why this has become necessary all of
# a sudden with VTK 5.8 and MSVC 2008 
# if I don't have this, I get unsatisfied link errors on bunches of 
# __imp_glEnable and __imp_gl* calls. I also tried ${OPENGL_gl_LIBRARY}
# and ${OPENGL_LIBRARIES}, but to no avail.
IF(MSVC)
    TARGET_LINK_LIBRARIES(vtktudossRendering "opengl32.lib")
ENDIF(MSVC)

# Create the vtktudossRendering Tcl library
IF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudossRenderingTCL RenderingTCL_SRCS
  "${Rendering_SRCS}" "")
  ADD_LIBRARY (vtktudossRenderingTCL ${RenderingTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossRenderingTCL vtkFilteringTCL vtkRenderingTCL vtktudossRendering)
ENDIF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)

# Create the vtktudossRendering Python library.
IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  VTK_WRAP_PYTHON3 (vtktudossRenderingPython RenderingPython_SRCS "${Rendering_SRCS}")
  ADD_LIBRARY (vtktudossRenderingPythonD ${RenderingPython_SRCS})
  ADD_LIBRARY (vtktudossRenderingPython MODULE
  vtktudossRenderingPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudossRenderingPython vtktudossRenderingPythonD)
  TARGET_LINK_LIBRARIES (vtktudossRenderingPythonD
			 vtkRenderingPythonD
                         vtktudossRendering)
  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(vtktudossRenderingPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)

# Create the vtktudossRendering Java library.
IF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudossRenderingJava ImagingJava_SRCS "${Rendering_SRCS}")
  ADD_LIBRARY (vtktudossRenderingJava SHARED ${Renderingava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossRenderingJava vtktudossRendering)
ENDIF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
