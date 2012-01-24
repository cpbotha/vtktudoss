#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (CompVis_SRCS
   vtkCVFixedPointVolumeRayCastMapper
   vtkFixedPointVolumeRayCastCVHelper
   vtkCVImageMapToColors
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (vtktudossCompVis ${CompVis_SRCS})
TARGET_LINK_LIBRARIES(vtktudossCompVis vtkImaging vtkVolumeRendering)

# Create the vtktudossCompVis Tcl library
IF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudossCompVisTCL CompVisTCL_SRCS
  "${CompVis_SRCS}" "")
  ADD_LIBRARY (vtktudossCompVisTCL ${CompVisTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossCompVisTCL vtkFilteringTCL vtkCompVisTCL vtktudossCompVis)
ENDIF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)

# Create the vtktudossCompVis Python library.
IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  VTK_WRAP_PYTHON3 (vtktudossCompVisPython CompVisPython_SRCS "${CompVis_SRCS}")
  ADD_LIBRARY (vtktudossCompVisPythonD ${CompVisPython_SRCS})
  ADD_LIBRARY (vtktudossCompVisPython MODULE
  vtktudossCompVisPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudossCompVisPython vtktudossCompVisPythonD)
  TARGET_LINK_LIBRARIES (vtktudossCompVisPythonD
			             vtkRenderingPythonD
                         vtkGraphicsPythonD 
                         vtkVolumeRenderingPythonD
                         vtktudossCompVis)
  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(vtktudossCompVisPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)

# Create the vtktudossCompVis Java library.
IF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudossCompVisJava ImagingJava_SRCS "${CompVis_SRCS}")
  ADD_LIBRARY (vtktudossCompVisJava SHARED ${CompVisava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossCompVisJava vtktudossCompVis)
ENDIF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)

