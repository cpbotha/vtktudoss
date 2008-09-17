#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/packages)

SET (STLib_SRCS
   vtkCPTDistanceField.cxx
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (vtktudossSTLib ${STLib_SRCS})
TARGET_LINK_LIBRARIES(vtktudossSTLib vtkFiltering vtkGraphics vtkCommon)

# Create the vtktudossSTLib Tcl library
IF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudossSTLibTCL STLibTCL_SRCS
  "${STLib_SRCS}" "")
  ADD_LIBRARY (vtktudossSTLibTCL ${STLibTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossSTLibTCL vtktudossSTLib)
ENDIF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)

# Create the vtktudossSTLib Python library.
IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  VTK_WRAP_PYTHON3 (vtktudossSTLibPython STLibPython_SRCS "${STLib_SRCS}")
  ADD_LIBRARY (vtktudossSTLibPythonD ${STLibPython_SRCS})
  ADD_LIBRARY (vtktudossSTLibPython MODULE
  vtktudossSTLibPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudossSTLibPython vtktudossSTLibPythonD)
  TARGET_LINK_LIBRARIES (vtktudossSTLibPythonD
                         vtkFilteringPythonD
                         vtkGraphicsPythonD
                         vtkCommonPythonD
                         vtktudossSTLib)
  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(vtktudossSTLibPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)

# Create the vtktudossSTLib Java library.
IF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudossSTLibJava ImagingJava_SRCS "${STLib_SRCS}")
  ADD_LIBRARY (vtktudossSTLibJava SHARED ${STLibava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossSTLibJava vtktudossSTLib)
ENDIF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
