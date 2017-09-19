#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Filtering_SRCS
   vtkCKdTree.cxx
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossFiltering C++ library
ADD_LIBRARY (vtktudossFiltering ${Filtering_SRCS})
TARGET_LINK_LIBRARIES(vtktudossFiltering vtkFiltering) 

# Create the vtktudossFiltering Tcl library
IF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)
  VTK_WRAP_TCL3 (vtktudossFilteringTCL FilteringTCL_SRCS
  "${Filtering_SRCS}" "")
  ADD_LIBRARY (vtktudossFilteringTCL ${FilteringTCL_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossFilteringTCL vtkFilteringTCL vtktudossFiltering)
ENDIF (VTK_WRAP_TCL AND VTKTUDOSS_WRAP_TCL)

# Create the vtktudossFiltering Python library.
IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  VTK_WRAP_PYTHON3 (vtktudossFilteringPython FilteringPython_SRCS "${Filtering_SRCS}")
  ADD_LIBRARY (vtktudossFilteringPythonD ${FilteringPython_SRCS})
  ADD_LIBRARY (vtktudossFilteringPython MODULE
  vtktudossFilteringPythonInit.cxx)
  TARGET_LINK_LIBRARIES (vtktudossFilteringPython vtktudossFilteringPythonD)
  TARGET_LINK_LIBRARIES (vtktudossFilteringPythonD
			             vtkFilteringPythonD
                         vtktudossFiltering)
  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(vtktudossFilteringPython PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)

# Create the vtktudossFiltering Java library.
IF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
  VTK_WRAP_JAVA3 (vtktudossFilteringJava ImagingJava_SRCS "${Filtering_SRCS}")
  ADD_LIBRARY (vtktudossFilteringJava SHARED ${Filteringava_SRCS})
  TARGET_LINK_LIBRARIES (vtktudossFilteringJava vtktudossFiltering)
ENDIF (VTK_WRAP_JAVA AND VTKTUDOSS_WRAP_JAVA)
