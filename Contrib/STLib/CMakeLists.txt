#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/packages)

SET (STLib_SRCS
   vtkCPTDistanceField.cxx)

# PROJECT_NAME = VTKTUDOSS
SET (KIT_NAME vtktudossSTLib)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (${KIT_NAME} ${STLib_SRCS})
TARGET_LINK_LIBRARIES(${KIT_NAME} ${VTK_LIBRARIES})

IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  VTK_WRAP_PYTHON3 (${KIT_NAME}Python STLibPython_SRCS "${STLib_SRCS}")
  ADD_LIBRARY (${KIT_NAME}PythonD ${STLibPython_SRCS})
  TARGET_LINK_LIBRARIES(${KIT_NAME}PythonD ${KIT_NAME} vtkImagingCorePythonD)

  # I want to use this new python_add_module() macro instead, but could
  # not yet get it to build without errors.
  # python_add_module(${KIT_NAME}Python ${KIT_NAME}PythonInit.cxx)
  add_library(${KIT_NAME}Python MODULE ${KIT_NAME}PythonInit.cxx)
  target_link_libraries(${KIT_NAME}Python ${KIT_NAME}PythonD)

  # it should be ${KIT_NAME}Python and NOT libvtk...
  # else it will try to find PyInit_lib${KIT_NAME}Python when we import
  SET_TARGET_PROPERTIES(${KIT_NAME}Python PROPERTIES PREFIX "")

  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(${KIT_NAME}Python PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
