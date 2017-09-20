#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (KIT_SRCS
   vtkEllipseSource.cxx
   )

SET (KIT Graphics)
SET (BASENAME ${PROJECT_NAME}${KIT})

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossGraphics C++ library
ADD_LIBRARY (${BASENAME} ${KIT_SRCS})
TARGET_LINK_LIBRARIES(${BASENAME} ${VTK_LIBRARIES})

# install the shared library
INSTALL(TARGETS ${BASENAME}
  RUNTIME DESTINATION ${BIN_DESTINATION}
  ARCHIVE DESTINATION ${LIB_DESTINATION}
  LIBRARY DESTINATION ${LIB_DESTINATION})

IF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  VTK_WRAP_PYTHON3 (${BASENAME}Python KITPython_SRCS "${KIT_SRCS}")
  ADD_LIBRARY (${BASENAME}PythonD ${KITPython_SRCS})
  # We need _PyvtkPolyDataAlgorithm_ClassNew from vtkCommonExecutionModelPythonD
  TARGET_LINK_LIBRARIES(${BASENAME}PythonD ${BASENAME} vtkCommonExecutionModelPythonD)

  # I want to use this new python_add_module() macro instead, but could
  # not yet get it to build without errors.
  # python_add_module(${BASENAME}Python ${BASENAME}PythonInit.cxx)
  add_library(${BASENAME}Python MODULE ${BASENAME}PythonInit.cxx)
  target_link_libraries(${BASENAME}Python ${BASENAME}PythonD)

  # it should be ${BASENAME}Python and NOT libvtk...
  # else it will try to find PyInit_lib${BASENAME}Python when we import
  SET_TARGET_PROPERTIES(${BASENAME}Python PROPERTIES PREFIX "")

  # Python extension modules on Windows must have the extension ".pyd"
  # instead of ".dll" as of Python 2.5.  Older python versions do support
  # this suffix.
  IF(WIN32 AND NOT CYGWIN)
    SET_TARGET_PROPERTIES(${BASENAME}Python PROPERTIES SUFFIX ".pyd")
  ENDIF(WIN32 AND NOT CYGWIN)

  # the python backend shared library ()
  INSTALL(TARGETS "${BASENAME}PythonD"
    RUNTIME DESTINATION ${BIN_DESTINATION}
    ARCHIVE DESTINATION ${LIB_DESTINATION}
    LIBRARY DESTINATION ${LIB_DESTINATION})

  # the frontend (importable) python lib
  INSTALL(TARGETS "${BASENAME}Python"
    LIBRARY DESTINATION "${INSTALL_PYTHON_MODULE_DIR}/${PROJECT_NAME}")
  
ENDIF (VTK_WRAP_PYTHON AND VTKTUDOSS_WRAP_PYTHON)
