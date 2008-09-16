# CMakeLists.txt file for VTKTUDOSS
# copyright TU Delft Data Visualization Group
# Berkeley Licensing

#
# Set the project name.
#

PROJECT (VTKTUDOSS)

#
# Assume everything is set up correctly for build.
#
# If some requirements are not met, set it to 0.
# In that case, the build won't be done.
#

SET (VTKTUDOSS_CAN_BUILD 1)

#
# Load CMake commands that you probably should not modify.
#

INCLUDE (${VTKTUDOSS_SOURCE_DIR}/CMakeOptions.cmake)

#
# If everything is OK, then go into the sub directories and build.
#

IF (VTKTUDOSS_CAN_BUILD)

  #
  # This is the list of kits available in VTKTUDOSS. Each kit corresponds to 
  # a sub-directory with the same name holding your local classes. Sorting 
  # classes by kit like VTK does (Common, Rendering, Filtering, Imaging, 
  # IO, etc.) is a good thing and prevents numerous dependencies problems.
  #

  # Core (non-optional kits)
  SET(VTKTUDOSS_KITS_CORE
    Graphics
    Widgets
  )

  # Optional kits
  # Syntax: Dirname "Description"
  SET(VTKTUDOSS_KITS_OPTIONAL
  )

  # Contrib (3rd party) kits
  # Syntax: Dirname "Description"
  SET(VTKTUDOSS_KITS_CONTRIB
      STLib "Mauch's STLib (CPT)"
  )

  #
  # --- You shouldn't need to modify anything below this line ---
  #

  # Keep lists of the kits built / not built in the current configuration
  SET(VTKTUDOSS_KITS_BUILT ${VTKTUDOSS_KITS_CORE})
  SET(VTKTUDOSS_KITS_NOTBUILT)

  # Process Core kits
  FOREACH(kit ${VTKTUDOSS_KITS_CORE})
    SUBDIRS(${kit})
    # Add it to the include dirs, otherwise building the wrappers may fail
    INCLUDE_DIRECTORIES(${VTKTUDOSS_SOURCE_DIR}/${kit})
  ENDFOREACH(kit)

  # Process optional kits
  LIST(LENGTH VTKTUDOSS_KITS_OPTIONAL nkits)
  IF(LENGTH) # we can only do this if there is at least 1 optional kit
  MATH(EXPR nkits "(${nkits} / 2) - 1" )
  FOREACH(i RANGE ${nkits})
    MATH(EXPR j "${i} * 2")
    MATH(EXPR k "${j} + 1")
    LIST(GET VTKTUDOSS_KITS_OPTIONAL ${j} kit)
    LIST(GET VTKTUDOSS_KITS_OPTIONAL ${k} kit_desc)
    STRING(TOUPPER ${kit} kit_uc)
    
    # Create an option for this kit and build it if selected
    OPTION(BUILD_${kit_uc} "${kit_desc}" OFF)
    IF (BUILD_${kit_uc})
      SUBDIRS(${kit})
      # Add it to the include dirs, otherwise building the wrappers may fail
      INCLUDE_DIRECTORIES(${VTKTUDOSS_SOURCE_DIR}/${kit})
      # Add to the list of kits built
      SET(VTKTUDOSS_KITS_BUILT ${VTKTUDOSS_KITS_BUILT} ${kit})
    ELSE (BUILD_${kit_uc})
      SET(VTKTUDOSS_KITS_NOTBUILT ${VTKTUDOSS_KITS_NOTBUILT} ${kit})
    ENDIF (BUILD_${kit_uc})
  ENDFOREACH(i)
  ENDIF(LENGTH)

  # Process Contrib kits
  OPTION(BUILD_CONTRIB "Build contributed code." OFF)
  IF(BUILD_CONTRIB)
    LIST(LENGTH VTKTUDOSS_KITS_CONTRIB nkits)
    MATH(EXPR nkits "(${nkits} / 2) - 1" )
    FOREACH(i RANGE ${nkits})
      MATH(EXPR j "${i} * 2")
      MATH(EXPR k "${j} + 1")
      LIST(GET VTKTUDOSS_KITS_CONTRIB ${j} kit)
      LIST(GET VTKTUDOSS_KITS_CONTRIB ${k} kit_desc)
      STRING(TOUPPER ${kit} kit_uc)
      
      # Create an option for this kit and build it if selected
      OPTION(BUILD_CONTRIB_${kit_uc} "Contrib: build ${kit_desc}." OFF)
      IF(BUILD_CONTRIB_${kit_uc})
        SUBDIRS(Contrib/${kit})
        # Add it to the include dirs, required for building wrappers
        INCLUDE_DIRECTORIES(${VTKTUDOSS_SOURCE_DIR}/Contrib/${kit})
        # Add to the list of kits built
        SET(VTKTUDOSS_KITS_BUILT ${VTKTUDOSS_KITS_BUILT} "${kit}")
      ELSE(BUILD_CONTRIB_${kit_uc})
        SET(VTKTUDOSS_KITS_NOTBUILT ${VTKTUDOSS_KITS_NOTBUILT} "${kit}")
      ENDIF(BUILD_CONTRIB_${kit_uc})
    ENDFOREACH(i)
  ENDIF(BUILD_CONTRIB)

  #
  # Build examples too ?
  #

  # uncomment when you start adding Examples!
  # OPTION(BUILD_EXAMPLES "Build examples." ON)
  # IF (BUILD_EXAMPLES)
  #   SUBDIRS(Examples)
  # ENDIF (BUILD_EXAMPLES)


  #
  # Build Utilities (Doxygen)
  #

  #SUBDIRS(Utilities)

  INCLUDE_DIRECTORIES(${VTKTUDOSS_BINARY_DIR})

  #
  # Generate files for using VTKTUDOSS from another cmake project
  #

  SET(VTKTUDOSS_INCLUDE_DIRS_CONFIG "${VTKTUDOSS_SOURCE_DIR}")
  SET(VTKTUDOSS_LIBRARY_DIRS_CONFIG "${LIBRARY_OUTPUT_PATH}")
  SET(VTKTUDOSS_USE_FILE "${VTKTUDOSS_BINARY_DIR}/UseVTKTUDOSS.cmake")
  SET(VTKTUDOSS_KITS_FORCONFIG)
  FOREACH(kit ${VTKTUDOSS_KITS_BUILT})
    STRING(TOUPPER ${kit} kit_uc)
    # Prettier formatting for Contrib kits
    IF(${kit} MATCHES "^Contrib")
      STRING(REGEX REPLACE "^CONTRIB(.*)" "CONTRIB_\\1" kit_uc ${kit_uc})
    ENDIF(${kit} MATCHES "^Contrib")
    SET(VTKTUDOSS_KITS_FORCONFIG ${VTKTUDOSS_KITS_FORCONFIG} ${kit_uc})
  ENDFOREACH(kit)
  CONFIGURE_FILE( 
    ${VTKTUDOSS_SOURCE_DIR}/UseVTKTUDOSS.cmake.in
    ${VTKTUDOSS_BINARY_DIR}/UseVTKTUDOSS.cmake
    COPYONLY IMMEDIATE
  )
  CONFIGURE_FILE(
    ${VTKTUDOSS_SOURCE_DIR}/VTKTUDOSSConfig.cmake.in
    ${VTKTUDOSS_BINARY_DIR}/VTKTUDOSSConfig.cmake
    @ONLY IMMEDIATE
  )

  #
  # Generate python import module
  #

  IF(VTKTUDOSS_WRAP_PYTHON)
    SET(VTKTUDOSS_KITS_FORPYTHON ${VTKTUDOSS_KITS_BUILT})
    CONFIGURE_FILE(
      "${VTKTUDOSS_SOURCE_DIR}/Wrapping/Python/vtktudoss.py.in" 
      "${VTKTUDOSS_BINARY_DIR}/Wrapping/Python/vtktudoss.py"
      @ONLY
    )
  ENDIF(VTKTUDOSS_WRAP_PYTHON)

  #
  # Disable deprecation warnings for standard C and STL functions in VS2005
  # and later
  #

  IF(CMAKE_COMPILER_2005)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_DEPRECATE)
  ENDIF(CMAKE_COMPILER_2005) 

ENDIF (VTKTUDOSS_CAN_BUILD)
