#
# Set the project name.
#

PROJECT (VTKTUD)

#
# Assume everything is set up correctly for build.
#
# If some requirements are not met, set it to 0.
# In that case, the build won't be done.
#

SET (VTKTUD_CAN_BUILD 1)

#
# Load CMake commands that you probably should not modify.
#

INCLUDE (${VTKTUD_SOURCE_DIR}/CMakeOptions.cmake)

#
# If everything is OK, then go into the sub directories and build.
#

IF (VTKTUD_CAN_BUILD)

  #
  # This is the list of kits available in VTKTUD. Each kit corresponds to 
  # a sub-directory with the same name holding your local classes. Sorting 
  # classes by kit like VTK does (Common, Rendering, Filtering, Imaging, 
  # IO, etc.) is a good thing and prevents numerous dependencies problems.
  #

  # Core (non-optional kits)
  SET(VTKTUD_KITS_CORE
    Common
    Imaging
    Graphics
  )

  # Optional kits
  # Syntax: Dirname "Description"
  SET(VTKTUD_KITS_OPTIONAL
    Collision "Opcode-based Collision"
    OpticFlow "Optic Flow"
    Comparative "Comparative Visualization"
    Widgets "Extra new-style (post 5.0) widgets"
  )

  # Contrib (3rd party) kits
  # Syntax: Dirname "Description"
  SET(VTKTUD_KITS_CONTRIB
    Curvature "curvature classes"
    SPCView "Weigle's SPCview curvature texturizer"
  )

  #
  # --- You shouldn't need to modify anything below this line ---
  #

  # Keep lists of the kits built / not built in the current configuration
  SET(VTKTUD_KITS_BUILT ${VTKTUD_KITS_CORE})
  SET(VTKTUD_KITS_NOTBUILT)

  # Process Core kits
  FOREACH(kit ${VTKTUD_KITS_CORE})
    SUBDIRS(${kit})
    # Add it to the include dirs, otherwise building the wrappers may fail
    INCLUDE_DIRECTORIES(${VTKTUD_SOURCE_DIR}/${kit})
  ENDFOREACH(kit)

  # Process optional kits
  LIST(LENGTH VTKTUD_KITS_OPTIONAL nkits)
  MATH(EXPR nkits "(${nkits} / 2) - 1" )
  FOREACH(i RANGE ${nkits})
    MATH(EXPR j "${i} * 2")
    MATH(EXPR k "${j} + 1")
    LIST(GET VTKTUD_KITS_OPTIONAL ${j} kit)
    LIST(GET VTKTUD_KITS_OPTIONAL ${k} kit_desc)
    STRING(TOUPPER ${kit} kit_uc)
    
    # Create an option for this kit and build it if selected
    OPTION(BUILD_${kit_uc} "${kit_desc}" OFF)
    IF (BUILD_${kit_uc})
      SUBDIRS(${kit})
      # Add it to the include dirs, otherwise building the wrappers may fail
      INCLUDE_DIRECTORIES(${VTKTUD_SOURCE_DIR}/${kit})
      # Add to the list of kits built
      SET(VTKTUD_KITS_BUILT ${VTKTUD_KITS_BUILT} ${kit})
    ELSE (BUILD_${kit_uc})
      SET(VTKTUD_KITS_NOTBUILT ${VTKTUD_KITS_NOTBUILT} ${kit})
    ENDIF (BUILD_${kit_uc})
  ENDFOREACH(i)

  # Process Contrib kits
  OPTION(BUILD_CONTRIB "Build contributed code." OFF)
  IF(BUILD_CONTRIB)
    LIST(LENGTH VTKTUD_KITS_CONTRIB nkits)
    MATH(EXPR nkits "(${nkits} / 2) - 1" )
    FOREACH(i RANGE ${nkits})
      MATH(EXPR j "${i} * 2")
      MATH(EXPR k "${j} + 1")
      LIST(GET VTKTUD_KITS_CONTRIB ${j} kit)
      LIST(GET VTKTUD_KITS_CONTRIB ${k} kit_desc)
      STRING(TOUPPER ${kit} kit_uc)
      
      # Create an option for this kit and build it if selected
      OPTION(BUILD_CONTRIB_${kit_uc} "Contrib: build ${kit_desc}." OFF)
      IF(BUILD_CONTRIB_${kit_uc})
        SUBDIRS(Contrib/${kit})
        # Add it to the include dirs, required for building wrappers
        INCLUDE_DIRECTORIES(${VTKTUD_SOURCE_DIR}/Contrib/${kit})
        # Add to the list of kits built
        SET(VTKTUD_KITS_BUILT ${VTKTUD_KITS_BUILT} "Contrib${kit}")
      ELSE(BUILD_CONTRIB_${kit_uc})
        SET(VTKTUD_KITS_NOTBUILT ${VTKTUD_KITS_NOTBUILT} "Contrib${kit}")
      ENDIF(BUILD_CONTRIB_${kit_uc})
    ENDFOREACH(i)
  ENDIF(BUILD_CONTRIB)

  #
  # Build examples too ?
  #

  OPTION(BUILD_EXAMPLES "Build examples." ON)
  IF (BUILD_EXAMPLES)
    SUBDIRS(Examples)
  ENDIF (BUILD_EXAMPLES)


  #
  # Build Utilities (Doxygen)
  #

  SUBDIRS(Utilities)

  INCLUDE_DIRECTORIES(${VTKTUD_BINARY_DIR})

  #
  # Generate files for using VTKTUD from another cmake project
  #

  SET(VTKTUD_INCLUDE_DIRS_CONFIG "${VTKTUD_SOURCE_DIR}")
  SET(VTKTUD_LIBRARY_DIRS_CONFIG "${LIBRARY_OUTPUT_PATH}")
  SET(VTKTUD_USE_FILE "${VTKTUD_BINARY_DIR}/UseVTKTUD.cmake")
  SET(VTKTUD_KITS_FORCONFIG)
  FOREACH(kit ${VTKTUD_KITS_BUILT})
    STRING(TOUPPER ${kit} kit_uc)
    # Prettier formatting for Contrib kits
    IF(${kit} MATCHES "^Contrib")
      STRING(REGEX REPLACE "^CONTRIB(.*)" "CONTRIB_\\1" kit_uc ${kit_uc})
    ENDIF(${kit} MATCHES "^Contrib")
    SET(VTKTUD_KITS_FORCONFIG ${VTKTUD_KITS_FORCONFIG} ${kit_uc})
  ENDFOREACH(kit)
  CONFIGURE_FILE( 
    ${VTKTUD_SOURCE_DIR}/UseVTKTUD.cmake.in
    ${VTKTUD_BINARY_DIR}/UseVTKTUD.cmake
    COPYONLY IMMEDIATE
  )
  CONFIGURE_FILE(
    ${VTKTUD_SOURCE_DIR}/VTKTUDConfig.cmake.in
    ${VTKTUD_BINARY_DIR}/VTKTUDConfig.cmake
    @ONLY IMMEDIATE
  )

  #
  # Generate python import module
  #

  IF(VTKTUD_WRAP_PYTHON)
    SET(VTKTUD_KITS_FORPYTHON ${VTKTUD_KITS_BUILT})
    CONFIGURE_FILE(
      "${VTKTUD_SOURCE_DIR}/Wrapping/Python/vtktud.py.in" 
      "${VTKTUD_BINARY_DIR}/Wrapping/Python/vtktud.py"
      @ONLY
    )
  ENDIF(VTKTUD_WRAP_PYTHON)

  #
  # Disable deprecation warnings for standard C and STL functions in VS2005
  # and later
  #

  IF(CMAKE_COMPILER_2005)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
    ADD_DEFINITIONS(-D_SCL_SECURE_NO_DEPRECATE)
  ENDIF(CMAKE_COMPILER_2005) 

ENDIF (VTKTUD_CAN_BUILD)
