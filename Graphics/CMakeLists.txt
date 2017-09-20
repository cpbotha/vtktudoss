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
TARGET_LINK_LIBRARIES(vtktudossGraphics ${VTK_LIBRARIES})

