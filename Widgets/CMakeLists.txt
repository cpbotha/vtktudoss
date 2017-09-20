#
# Source files
#
# Here is where you can add the name of your local imaging classes.
#

SET (Widgets_SRCS
   vtkAdvancedAngleWidget.cxx
   vtkEllipseRepresentation.cxx
   vtkEllipseWidget.cxx
)

# --------------------------------------------------------------------------
# You probably do not need to modify anything below this line

# Create the vtktudossWidgets C++ library
ADD_LIBRARY (vtktudossWidgets ${Widgets_SRCS})
TARGET_LINK_LIBRARIES(vtktudossWidgets vtktudossGraphics ${VTK_LIBRARIES})

