#---------------------------------------------------------------------------
#
# BTKConfig.cmake - BTK CMake configuration file for external projects.
#
# This file is configured by BTK and used by the UseBTK.cmake module
# to load BTK's settings for an external project.

# Compute the installation prefix from this BTKConfig.cmake file location.
GET_FILENAME_COMPONENT(BTK_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(BTK_INSTALL_PREFIX "${BTK_INSTALL_PREFIX}" PATH)
GET_FILENAME_COMPONENT(BTK_INSTALL_PREFIX "${BTK_INSTALL_PREFIX}" PATH)

# The BTK include file directories.
SET(BTK_INCLUDE_DIRS "${BTK_INSTALL_PREFIX}/include/btk-0.1;${BTK_INSTALL_PREFIX}/include/btk-0.1/BasicFilters;${BTK_INSTALL_PREFIX}/include/btk-0.1/Common;${BTK_INSTALL_PREFIX}/include/btk-0.1/IO;${BTK_INSTALL_PREFIX}/include/btk-0.1/Utilities;${BTK_INSTALL_PREFIX}/include/btk-0.1/Utilities/eigen2")

# The BTK library directories.
SET(BTK_LIBRARY_DIRS "${BTK_INSTALL_PREFIX}/lib/btk-0.1")

# The C and C++ flags added by BTK to the cmake-configured flags.
SET(BTK_REQUIRED_C_FLAGS "")
SET(BTK_REQUIRED_CXX_FLAGS "")
SET(BTK_REQUIRED_LINK_FLAGS "")

# The BTK version number
SET(BTK_VERSION_MAJOR "0")
SET(BTK_VERSION_MINOR "1")
SET(BTK_VERSION_PATCH "10")

# The location of the UseBTK.cmake file.
SET(BTK_USE_FILE "${BTK_INSTALL_PREFIX}/share/btk-0.1/UseBTK.cmake")

# The build settings file.
SET(BTK_BUILD_SETTINGS_FILE "${BTK_INSTALL_PREFIX}/share/btk-0.1/BTKBuildSettings.cmake")

# The library dependencies file.
SET(BTK_LIBRARY_DEPENDS_FILE "${BTK_INSTALL_PREFIX}/share/btk-0.1/BTKLibraryDepends.cmake")

# Whether BTK was built with shared libraries.
SET(BTK_BUILD_SHARED "OFF")

# Whether BTK was built with Wrapping support.
SET(BTK_WRAP_MATLAB "ON")

# Whether BTK was built with Wrapping support.
SET(BTK_WRAP_SCILAB "OFF")

# A list of all libraries for BTK.  Those listed here should
# automatically pull in their dependencies.
SET(BTK_LIBRARIES BTKIO BTKBasicFilters)

# The BTK library dependencies.
IF(NOT BTK_NO_LIBRARY_DEPENDS AND EXISTS "${BTK_LIBRARY_DEPENDS_FILE}")
  INCLUDE("${BTK_LIBRARY_DEPENDS_FILE}")
ENDIF(NOT BTK_NO_LIBRARY_DEPENDS AND EXISTS "${BTK_LIBRARY_DEPENDS_FILE}")

# Whether BTK was built using a system eigen2
SET(BTK_USE_SYSTEM_EIGEN2 "OFF")

# The EIGEN2_DIR setting used to build BTK.  Set if BTK_USE_SYSTEM_EIGEN2 is true.
SET(BTK_EIGEN2_INCLUDE_DIR "")
