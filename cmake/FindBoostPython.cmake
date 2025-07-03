#[========================================================================[.rst:
FindBoostPython
---------------

Finds the Boost Python library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

`BoostPython_FOUND`
  True if the system has the Boost Python library.
`BoostPython_INCLUDE_DIRS`
  Include directories needed to use Boost Python.
`BoostPython_LIBRARY`
  Path to Boost Python library.

#]========================================================================]

include(FindPackageHandleStandardArgs)

# If a conda environment is active, search for Boost Python within it.
if(DEFINED ENV{CONDA_PREFIX})
  # Conda Boost DOES ship with Boost Python,
  # but does not provide a CMake config for Boost Python.
  # Previously we could fall back from config mode to module mode FindBoost.
  # This was deprecated in CMake 4.0.
  # Our current best solution is to search for Boost Python manually.
  message(STATUS "Searching for Boost Python in "
    "active Anaconda environment ($ENV{CONDA_PREFIX})")
  # Search for Boost Python header files.
  find_path(BoostPython_INCLUDE_DIRS
    NAMES boost/python.hpp
    PATHS $ENV{CONDA_PREFIX}/include
    NO_DEFAULT_PATH)
  # Search for Boost Python library.
  file(GLOB BOOST_PYTHON_LIB_CANDIDATES
    "$ENV{CONDA_PREFIX}/lib/libboost_python*.so"    # GNU+Linux
    "$ENV{CONDA_PREFIX}/lib/libboost_python*.dylib" # macOS
    "$ENV{CONDA_PREFIX}/lib/*boost_python*.lib"     # Windows
  )
  list(LENGTH BOOST_PYTHON_LIB_CANDIDATES NUM_BOOST_PYTHON_LIB_CANDIDATES)
  if (NUM_BOOST_PYTHON_LIB_CANDIDATES EQUAL 1)
    # Set the output variables and print the status message.
    set(BoostPython_LIBRARY ${BOOST_PYTHON_LIB_CANDIDATES})
    find_package_handle_standard_args(BoostPython DEFAULT_MSG
      BoostPython_LIBRARY BoostPython_INCLUDE_DIRS)
    return()
  elseif (NUM_BOOST_PYTHON_LIB_CANDIDATES GREATER 1)
    message(FATAL_ERROR
      "Found multiple ambiguous Boost Python candidates: "
      "${BOOST_PYTHON_LIB_CANDIDATES}")
  endif()
endif()

# Search for Boost Python among the system libraries.
message(STATUS "Searching for Boost Python in system libraries")
find_package(Boost CONFIG COMPONENTS python)
# Set the output variables and print the status message.
if (TARGET Boost::python)
  set(BoostPython_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
  get_target_property(BoostPython_LIBRARY Boost::python LOCATION)
endif()
find_package_handle_standard_args(BoostPython DEFAULT_MSG
  BoostPython_LIBRARY BoostPython_INCLUDE_DIRS)