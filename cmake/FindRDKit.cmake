#[=======================================================================[.rst:
FindRDKit
-------

Finds the RDKit library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``RDKit_FOUND``
  True if the system has the RDKit library.
``RDKit_VERSION``
  The version of the RDKit library which was found.
``RDKit_INCLUDE_DIRS``
  Include directories needed to use RDKit.
``RDKit_LIBRARY_DIRS``
  Library directories needed to link to RDKit.

#]=======================================================================]

set(RDKit_FOUND FALSE)

if(RDKit_ROOT)
  # Recursively search for rdkit-config.cmake and rdkit-config-version.cmake.
  # The search is started in RDKit_ROOT. If successful find_package sets 
  # RDKit_FOUND, RDKit_DIR and RDKit_INCLUDE_DIRS.
  find_package(RDKit CONFIG 
    PATHS ${RDKit_ROOT}
    NO_DEFAULT_PATH)
endif()

# If an Anaconda environment is active check if the RDKit is installed there.
if(NOT RDKit_FOUND AND DEFINED ENV{CONDA_PREFIX})
  find_package(RDKit CONFIG 
    PATHS $ENV{CONDA_PREFIX}
    NO_DEFAULT_PATH)
  if(RDKit_FOUND)
    set(RDKit_ROOT $ENV{CONDA_PREFIX})
  endif()
endif()

# If the RDBASE environment variable is set search in said directory.
if(NOT RDKit_FOUND AND DEFINED ENV{RDBASE})
  find_package(RDKit CONFIG
    PATHS $ENV{RDBASE}
    NO_DEFAULT_PATH)
  if(RDKit_FOUND)
    set(RDKit_ROOT $ENV{RDBASE})
  endif()
endif()

if(RDKit_FOUND)
  # The RDKit config doesn't set RDKit_LIBRARY_DIRS so we set it manually. 
  # We cache it because RDKit_DIR (and in turn RDKit_INCLUDE_DIRS) is cached.
  set(RDKit_LIBRARY_DIRS ${RDKit_ROOT}/lib 
    CACHE PATH "Directory containing RDKit libraries")
  message(STATUS "Found RDKit: "
    "${RDKit_INCLUDE_DIRS}, " 
    "${RDKit_LIBRARY_DIRS} "
    "(version ${RDKit_VERSION})")
elseif(RDKit_FIND_REQUIRED)
  message(FATAL_ERROR "RDKit not found. "
    "Consider specifying RDKit_ROOT.")
endif()