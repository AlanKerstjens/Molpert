cmake_minimum_required(VERSION 3.18)

# Define the project.
project(Molpert
  LANGUAGES CXX
  DESCRIPTION "Library for graph-based perturbation of RDKit molecules")

# Import project's CMake modules.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(FromConda)

# Specify the C++ standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define the directories where the source code can be found and where the
# compiled libraries/executables should be stored.
set(Molpert_SOURCE_DIR ${PROJECT_SOURCE_DIR}/source)
set(Molpert_LIBRARY_DIR ${PROJECT_SOURCE_DIR}/lib)

# Search for the RDKit.
# If found this will also retrieve the Boost version used to build the RDKit.
find_package(RDKit REQUIRED)

# Define the library target.
add_library(Molpert INTERFACE)

target_include_directories(Molpert INTERFACE
  ${Molpert_SOURCE_DIR}
  ${Molpert_SOURCE_DIR}/addons
  ${Boost_INCLUDE_DIRS}
  ${RDKit_INCLUDE_DIRS})

target_link_directories(Molpert INTERFACE
  ${Boost_LIBRARY_DIRS}
  ${RDKit_LIBRARY_DIRS})

target_link_libraries(Molpert INTERFACE
  RDKitRDGeneral
  RDKitGraphMol
  RDKitFingerprints
  RDKitSmilesParse)


# Enable building Python bindings.
option(BUILD_PYTHON_BINDINGS "Enable building Python bindings" ON)

# Conditionally build Python bindings.
if(BUILD_PYTHON_BINDINGS)

  # Add preprocessor definitions for Python binding compilation.
  add_definitions(-DBUILD_PYTHON_BINDINGS)

  # Search for the Python module dependencies.
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
  find_package(BoostPython REQUIRED)

  # To avoid ABI incompatibility, Molpert needs to be built using the 
  # same libraries that were used to build the RDKit. Check if we found them.
  from_conda(RDKit_FROM_CONDA
    RDKit_INCLUDE_DIRS RDKit_LIBRARY_DIRS)
  if (${RDKit_FROM_CONDA})
    from_conda(REST_FROM_CONDA
      Boost_INCLUDE_DIRS Boost_LIBRARY_DIRS
      BoostPython_INCLUDE_DIRS BoostPython_LIBRARY
      Python3_INCLUDE_DIRS Python3_LIBRARY_DIRS Python3_LIBRARIES)
    if (NOT REST_FROM_CONDA)
      message(SEND_ERROR 
        "The RDKit was installed with conda, "
        "but some of the remaining dependencies were not. "
        "This is very likely to result in ABI incompatibility. "
        "Install the remaining dependencies through conda, "
        "or install the RDKit from source.")
    endif()
  endif()

  # Define the Python module target.
  add_library(pyMolpert MODULE
    ${Molpert_SOURCE_DIR}/wrap/pyMolpert.cpp)

  target_include_directories(pyMolpert PUBLIC
    ${Molpert_SOURCE_DIR}/wrap
    ${Python3_INCLUDE_DIRS}
    ${BoostPython_INCLUDE_DIRS})

  target_link_directories(pyMolpert PUBLIC
    ${Python3_LIBRARY_DIRS})

  target_link_libraries(pyMolpert PUBLIC
    Molpert
    ${Python3_LIBRARIES}
    ${BoostPython_LIBRARY})

  # Define and set the RPATH.
  set(LINK_DIRECTORIES
    ${Boost_LIBRARY_DIRS}
    ${RDKit_LIBRARY_DIRS}
    ${Python3_LIBRARY_DIRS})
  list(REMOVE_DUPLICATES LINK_DIRECTORIES)
  string(JOIN ":" RPATH ${LINK_DIRECTORIES})

  set_target_properties(pyMolpert PROPERTIES
    INSTALL_RPATH ${RPATH})

  # Specify the Python module name and installation location.
  set_target_properties(pyMolpert PROPERTIES
    PREFIX ""
    IMPORT_PREFIX ""
    OUTPUT_NAME "molpert")

  install(
    TARGETS
      pyMolpert
    DESTINATION ${Molpert_LIBRARY_DIR})

endif()