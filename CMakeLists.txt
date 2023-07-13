cmake_minimum_required(VERSION 3.18)

# Define the project.
project(Molpert
  LANGUAGES CXX
  DESCRIPTION "Library for graph-based perturbation of RDKit molecules")

# Specify the C++ standard.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define the directories where the source code can be found and where the
# compiled libraries/executables should be stored.
set(MOLPERT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/source)
set(MOLPERT_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

# Include the library directory in the default RPATH of the targets.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${MOLPERT_LIB_DIR})

# Search for the RDKit.
# If found this will also retrieve the Boost version used to build the RDKit.
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(RDKit REQUIRED)

# Define the library targets.
add_library(Molpert INTERFACE)

target_include_directories(Molpert INTERFACE
  ${MOLPERT_SOURCE_DIR}
  ${MOLPERT_SOURCE_DIR}/addons
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

# Enable conditional compilation of Python bindings
option(BUILD_PYTHON_BINDINGS "Enable Python binding compilation" ON)
if(BUILD_PYTHON_BINDINGS)

  add_definitions(-DBUILD_PYTHON_BINDINGS)

  find_package(Python REQUIRED COMPONENTS Interpreter Development)
  find_package(Boost REQUIRED COMPONENTS python)

  add_library(pyMolpert MODULE
    ${MOLPERT_SOURCE_DIR}/wrap/pyMolpert.cpp)

  target_include_directories(pyMolpert PUBLIC
    ${MOLPERT_SOURCE_DIR}/wrap
    ${Python_INCLUDE_DIRS})

  target_link_directories(pyMolpert PUBLIC
    ${Python_LIBRARY_DIRS})

  target_link_libraries(pyMolpert PUBLIC
    ${Python_LIBRARIES}
    Boost::python
    Molpert)

  set_target_properties(pyMolpert PROPERTIES
    PREFIX ""
    IMPORT_PREFIX ""
    OUTPUT_NAME "molpert")

  install(
    TARGETS
      pyMolpert
    DESTINATION ${MOLPERT_LIB_DIR})

endif()