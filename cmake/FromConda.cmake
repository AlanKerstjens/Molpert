#[========================================================================[.rst:
FromConda
---------

Checks if a list of input variables resolve to Anaconda environment paths.

Parameters
^^^^^^^^^^

`OUTPUT_FROM_CONDA`
  Name of the variable storing the function result.
`INPUT_PATH_VARS`
  (Implicit) List of evaluated variable names. May be of variable size.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

`OUTPUT_FROM_CONDA`
  True if all input variables resolve to Anaconda environment paths.
  False otherwise.

#]========================================================================]

function(from_conda OUTPUT_FROM_CONDA)
  # If we got no input (other than the required output variable) return false.
  if (${ARGC} EQUAL 1)
    set(${OUTPUT_FROM_CONDA} FALSE PARENT_SCOPE)
    return()
  endif()
  # If no conda environment is active return false.
  if (NOT DEFINED ENV{CONDA_PREFIX})
    set(${OUTPUT_FROM_CONDA} FALSE PARENT_SCOPE)
    return()
  endif()
  # Iterate over the (non-empty) input path variables.
  foreach(INPUT_PATH_VAR IN LISTS ARGN)
    set(INPUT_PATH ${${INPUT_PATH_VAR}})
    if(NOT INPUT_PATH STREQUAL "")
      # Search for the conda prefix in the input path.
      string(FIND "${INPUT_PATH}" "$ENV{CONDA_PREFIX}" POSITION)
      # If it isn't present the input path is not from a conda environment.
      # Return false.
      if (POSITION EQUAL -1)
        set(${OUTPUT_FROM_CONDA} FALSE PARENT_SCOPE)
        return()
      endif()
    endif()
  endforeach()
  # If all input paths are from a conda environment return true.
  set(${OUTPUT_FROM_CONDA} TRUE PARENT_SCOPE)
endfunction()