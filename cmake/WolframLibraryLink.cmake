include_guard(DIRECTORY)

set(
  WOLFRAM_LIBRARYLINK_INCLUDE_DIR
  ""
  CACHE PATH
  "Optional override for the directory containing WolframLibrary.h"
)

if(WOLFRAM_LIBRARYLINK_INCLUDE_DIR)
  set(_wolfram_librarylink_include_dir "${WOLFRAM_LIBRARYLINK_INCLUDE_DIR}")
else()
  set(_wolfram_librarylink_include_dir "")

  foreach(_wolfram_executable_name IN ITEMS WolframKernel MathKernel wolframscript)
    unset(_wolfram_executable CACHE)
    unset(_wolfram_executable)
    find_program(_wolfram_executable NAMES "${_wolfram_executable_name}")

    if(_wolfram_executable)
      get_filename_component(_wolfram_executable_real "${_wolfram_executable}" REALPATH)
      get_filename_component(_wolfram_executables_dir "${_wolfram_executable_real}" DIRECTORY)
      get_filename_component(_wolfram_installation_dir "${_wolfram_executables_dir}" DIRECTORY)
      set(_wolfram_candidate "${_wolfram_installation_dir}/SystemFiles/IncludeFiles/C")

      if(EXISTS "${_wolfram_candidate}/WolframLibrary.h")
        set(_wolfram_librarylink_include_dir "${_wolfram_candidate}")
        break()
      endif()
    endif()
  endforeach()
endif()

if(NOT EXISTS "${_wolfram_librarylink_include_dir}/WolframLibrary.h")
  message(
    FATAL_ERROR
    "Could not locate Wolfram LibraryLink headers. Ensure WolframKernel is on PATH "
    "or configure with -DWOLFRAM_LIBRARYLINK_INCLUDE_DIR=/path/to/SystemFiles/IncludeFiles/C."
  )
endif()

# Keep automatically detected paths out of the cache so that a later Wolfram
# upgrade is picked up when CMake is configured again.
set(WOLFRAM_LIBRARYLINK_INCLUDE_DIR "${_wolfram_librarylink_include_dir}")
message(STATUS "Wolfram LibraryLink headers: ${WOLFRAM_LIBRARYLINK_INCLUDE_DIR}")
