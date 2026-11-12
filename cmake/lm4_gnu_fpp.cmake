# Preprocess a LM4 source file with an ANSI C preprocessor for gfortran.
#
# cpp would read the Fortran concatenation operator // as a C++ comment, so the
# source file and all include files are copied to a seperate work tree with //
# replaced by a placeholder, preprocessed there, and the placeholder is restored.
#
# Usage:
#   cmake -DSRC=<source> -DOUT=<output .f90> -DROOT=<LM4-driver source dir>
#         -DCPP=<C compiler> [-DCPP_DEFS=a|b=1] [-DCPP_INCS=/dir1|/dir2]
#         -P lm4_gnu_fpp.cmake

cmake_minimum_required(VERSION 3.19)

foreach(_var SRC OUT ROOT CPP)
  if(NOT DEFINED ${_var})
    message(FATAL_ERROR "lm4_gnu_fpp.cmake: ${_var} is not set")
  endif()
endforeach()

set(_marker "__LM4_FORTRAN_CONCAT_OP__")
set(_work "${OUT}.work")

# Mirror the source file and all include files with // protected, keeping the
# directory layout so that relative includes ("../shared/debug.inc") resolve.
file(RELATIVE_PATH _rel "${ROOT}" "${SRC}")
file(GLOB_RECURSE _incs RELATIVE "${ROOT}" "${ROOT}/*.inc" "${ROOT}/*.h")
file(REMOVE_RECURSE "${_work}")
foreach(_f IN LISTS _rel _incs)
  file(READ "${ROOT}/${_f}" _txt)
  string(REPLACE "//" "${_marker}" _txt "${_txt}")
  file(WRITE "${_work}/${_f}" "${_txt}")
endforeach()

set(_args -E -undef -x assembler-with-cpp)
if(CPP_DEFS)
  string(REPLACE "|" ";" _defs "${CPP_DEFS}")
  foreach(_d IN LISTS _defs)
    list(APPEND _args "-D${_d}")
  endforeach()
endif()
if(CPP_INCS)
  string(REPLACE "|" ";" _incdirs "${CPP_INCS}")
  foreach(_i IN LISTS _incdirs)
    list(APPEND _args "-I${_i}")
  endforeach()
endif()

execute_process(
  COMMAND "${CPP}" ${_args} "${_work}/${_rel}"
  OUTPUT_VARIABLE _out
  ERROR_VARIABLE  _err
  RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "lm4_gnu_fpp.cmake: preprocessing ${SRC} failed:\n${_err}")
endif()

# Restore //, and point line markers back at the real source files so that
# compiler diagnostics refer to the original files.
string(REPLACE "${_marker}" "//" _out "${_out}")
string(REPLACE "${_work}/" "${ROOT}/" _out "${_out}")

# Only touch the output when it changes, to avoid needless recompilation.
set(_old "")
if(EXISTS "${OUT}")
  file(READ "${OUT}" _old)
endif()
if(NOT _old STREQUAL _out)
  file(WRITE "${OUT}" "${_out}")
endif()
file(REMOVE_RECURSE "${_work}")
