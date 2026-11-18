### Dependent
# Is C needed?
find_package(HDF5 COMPONENTS C Fortran REQUIRED)

### LM4 Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Waliasing -fcray-pointer -fconvert=big-endian -ffree-line-length-none -fno-range-check -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -sox")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv")
  set(CMAKE_Fortran_LINK_FLAGS "")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${LM4_Extra_FORTRAN_FLAGS}")

# Add preprocessor flag for namelist
add_definitions(-DINTERNAL_FILE_NML)

# Set MPI flags
if(MPI)
  add_definitions(-DMPI)
endif()

# LM4 src files
include("lm4_src_files.cmake")

# gfortran's preprocessor doesn't support some of the C preprocessing of LM4 macros
# Run the sources that use such macros (directly, or through an include file) through
# a C preprocessor first
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(_ansi_macro_regex "^#[ \t]*define[^\n]*([^#]#[A-Za-z_]|##)")

  # include files defining such macros
  file(GLOB_RECURSE _lm4_headers ${CMAKE_CURRENT_SOURCE_DIR}/*.inc ${CMAKE_CURRENT_SOURCE_DIR}/*.h)
  set(_ansi_headers)
  foreach(_hdr IN LISTS _lm4_headers)
    file(STRINGS ${_hdr} _hits REGEX "${_ansi_macro_regex}")
    if(_hits)
      get_filename_component(_name ${_hdr} NAME)
      list(APPEND _ansi_headers ${_name})
    endif()
  endforeach()

  get_directory_property(_lm4_cpp_defs COMPILE_DEFINITIONS)
  list(JOIN _lm4_cpp_defs "|" _lm4_cpp_defs)

  foreach(_src IN LISTS lm4_src_files)
    if(NOT _src MATCHES "\\.F90$")
      continue()
    endif()
    set(_needs_ansi_cpp FALSE)
    file(STRINGS ${_src} _hits REGEX "${_ansi_macro_regex}")
    if(_hits)
      set(_needs_ansi_cpp TRUE)
    else()
      file(STRINGS ${_src} _includes REGEX "^#[ \t]*include")
      foreach(_inc IN LISTS _includes)
        foreach(_name IN LISTS _ansi_headers)
          if(_inc MATCHES "[\"</]${_name}[\">]")
            set(_needs_ansi_cpp TRUE)
          endif()
        endforeach()
      endforeach()
    endif()
    if(NOT _needs_ansi_cpp)
      continue()
    endif()

    string(REGEX REPLACE "\\.F90$" ".f90" _out ${CMAKE_CURRENT_BINARY_DIR}/gnu_fpp/${_src})
    add_custom_command(
      OUTPUT  ${_out}
      COMMAND ${CMAKE_COMMAND}
              -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/${_src}
              -DOUT=${_out}
              -DROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DCPP=${CMAKE_C_COMPILER}
              -DCPP_DEFS=${_lm4_cpp_defs}
              "-DCPP_INCS=$<JOIN:$<TARGET_PROPERTY:lm4,INCLUDE_DIRECTORIES>,|>"
              -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lm4_gnu_fpp.cmake
      DEPENDS ${_src} ${_lm4_headers} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/lm4_gnu_fpp.cmake
      COMMENT "Preprocessing ${_src} with ANSI cpp for gfortran"
      VERBATIM)
    set_source_files_properties(${_out} PROPERTIES Fortran_PREPROCESS OFF)
    list(TRANSFORM lm4_src_files REPLACE "^${_src}$" ${_out})
  endforeach()
endif()

add_library(lm4 OBJECT ${lm4_src_files})
set_target_properties(lm4 PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
#target_compile_definitions(lm4 PUBLIC "${_lm4_defs}")
target_include_directories(lm4 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)

find_package(HDF5 COMPONENTS C Fortran REQUIRED)
target_link_libraries(lm4 PUBLIC ESMF::ESMF fms ${HDF5_Fortran_LIBRARIES})

###############################################################################
### Install
###############################################################################

install(
  TARGETS lm4
  EXPORT  lm4-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/LM4)

install(EXPORT      lm4-config
        DESTINATION lib/cmake)
