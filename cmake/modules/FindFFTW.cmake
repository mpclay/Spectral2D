# Find the FFTW library.
#
# Usage:
#   FIND_PACKAGE(FFTW [REQUIRED] [QUIET])
#
# It sets the following variables:
#   FFTW_FOUND              ... true if fftw is found on the system
#   FFTW_LIBRARIES          ... full path to fftw library
#   FFTW_INCLUDES           ... fftw include directory
#
# The following variables will be checked by the function
#   FFTW_FORTRAN_ENABLED    ... if true, Fortran includes will be looked for
#   FFTW_USE_STATIC_LIBS    ... if true, only static libraries are found
#   FFTW_USE_MPI            ... if true, the MPI versions of FFTW are located
#   FFTW_SINGLE_PRECISION   ... if true, the floating point libraries will be
#                               search for
#   FFTW_QUAD_PRECISION     ... if true, the quadruple precision libraries will
#                               be searched for
#   FFTW_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   FFTW_LIBRARY            ... fftw library to use
#   FFTW_INCLUDE_DIR        ... fftw include directory
#

# The FFTW_ROOT environment variable is recommended to be set.
SET(FFTW_ROOT $ENV{FFTW_ROOT})

# Warn the user if FFTW_ROOT is not found.
IF(NOT FFTW_ROOT)
  MESSAGE("FFTW_ROOT environment variable not set. Using PkgConfig.")
ENDIF()

# Check if we can use PkgConfig.
FIND_PACKAGE(PkgConfig)

# Determine from PKG (this is a backup).
IF(PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
  pkg_check_modules(PKG_FFTW QUIET "fftw3")
ENDIF()

# Check whether to search static or dynamic libs.
SET(CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES})
IF(${FFTW_USE_STATIC_LIBS})
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
ELSE()
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDIF()

# Locate the FFTW libraries and includes with the FFTW_ROOT path.

# Double precision version of FFTW.
FIND_LIBRARY(FFTW_LIB
             NAMES "fftw3"
             PATHS ${FFTW_ROOT}
             PATH_SUFFIXES "lib" "lib64"
             NO_DEFAULT_PATH)

IF(FFTW_SINGLE_PRECISION)
  # Single precision version of FFTW.
  FIND_LIBRARY(FFTWF_LIB
               NAMES "fftw3f"
               PATHS ${FFTW_ROOT}
               PATH_SUFFIXES "lib" "lib64"
               NO_DEFAULT_PATH)
ENDIF()

IF(FFTW_QUAD_PRECISION)
  # Quadruple precision version of FFTW.
  FIND_LIBRARY(FFTWL_LIB
               NAMES "fftw3l"
               PATHS ${FFTW_ROOT}
               PATH_SUFFIXES "lib" "lib64"
               NO_DEFAULT_PATH)
ENDIF()

# Find the includes.
FIND_PATH(FFTW_INCLUDES
          NAMES "fftw3.h"
          PATHS ${FFTW_ROOT}
          PATH_SUFFIXES "include"
          NO_DEFAULT_PATH)

# Find the includes for Fortran.
IF(FFTW_FORTRAN_ENABLED)
  FIND_PATH(FFTW_Fortran_INCLUDES
            NAMES "fftw3.f03"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH)
ENDIF()

# MPI parallelized version of FFTW.
IF(FFTW_USE_MPI)
  # Double precision version of FFTW MPI.
  FIND_LIBRARY(FFTW_MPI_LIB
               NAMES "fftw3_mpi"
               PATHS ${FFTW_ROOT}
               PATH_SUFFIXES "lib" "lib64"
               NO_DEFAULT_PATH)

  IF(FFTW_SINGLE_PRECISION)
    # Single precision version of FFTW MPI.
    FIND_LIBRARY(FFTWF_MPI_LIB
                 NAMES "fftw3f_mpi"
                 PATHS ${FFTW_ROOT}
                 PATH_SUFFIXES "lib" "lib64"
                 NO_DEFAULT_PATH)
  ENDIF()

  IF(FFTW_QUAD_PRECISION)
    # Quadruple precision version of FFTW MPI.
    FIND_LIBRARY(FFTWL_MPI_LIB
                 NAMES "fftw3l_mpi"
                 PATHS ${FFTW_ROOT}
                 PATH_SUFFIXES "lib" "lib64"
                 NO_DEFAULT_PATH)
  ENDIF()

  # Find the includes.
  FIND_PATH(FFTW_MPI_INCLUDES
            NAMES "fftw3-mpi.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH)

  # Find the includes for Fortran.
  IF(FFTW_FORTRAN_ENABLED)
    FIND_PATH(FFTW_MPI_Fortran_INCLUDES
              NAMES "fftw3.f03"
              PATHS ${FFTW_ROOT}
              PATH_SUFFIXES "include"
              NO_DEFAULT_PATH)
  ENDIF()
ENDIF()

# Set the main output variables to used by cmake.
SET(FFTW_LIBRARIES
    ${FFTW_MPI_LIB}
    ${FFTWF_MPI_LIB}
    ${FFTWL_MPI_LIB}
    ${FFTW_LIB}
    ${FFTWF_LIB}
    ${FFTWL_LIB})
SET(FFTW_INCLUDES
    ${FFTW_MPI_INCLUDES}
    ${FFTW_MPI_Fortran_INCLUDES}
    ${FFTW_INCLUDES}
    ${FFTW_Fortran_INCLUDES})

MARK_AS_ADVANCED(FFTW_INCLUDES FFTW_LIBRARIES)

