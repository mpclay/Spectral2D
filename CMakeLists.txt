# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>

# Main cmake setup for the project.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.9)
PROJECT(Spectral2D)
ENABLE_LANGUAGE(Fortran)

# Main source directory.
ADD_SUBDIRECTORY(src)

# Since FFTW doesn't have a FindFFTW.cmake, we include it in the code. The
# user must have the environment variable FFTW_ROOT set properly.
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Find the MPI libraries.
FIND_PACKAGE(MPI REQUIRED)
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
ADD_DEFINITIONS(${MPI_Fortran_COMPILE_FLAGS})

# Find the FFTW libraries.
SET(FFTW_USE_STATIC_LIBS TRUE)
SET(FFTW_FORTRAN_ENABLED TRUE)
SET(FFTW_USE_MPI TRUE)
FIND_PACKAGE(FFTW REQUIRED)
INCLUDE_DIRECTORIES(${FFTW_INCLUDES})

# Find the HDF5 library.
SET(HDF5_USE_STATIC_LIBRARIES FALSE)
FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR_FORTRAN})
LINK_DIRECTORIES(${HDF5_LIBRARY_DIRS})

# Build the main program.
ADD_EXECUTABLE(spectral2D
               ${CORE_SOURCE_FILES}
               ${SPECTRAL2D_SOURCE_FILES})
TARGET_LINK_LIBRARIES(spectral2D ${MPI_Fortran_LIBRARIES})
TARGET_LINK_LIBRARIES(spectral2D ${FFTW_LIBRARIES})
TARGET_LINK_LIBRARIES(spectral2D ${HDF5_LIBRARIES})
INSTALL(TARGETS spectral2D DESTINATION ${CMAKE_INSTALL_PREFIX})

