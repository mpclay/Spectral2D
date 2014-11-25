! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! Copyright (C) 2014 by Matthew Clay <mpclay@gmail.com>
!
!> @file Alloc.F90
!> @author Matthew Clay
!> @brief Module to handle memory allocation/deallocation for the program.
!!
!! This module provides routines to allocate and free memory using underlying
!! FFTW procedures.
MODULE Alloc_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPC, IWPF, RWPC, RWPF, CWPC

   IMPLICIT NONE

   ! Main interface to allocation routines.
   INTERFACE Alloc
      MODULE PROCEDURE Alloc2DComplexFFTW
      MODULE PROCEDURE Alloc2DComplexFortran
      MODULE PROCEDURE Alloc4DComplexFFTW
   END INTERFACE

   ! Module procedures.
   PUBLIC :: GetAllocSize, Dealloc
   PRIVATE :: Alloc2DComplexFFTW, Alloc2DComplexFortran, Alloc4DComplexFFTW

CONTAINS

   !> Procedure to determine size of FFTW working arrays.
   !!
   !> @param[in] nxW Working number of grid points in the x direction.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[out] j1 Starting j index for this process.
   !> @param[out] j2 Ending j index for this process.
   !> @param[out] jSize Number of j columns owned by this MPI process.
   !> @param[out] localJ Offset for j based on FFTW MPI decomposition.
   !> @param[out] rISize Size of i-dimension for real data array.
   !> @param[out] cISize Size of i-dimension for complex data array.
   SUBROUTINE GetAllocSize(nxW, nyW, j1, j2, jSize, localJ, rISize, cISize)
      ! Required modules.
      USE ISO_C_BINDING
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxW, nyW
      INTEGER(KIND=IWPF),INTENT(OUT) :: j1, j2, jSize, localJ, rISize, cISize
      ! Local variables.
      ! Size of working array for this process as determined by FFTW.
      INTEGER(KIND=IWPC) :: allocLocal
      ! Local starting index for this task's portion of the data.
      INTEGER(KIND=IWPC) :: localJ_
      ! Number of columns (starting from localJ) this process owns.
      INTEGER(KIND=IWPC) :: jSize_

      ! Determine number of columns this MPI task gets for the FFTs. The
      ! dimensions are reversed because a C routine is being called.
      allocLocal = FFTW_MPI_LOCAL_SIZE_2D(INT(nyW, IWPC), &
                                          INT(nxW, IWPC), &
                                          MPI_COMM_WORLD, &
                                          jSize_, localJ_)
      !
      ! The j array bounds for this process.
      j1 = INT(localJ_, IWPF) + 1_IWPF
      j2 = INT(localJ_ + jSize_, IWPF)
      !
      ! Determine the size of the arrays for in-place transforms.
      rISize = 2_IWPF*(nxW/2_IWPF + 1_IWPF)
      cISize = nxW/2_IWPF + 1_IWPF

      ! Subroutine outputs back to the calling code in the proper type.
      jSize = INT(jSize_, IWPF)
      localJ = INT(localJ_, IWPF)
   END SUBROUTINE GetAllocSize

   !> Procedure to allocate 2D working arrays.
   !!
   !> @param[in] cISize Desired size of the i-dimension for complex numbers.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[in,out] array Data array for the simulation.
   SUBROUTINE Alloc2DComplexFFTW(cISize, nyW, array)
      ! Required modules.
      USE ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyW
      TYPE(C_PTR),INTENT(INOUT) :: array

      ! Main memory allocation.
      !
      ! Each array is allocated by FFTW for complex data. Each array has
      ! dimensions:
      !
      !     (cISize=size of the i dimension in spectral space) X
      !     (nyW=size of the j dimension in physical/spectral space)
      array = FFTW_ALLOC_COMPLEX(INT(cISize*nyW, C_SIZE_T))
   END SUBROUTINE Alloc2DComplexFFTW

   !> Procedure to allocate 2D working arrays.
   !!
   !> @param[in] cISize Desired size of the i-dimension for complex numbers.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[in] iStart Starting index for i dimension.
   !> @param[in] jStart Starting index for j dimension.
   !> @param[in,out] array Data array being allocated.
   SUBROUTINE Alloc2DComplexFortran(cISize, nyW, iStart, jStart, array)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyW, iStart, jStart
      COMPLEX(KIND=CWPC),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) :: array
      ! Local variables.
      ! Ending indices for the allocated array.
      INTEGER(KIND=IWPF) :: iEnd, jEnd

      ! Free memory if it is already in use.
      IF (ALLOCATED(array)) THEN
         DEALLOCATE(array)
      END IF

      ! Allocate memory with desired starting and ending indices.
      iEnd = iStart + cISize - 1_IWPF
      jEnd = jStart + nyW - 1_IWPF
      ALLOCATE(array(iStart:iEnd,jStart:jEnd))
   END SUBROUTINE Alloc2DComplexFortran

   !> Procedure to allocate 4D working arrays.
   !!
   !> @param[in] cISize Desired size of the i-dimension for complex numbers.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[in] nVar Number of variables in the simulation.
   !> @param[in] nStr Required number of storage arrays for time integration.
   !> @param[in,out] array Data array for the simulation.
   SUBROUTINE Alloc4DComplexFFTW(cISize, nyW, nVar, nStr, array)
      ! Required modules.
      USE ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyW, nVar, nStr
      TYPE(C_PTR),INTENT(INOUT) :: array

      ! Main memory allocation.
      !
      ! Each array is allocated by FFTW for complex data. Each array has
      ! dimensions:
      !
      !     (cISize=size of the i dimension in spectral space) X
      !     (nyW=size of the j dimension in physical/spectral space) X
      !     (nVar=number of variables in the simulation) X
      !     (nStr=number of storage sites for time integration).
      array = FFTW_ALLOC_COMPLEX(INT(cISize*nyW*nVar*nStr, C_SIZE_T))
   END SUBROUTINE Alloc4DComplexFFTW

   !> Routine to free memory for the simulation.
   !!
   !> @param[in,out] array Working array in the simulation.
   SUBROUTINE Dealloc(array)
      ! Required modules.
      USE ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      TYPE(C_PTR),INTENT(INOUT) :: array

      ! Free memory using FFTW routines.
      CALL FFTW_FREE(array)
   END SUBROUTINE Dealloc

END MODULE Alloc_m

