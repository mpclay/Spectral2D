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
!! FFTW and Fortran procedures.
MODULE Alloc_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPC, IWPF, RWPC, RWPF, CWPC

   IMPLICIT NONE

   ! Main interface to allocation routines.
   INTERFACE Alloc
      MODULE PROCEDURE AllocComplexFFTW
      MODULE PROCEDURE Alloc2DComplexFortran
   END INTERFACE

   ! Module procedures.
   PUBLIC :: GetAllocSize, Dealloc
   PRIVATE :: AllocComplexFFTW, Alloc2DComplexFortran

CONTAINS

   !> Procedure to determine size of FFTW working arrays.
   !!
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] kxG Global array of wavenumbers for the x direction.
   !> @param[in] kyG Global array of wavenumbers for the y direction.
   !> @param[in] kxMU Maximum useful x wavenumber.
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[out] nxP Number of x grid points on this process.
   !> @param[out] nyP Number of y grid points on this process.
   !> @param[out] j1 Starting j index for this process.
   !> @param[out] j2 Ending j index for this process.
   !> @param[out] localJ Offset for j based on FFTW MPI decomposition.
   !> @param[out] allocLocal Size of data array required by FFTW.
   !> @param[out] kxLP Lowest x wavenumber for this process.
   !> @param[out] kyLP Lowest y wavenumber for this process.
   !> @param[out] kxMP Largest x wavenumber for this process.
   !> @param[out] kyMP Largest y wavenumber for this process.
   !> @param[out] kxP Array of x wavenumbers for this process.
   !> @param[out] kyP Array of y wavenumbers for this process.
   SUBROUTINE GetAllocSize(nxG, nyG, kxG, kyG, kxMU, cISize, nxP, nyP, j1, j2, &
                           localJ, allocLocal, kxLP, kyLP, kxMP, kyMP, kxP, kyP)
      ! Required modules.
      USE ISO_C_BINDING
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, kxMU, cISize
      INTEGER(KIND=IWPF),DIMENSION(:),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      INTEGER(KIND=IWPF),INTENT(OUT) :: nxP, nyP, j1, j2, localJ, allocLocal
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLP, kyLP, kxMP, kyMP
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: kxP, kyP
      ! Local variables.
      ! Size of complex data arrays required by FFTW.
      INTEGER(KIND=IWPC) :: allocLocal_
      ! Local starting index for this task's portion of the data.
      INTEGER(KIND=IWPC) :: localJ_
      ! Number of columns (starting from localJ) this process owns.
      INTEGER(KIND=IWPC) :: jSize_
      ! Variable for locally accessing the kyP array.
      INTEGER(KIND=IWPF) :: jInd
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j

      ! Determine number of columns this MPI task gets for the FFTs. The
      ! dimensions are reversed because a C routine is being called.
      allocLocal_ = FFTW_MPI_LOCAL_SIZE_2D(INT(nyG, IWPC), &
                                           INT(cISize, IWPC), &
                                           MPI_COMM_WORLD, &
                                           jSize_, localJ_)
      allocLocal = INT(allocLocal_, IWPF)
      !
      ! The j array bounds for this process.
      j1 = INT(localJ_, IWPF) + 1_IWPF
      j2 = INT(localJ_ + jSize_, IWPF)
      !
      ! Set the size of grid points for this process.
      nxP = nxG
      nyP = INT(jSize_, IWPF)
      localJ = INT(localJ_, IWPF)
      !
      ! Set the min/max wavenumbers in each direction for this process.
      kxLP = 0_IWPF
      kxMP = kxMU
      kyLP = MINVAL(kyG(j1:j2))
      kyMP = MAXVAL(kyG(j1:j2))
      !
      ! Set the wavenumber arrays for this process.
      ALLOCATE(kxP(SIZE(kxG)))
      DO i = 1, SIZE(kxG)
         kxP(i) = kxG(i)
      END DO
      ALLOCATE(kyP(nyP))
      jInd = 1_IWPF
      DO j = j1, j2
         kyP(jInd) = kyG(j)
         jInd = jInd + 1_IWPF
      END DO
   END SUBROUTINE GetAllocSize

   !> Procedure to allocate working arrays with FFTW.
   !!
   !! Keep in mind that allocLocal might be larger than the acutal useful
   !! memory if FFTW MPI routines need additional space for transpositions.
   !! The useful portions of the array are at the front, and the rest of the
   !! array can be ignored for all practicaly purposes.
   !!
   !> @param[in] allocLocal Size of the data to be allocated.
   !> @param[in,out] array Data array for the simulation.
   SUBROUTINE AllocComplexFFTW(allocLocal, array)
      ! Required modules.
      USE ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: allocLocal
      TYPE(C_PTR),INTENT(INOUT) :: array

      ! Main memory allocation.
      array = FFTW_ALLOC_COMPLEX(INT(allocLocal, C_SIZE_T))
   END SUBROUTINE AllocComplexFFTW

   !> Procedure to allocate 2D working arrays.
   !!
   !> @param[in] cISize Desired size of the i-dimension for complex numbers.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] iStart Starting index for i dimension.
   !> @param[in] jStart Starting index for j dimension.
   !> @param[in,out] array Data array being allocated.
   SUBROUTINE Alloc2DComplexFortran(cISize, nyP, iStart, jStart, array)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyP, iStart, jStart
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
      jEnd = jStart + nyP - 1_IWPF
      ALLOCATE(array(iStart:iEnd,jStart:jEnd))
   END SUBROUTINE Alloc2DComplexFortran

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

