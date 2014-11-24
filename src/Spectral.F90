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
!> @file Spectral.F90
!> @author Matthew Clay
!> @brief Module to solve spatial part of the governing equations.
MODULE Spectral_m

   ! Required modules.
   USE ISO_C_BINDING
   USE Parameters_m,ONLY: IWPF, RWPF, IWPC, RWPC

   IMPLICIT NONE

   ! FFTW procedure definitions.
   INCLUDE 'fftw3-mpi.f03'

   ! Methods available to handle aliasing.
   !
   !> No dealiasing at all.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: NO_DEALIASING = 0_IWPF
   !> Set size of grid based on 3/2 rule.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: DEALIAS_3_2_RULE = 1_IWPF
   !> Set size of grid based on 2/3 rule.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: DEALIAS_2_3_RULE = 2_IWPF
   !> Multiple grid shifts to eliminate aliasing.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: DEALIAS_GRID_SHIFT = 3_IWPF

   ! Local variables.
   !
   !> Dealiasing technique set by the calling code.
   INTEGER(KIND=IWPF),PRIVATE :: dealias

   ! Variables for working with the FFTW library.
   !
   !> FFT plan for computing the forward (R to C) DFT.
   TYPE(C_PTR),PRIVATE :: r2cPlan
   !> FFT plan for computing the reverse (C to R) DFT.
   TYPE(C_PTR),PRIVATE :: c2rPlan

   ! Module procedures.
   PUBLIC :: GetGridSize, SpectralSetup, SpectralFinalize

CONTAINS

   !> Adjust the working grid size based on the dealiasing technique.
   !!
   !> @param[in] nx Number of grid points in the x direction.
   !> @param[in] ny Number of grid points in the y direction.
   !> @param[in] dealias_ Desired dealiasing technique.
   !> @param[out] nxW Working number of grid points in the x direction.
   !> @param[out] nyW Working number of grid points in the y direction.
   SUBROUTINE GetGridSize(nx, ny, dealias_, nxW, nyW)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nx, ny, dealias_
      INTEGER(KIND=IWPF),INTENT(OUT) :: nxW, nyW

      ! Store the desired dealiasing technique.
      dealias = dealias_
      !
      ! Adjust the working grid size as necessary.
      SELECT CASE (dealias)
         CASE (NO_DEALIASING)
            ! Don't do anything to the grid size.
            nxW = nx
            nyW = ny
         CASE (DEALIAS_3_2_RULE)
         CASE (DEALIAS_2_3_RULE)
         CASE (DEALIAS_GRID_SHIFT)
         CASE DEFAULT
            ! Default will be no dealiasing, but with a warning.
            nxW = nx
            nyW = ny
      END SELECT
   END SUBROUTINE GetGridSize

   !> Initialize the spectral module.
   SUBROUTINE SpectralSetup()
      IMPLICIT NONE
   END SUBROUTINE SpectralSetup

   !> Finalize the spectral module.
   SUBROUTINE SpectralFinalize()
      IMPLICIT NONE

      ! Clean up the FFTW plans.
      CALL FFTW_DESTROY_PLAN(r2cPlan)
      CALL FFTW_DESTROY_PLAN(c2rPlan)
      CALL FFTW_MPI_CLEANUP()
      CALL FFTW_CLEANUP()
   END SUBROUTINE SpectralFinalize

END MODULE Spectral_m

