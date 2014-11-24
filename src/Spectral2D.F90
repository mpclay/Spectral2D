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
!> @file Spectral2D.F90
!> @author Matthew Clay
!> @brief 2D Navier-Stokes solver using a pseudo-spectral method.
PROGRAM Spectral2D_p

   ! Required modules.
   USE ISO_C_BINDING
   USE MPI
   USE Parameters_m,ONLY: IWPF, RWPF, IWPC, RWPC, PI
   USE Spectral_m,ONLY: GetGridSize, SpectralSetup, SpectralFinalize, &
                        NO_DEALIASING, DEALIAS_3_2_RULE, DEALIAS_2_3_RULE, &
                        DEALIAS_GRID_SHIFT

   IMPLICIT NONE

   ! FFTW MPI procedure declarations.
   INCLUDE 'fftw3-mpi.f03'

   ! User inputs to the simulation.
   !
   !> Number of grid points in the x direction.
   INTEGER(KIND=IWPF),PARAMETER :: nx = 256_IWPF
   !> Number of grid points in the y direction.
   INTEGER(KIND=IWPF),PARAMETER :: ny = 256_IWPF
   !> Dealiasing technique.
   INTEGER(KIND=IWPF),PARAMETER :: dealias = NO_DEALIASING

   ! MPI related variables.
   !
   !> Process ID.
   INTEGER(KIND=IWPF) :: rank
   !> Number of processes in MPI_COMM_WORLD.
   INTEGER(KIND=IWPF) :: nprc

   ! Variables required for simulation.
   !
   !> Working number of grid points in the x direction.
   INTEGER(KIND=IWPF) :: nxW
   !> Working number of grid points in the y direction.
   INTEGER(KIND=IWPF) :: nyW
   !> Number of variables in the Q array for the simulation.
   INTEGER(KIND=IWPF) :: nVar
   !> Number of storage locations for time integration.
   INTEGER(KIND=IWPF) :: nStr
   !> Number of j columns given to this process by FFTW.
   INTEGER(KIND=IWPF) :: jSize
   !> Local j offset based on FFTW decomposition.
   INTEGER(KIND=IWPF) :: localJ
   !> Data array for the simulation.
   TYPE(C_PTR) :: Q

   ! Extraneous variables.
   !
   !> Looping indices.
   INTEGER(KIND=IWPF) :: i, j, m, n
   !> MPI error handling.
   INTEGER(KIND=IWPF) :: ierr

   ! Initialize MPI.
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprc, ierr)

   ! Initialize the FFTW MPI interface.
   CALL FFTW_MPI_INIT()

   ! Determine the working number of grid points based on the dealiasing
   ! technique being used in the simulation.
   CALL GetGridSize(nx, ny, dealias, nxW, nyW)

   ! Allocate memory for the simulation, and get the local grid extents.
   CALL Alloc(nxW, nyW, nVar, nStr, jSize, localJ, Q)

   ! Finalize the program.
   CALL SpectralFinalize()
   CALL MPI_FINALIZE(ierr)

CONTAINS

   !> Procedure to allocate working arrays.
   !!
   !> @param[in] nxW Working number of grid points in the x direction.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[in] nVar Number of variables in the simulation.
   !> @param[in] nStr Required number of storage arrays for time integration.
   !> @param[out] jSize Number of j columns owned by this MPI process.
   !> @param[out] localJ Offset for j based on FFTW MPI decomposition.
   !> @param[in,out] Q Data array for the simulation.
   SUBROUTINE Alloc(nxW, nyW, nVar, nStr, jSize, localJ, Q)
      ! Required modules.
      USE ISO_C_BINDING
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxW, nyW, nVar, nStr
      INTEGER(KIND=IWPF),INTENT(OUT) :: jSize, localJ
      TYPE(C_PTR),INTENT(INOUT) :: Q
      ! Local variables.
      ! Size of working array for this process as determined by FFTW.
      INTEGER(KIND=IWPC) :: allocLocal
      ! Local starting index for this task's portion of the data.
      INTEGER(KIND=IWPC) :: localJ_
      ! Number of columns (starting from localX) this process owns.
      INTEGER(KIND=IWPC) :: jSize_

      ! Determine number of columns this MPI task gets for the FFTs.
      allocLocal = FFTW_MPI_LOCAL_SIZE_2D(INT(nxW, IWPC), &
                                          INT(nyW, IWPC), &
                                          MPI_COMM_WORLD, &
                                          localJ_, jSize_)
      !
      ! Increase allocLocal according to the number of variables in the
      ! simulation and the number of storage locations for time int.
      allocLocal = allocLocal*INT(nVar, IWPC)*INT(nStr, IWPC)

      ! Subroutine outputs back to the calling code in the proper type.
      jSize = INT(jSize_, IWPF)
      localJ = INT(localJ_, IWPF)
   END SUBROUTINE Alloc

END PROGRAM Spectral2D_p

