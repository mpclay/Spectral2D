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
   USE ISO_FORTRAN_ENV,ONLY: OUTPUT_UNIT
   USE MPI
   USE Parameters_m,ONLY: IWPF, RWPF, IWPC, RWPC, CWPC, PI
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
   INTEGER(KIND=IWPF) :: nVar = 2_IWPF
   !> Number of storage locations for time integration.
   INTEGER(KIND=IWPF) :: nStr = 2_IWPF
   !> Number of j columns given to this process by FFTW.
   INTEGER(KIND=IWPF) :: jSize
   !> Local j offset based on FFTW decomposition.
   INTEGER(KIND=IWPF) :: localJ
   !> Starting j index for this process.
   INTEGER(KIND=IWPF) :: j1
   !> Ending j index for this process.
   INTEGER(KIND=IWPF) :: j2
   !> Size of the i dimension for real data arrays.
   INTEGER(KIND=IWPF) :: rISize
   !> Size of the i dimension for complex data arrays.
   INTEGER(KIND=IWPF) :: cISize
   !> Data array for the simulation.
   TYPE(C_PTR) :: Q
   !> Real cast of Q for working in physical space.
   REAL(KIND=RWPC),DIMENSION(:,:,:,:),POINTER :: Qr
   !> Complex cast of Q for working in spectral space.
   COMPLEX(KIND=CWPC),DIMENSION(:,:,:,:),POINTER :: Qc

   ! Extraneous variables.
   !
   !> Looping indices.
   INTEGER(KIND=IWPF) :: i, j, m, n, p
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
   CALL Alloc(nxW, nyW, nVar, nStr, j1, j2, jSize, localJ, rISize, cISize, Q)
   !
   ! Create real and complex casts of Q for use in physical/spectral space.
   CALL C_F_POINTER(Q, Qr, [rISize,nyW,nVar,nStr])
   CALL C_F_POINTER(Q, Qc, [cISize,nyW,nVar,nStr])
   Qc(:,:,:,:) = (0.0_RWPC, 0.0_RWPC)

   ! Initialize the FFT module, which creates the FFTW plans.
   CALL SpectralSetup(nxW, nyW, rISize, cISize, nVar, nStr, Qr, Qc)

   ! Finalize the program.
   CALL Dealloc(Q)
   CALL SpectralFinalize()
   CALL MPI_FINALIZE(ierr)

CONTAINS

   !> Procedure to allocate working arrays.
   !!
   !> @param[in] nxW Working number of grid points in the x direction.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[in] nVar Number of variables in the simulation.
   !> @param[in] nStr Required number of storage arrays for time integration.
   !> @param[out] j1 Starting j index for this process.
   !> @param[out] j2 Ending j index for this process.
   !> @param[out] jSize Number of j columns owned by this MPI process.
   !> @param[out] localJ Offset for j based on FFTW MPI decomposition.
   !> @param[out] rISize Size of i-dimension for real data array.
   !> @param[out] cISize Size of i-dimension for complex data array.
   !> @param[in,out] Q Data array for the simulation.
   SUBROUTINE Alloc(nxW, nyW, nVar, nStr, j1, j2, jSize, localJ, &
                    rISize, cISize, Q)
      ! Required modules.
      USE ISO_C_BINDING
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxW, nyW, nVar, nStr
      INTEGER(KIND=IWPF),INTENT(OUT) :: j1, j2, jSize, localJ, rISize, cISize
      TYPE(C_PTR),INTENT(INOUT) :: Q
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

      ! Main memory allocation.
      !
      ! Each array is allocated by FFTW for complex data. Each array has
      ! dimensions:
      !
      !     (cISize=size of the i dimension in spectral space) X
      !     (nyW=size of the j dimension in physical/spectral space) X
      !     (nVar=number of variables in the simulation) X
      !     (nStr=number of storage sites for time integration).
      Q = FFTW_ALLOC_COMPLEX(INT(cISize*nyW*nVar*nStr, C_SIZE_T))

      ! Subroutine outputs back to the calling code in the proper type.
      jSize = INT(jSize_, IWPF)
      localJ = INT(localJ_, IWPF)
   END SUBROUTINE Alloc

   !> Routine to free memory for the simulation.
   !!
   !> @param[in,out] Q Working array for the simulation.
   SUBROUTINE Dealloc(Q)
      ! Required modules.
      USE ISO_C_BINDING
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'
      ! Calling arguments.
      TYPE(C_PTR),INTENT(INOUT) :: Q

      ! Free memory using FFTW routines.
      CALL FFTW_FREE(Q)
   END SUBROUTINE Dealloc

END PROGRAM Spectral2D_p

