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
   USE Alloc_m,ONLY: GetAllocSize, Alloc, Dealloc
   USE Spectral_m,ONLY: GetGridSize, SpectralSetup, SpectralFinalize, &
                        NO_DEALIASING, DEALIAS_3_2_RULE, DEALIAS_2_3_RULE, &
                        DEALIAS_GRID_SHIFT
   USE TimeIntegration_m,ONLY: TimeIntegrationSetup, IntegrateOneStep, &
                               ComputeTimeStep, RK3_TVD_SHU
   USE SetIC_m,ONLY: SetICCosineShearX, COSINE_SHEAR_X
   USE IO_m,ONLY: FileName, WriteGridPlot3d, WriteRestartPlot3d, &
                  FILE_NAME_LENGTH

   IMPLICIT NONE

   ! FFTW MPI procedure declarations.
   INCLUDE 'fftw3-mpi.f03'

   ! User inputs to the simulation.
   !
   !> Number of grid points in the x direction.
   INTEGER(KIND=IWPF),PARAMETER :: nx = 128_IWPF
   !> Number of grid points in the y direction.
   INTEGER(KIND=IWPF),PARAMETER :: ny = 128_IWPF
   !> Dealiasing technique.
   INTEGER(KIND=IWPF),PARAMETER :: dealias = NO_DEALIASING
   !> Physical viscosity.
   REAL(KIND=RWPC),PARAMETER :: nu = 1.0e-3_RWPC
   !> Time integration scheme.
   INTEGER(KIND=IWPF),PARAMETER :: timeScheme = RK3_TVD_SHU
   !> End time for the simulation.
   REAL(KIND=RWPC),PARAMETER :: tEnd = 0.01_RWPC
   !> Time period after which to write a data file.
   REAL(KIND=RWPC),PARAMETER :: writePeriod = 0.001_RWPC
   !> Time period after which to print information to the user.
   REAL(KIND=RWPC),PARAMETER :: printPeriod = 0.001_RWPC
   !> Desired initial conditions.
   INTEGER(KIND=IWPF),PARAMETER :: ics = COSINE_SHEAR_X

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
   INTEGER(KIND=IWPF) :: nStg
   !> Index for current stage in time in the Q array.
   INTEGER(KIND=IWPF) :: cS
   !> Index for the next stage in time in the Q array.
   INTEGER(KIND=IWPF) :: nS
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
   !> Lowest x wavenumber in the total allocated memory (may be truncated).
   INTEGER(KIND=IWPF) :: kxLT
   !> Lowest y wavenumber in the total allocated memory (may be truncated).
   INTEGER(KIND=IWPF) :: kyLT
   !> Maximum x wavenumber in the total allocated memory (may be truncated).
   INTEGER(KIND=IWPF) :: kxMT
   !> Maximum y wavenumber in the total allocated memory (may be truncated).
   INTEGER(KIND=IWPF) :: kyMT
   !> Lowest x wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kxLU
   !> Lowest y wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kyLU
   !> Maximum x wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kxMU
   !> Maximum y wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kyMU
   !> Data array for the simulation.
   TYPE(C_PTR) :: Q
   !> Real cast of Q for working in physical space.
   REAL(KIND=RWPC),DIMENSION(:,:,:,:),POINTER :: Qr
   !> Complex cast of Q for working in spectral space.
   COMPLEX(KIND=CWPC),DIMENSION(:,:,:,:),POINTER :: Qc
   !> Data array for the u velocity and nonlinear term.
   TYPE(C_PTR) :: u
   !> Real cast of u for working in physical space.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: uR
   !> Complex cast of u for working in spectral space.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: uC
   !> Data array for the v velocity and nonlinear term.
   TYPE(C_PTR) :: v
   !> Real cast of v for working in physical space.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: vR
   !> Complex cast of v for working in spectral space.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: vC

   ! Time stepping variables.
   !
   !> Current simulation time.
   REAL(KIND=RWPF) :: time = 0.0_RWPC
   !> Simulation time step.
   REAL(KIND=RWPC) :: dt = 1.0e-6_RWPC
   !> Current simluation step.
   INTEGER(KIND=IWPF) :: nadv = 0_IWPF

   ! Extraneous variables.
   !
   !> Looping indices.
   INTEGER(KIND=IWPF) :: i, j, m, n, p
   !> MPI error handling.
   INTEGER(KIND=IWPF) :: ierr
   !> Logic to check when we exit time stepping.
   LOGICAL :: loopBool, exitThisStep
   !> Logic to determine when to write data and print info to the user.
   LOGICAL :: printBool, writeBool
   !> Time after which a data file will be written.
   REAL(KIND=RWPC) :: writeTime = writePeriod
   !> Time after which information will be printed to the user.
   REAL(KIND=RWPC) :: printTime = printPeriod
   !> Buffer for output file names.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname

   ! Initialize MPI.
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprc, ierr)

   ! Initialize the FFTW MPI interface.
   CALL FFTW_MPI_INIT()

   ! Determine the working number of grid points based on the dealiasing
   ! technique being used in the simulation.
   CALL GetGridSize(nx, ny, dealias, nxW, nyW, &
                    kxLT, kyLT, kxMT, kyMT, &
                    kxLU, kyLU, kxMU, kyMU)

   ! Get the allocation sizes as determined by FFTW.
   CALL GetAllocSize(nxW, nyW, j1, j2, jSize, localJ, rISize, cISize)

   ! Initialize the time integration module. This also sets nStg for Q.
   CALL TimeIntegrationSetup(timeScheme, rISize, cISize, nyW, &
                             kxLU, kyLU, kxMU, kyMU, nStg, cS, nS)

   ! Allocate memory for the simulation.
   CALL Alloc(cISize, nyW, nVar, nStg, Q)
   CALL Alloc(cISize, nyW, u)
   CALL Alloc(cISize, nyW, v)
   !
   ! Create real and complex casts of Q for use in physical/spectral space.
   CALL C_F_POINTER(Q, Qr, [rISize,nyW,nVar,nStg])
   CALL C_F_POINTER(Q, Qc, [cISize,nyW,nVar,nStg])
   Qc(:,:,:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(u, uR, [rISize,nyW])
   CALL C_F_POINTER(u, uC, [cISize,nyW])
   uC(:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(v, vR, [rISize,nyW])
   CALL C_F_POINTER(v, vC, [cISize,nyW])
   vC(:,:) = (0.0_RWPC, 0.0_RWPC)

   ! Initialize the FFT module, which creates the FFTW plans.
   CALL SpectralSetup(nxW, nyW, rISize, cISize, nVar, nStg, Qr, Qc)

   ! Set the initial conditions.
   SELECT CASE (ics)
      CASE (COSINE_SHEAR_X)
         CALL SetICCosineShearX(rISize, cISize, nxW, nyW, &
                                kxLT, kyLT, kxMT, kyMT, &
                                kxLU, kyLU, kxMU, kyMU, &
                                nVar, nStg, &
                                uC, uR, vC, vR, Qc, Qr)
      CASE DEFAULT
   END SELECT
   !
   ! Write out the grid and the initial conditions.
   CALL WriteGridPlot3d(nxW, nyW)
   CALL FileName('REST', 0, 'f', fname)
   CALL WriteRestartPlot3d(rISize, cISize, nxW, nyW, &
                           kxLT, kyLT, kxMT, kyMT, &
                           kxLU, kyLU, kxMU, kyMU, &
                           nVar, nStg, cS, &
                           uC, uR, vC, vR, Qc, Qr, fname)
   STOP

   ! Enter the main time stepping loop.
   loopBool = .TRUE.
   exitThisStep = .FALSE.
   tloop: DO WHILE (loopBool)
      ! Increment the Q array by 1 time step.
      !
      ! NOTE: The only arrays that are correct after this step are Qr and Qc.
      CALL IntegrateOneStep(rISize, cISize, nxW, nyW, &
                            kxLT, kyLT, kxMT, kyMT, &
                            kxLU, kyLU, kxMU, kyMU, &
                            nu, dt, &
                            nVar, nStg, cS, nS, &
                            uC, uR, vC, vR, Qc, Qr)

      ! Increment the time and step counters.
      time = time + dt
      nadv = nadv + 1_IWPF

      ! Compute the next time step.
      CALL ComputeTimeStep()

      ! Check whether or not we will write out data to file.
      IF (writeBool) THEN
         writeTime = writeTime + writePeriod
         writeBool = .FALSE.
      ELSE
         IF (time + dt >= writeTime) THEN
            writeBool = .TRUE.
         END IF
      END IF

      ! Check whether or not we will print information to the user.
      IF (printBool .OR. exitThisStep) THEN
         printTime = printTime + printPeriod
         printBool = .FALSE.
         WRITE(*,500) 'Simulation step number: ', nadv, &
                      '; Simulation time: ', time, &
                      '; Max W: ', MAXVAL(ABS(Qc(:,:,1,cS))), &
                      '; Min W: ', MINVAL(ABS(Qc(:,:,1,cS))), &
                      '; Max Psi: ', MAXVAL(ABS(Qc(:,:,2,cS))), &
                      '; Min Psi: ', MINVAL(ABS(Qc(:,:,2,cS)))
         500 FORMAT (A,I8.8,A,ES15.8,A,ES15.8,A,ES15.8,A,ES15.8,A,ES15.8)
      ELSE
         IF (time + dt >= printTime) THEN
            printBool = .TRUE.
         END IF
      END IF

      ! Evaluate the exit conditions.
      IF (exitThisStep) THEN
         EXIT tloop
      END IF
      !
      ! Adjust the time step so the desired end time is reached.
      IF (time + dt > tend) THEN
         dt = tend - time
         exitThisStep = .TRUE.
      END IF
   END DO tloop

   ! Finalize the program.
   CALL Dealloc(Q)
   CALL SpectralFinalize()
   CALL MPI_FINALIZE(ierr)

END PROGRAM Spectral2D_p

