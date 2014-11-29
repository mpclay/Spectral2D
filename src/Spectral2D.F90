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
                        TransformR2C, TransformC2R, ComputeVelocity, &
                        NO_DEALIASING, DEALIAS_3_2_RULE, DEALIAS_2_3_RULE
   USE TimeIntegration_m,ONLY: TimeIntegrationSetup, IntegrateOneStep, &
                               ComputeTimeStep, RK3_TVD_SHU
   USE SetIC_m,ONLY: SetICCosineShearX, COSINE_SHEAR_X
   USE IO_m,ONLY: FileName, WriteGridHDF5, WriteRestart, FILE_NAME_LENGTH, &
                  HDF5_OUTPUT

   IMPLICIT NONE

   ! FFTW MPI procedure declarations.
   INCLUDE 'fftw3-mpi.f03'

   ! User inputs to the simulation.
   !
   !> Number of grid points in the x direction.
   INTEGER(KIND=IWPF),PARAMETER :: nx = 32_IWPF
   !> Number of grid points in the y direction.
   INTEGER(KIND=IWPF),PARAMETER :: ny = 32_IWPF
   !> Dealiasing technique.
   INTEGER(KIND=IWPF),PARAMETER :: dealias = NO_DEALIASING
   !> Physical viscosity.
   REAL(KIND=RWPC),PARAMETER :: nu = 1.0e-2_RWPC
   !> Time integration scheme.
   INTEGER(KIND=IWPF),PARAMETER :: timeScheme = RK3_TVD_SHU
   !> End time for the simulation.
   REAL(KIND=RWPC),PARAMETER :: tEnd = 100.0_RWPC
   !> Time period after which to write a data file.
   REAL(KIND=RWPC),PARAMETER :: writePeriod = 10.0_RWPC
   !> Time period after which to print information to the user.
   REAL(KIND=RWPC),PARAMETER :: printPeriod = 0.1_RWPC
   !> Desired initial conditions.
   INTEGER(KIND=IWPF),PARAMETER :: ics = COSINE_SHEAR_X

   ! MPI related variables.
   !
   !> Process ID.
   INTEGER(KIND=IWPF) :: rank
   !> Number of processes in MPI_COMM_WORLD.
   INTEGER(KIND=IWPF) :: nprc

   ! Variables related to number of variables, storage size, etc.
   !
   !> Number of variables in the Q array for the simulation.
   INTEGER(KIND=IWPF) :: nVar = 2_IWPF
   !> Number of storage locations for time integration.
   INTEGER(KIND=IWPF) :: nStg
   !> Index for current stage in time in the Q array.
   INTEGER(KIND=IWPF) :: cS
   !> Index for the next stage in time in the Q array.
   INTEGER(KIND=IWPF) :: nS
   !
   ! Global grid and wavenumber variables.
   !
   !> Total number of grid points in the x direction.
   INTEGER(KIND=IWPF) :: nxG
   !> Total number of grid points in the y direction.
   INTEGER(KIND=IWPF) :: nyG
   !> Lowest x wavenumber in the total allocated memory.
   INTEGER(KIND=IWPF) :: kxLG
   !> Lowest y wavenumber in the total allocated memory.
   INTEGER(KIND=IWPF) :: kyLG
   !> Maximum x wavenumber in the total allocated memory.
   INTEGER(KIND=IWPF) :: kxMG
   !> Maximum y wavenumber in the total allocated memory.
   INTEGER(KIND=IWPF) :: kyMG
   !> Lowest x wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kxLU
   !> Lowest y wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kyLU
   !> Maximum x wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kxMU
   !> Maximum y wavenumber actually being used.
   INTEGER(KIND=IWPF) :: kyMU
   !> Array of wavenumbers in the x direction.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kxG
   !> Array of wavenumbers in the y direction.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kyG
   !
   ! Local grid and wavenumber variables.
   !
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
   !
   ! Main working arrays for the simulation.
   !
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
   !
   ! Time stepping variables.
   !
   !> Current simulation time.
   REAL(KIND=RWPF) :: time = 0.0_RWPC
   !> Simulation time step.
   REAL(KIND=RWPC) :: dt = 1.0e-5_RWPC
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
   CALL GetGridSize(nx, ny, dealias, nxG, nyG, &
                    kxLG, kyLG, kxMG, kyMG, &
                    kxLU, kyLU, kxMU, kyMU, &
                    rISize, cISize, kxG, kyG)

   ! Write out the grid.
   CALL FileName('GRID', 'h5', fname)
   CALL WriteGridHDF5(MPI_COMM_WORLD, rank, fname, nxG, nyG, nxG, nyG, &
                      1_IWPF, nxG, 1_IWPF, nyG)

   ! Get the allocation sizes as determined by FFTW.
   CALL GetAllocSize(nxG, nyG, j1, j2, jSize, localJ)

   ! Initialize the time integration module. This also sets nStg for Q.
   CALL TimeIntegrationSetup(timeScheme, rISize, cISize, nyG, nStg, cS, nS)

   ! Allocate memory for the simulation.
   CALL Alloc(cISize, nyG, nVar, nStg, Q)
   CALL Alloc(cISize, nyG, u)
   CALL Alloc(cISize, nyG, v)
   !
   ! Create real and complex casts of Q for use in physical/spectral space.
   CALL C_F_POINTER(Q, Qr, [rISize,nyG,nVar,nStg])
   CALL C_F_POINTER(Q, Qc, [cISize,nyG,nVar,nStg])
   Qc(:,:,:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(u, uR, [rISize,nyG])
   CALL C_F_POINTER(u, uC, [cISize,nyG])
   uC(:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(v, vR, [rISize,nyG])
   CALL C_F_POINTER(v, vC, [cISize,nyG])
   vC(:,:) = (0.0_RWPC, 0.0_RWPC)

   ! Initialize the FFT module, which creates the FFTW plans.
   CALL SpectralSetup(nxG, nyG, rISize, cISize, nVar, nStg, Qr, Qc)

   ! Set the initial conditions.
   SELECT CASE (ics)
      CASE (COSINE_SHEAR_X)
         CALL SetICCosineShearX(nxG, nyG, rISize, cISize, kxG, kyG, &
                                uC, uR, vC, vR, Qc(:,:,1,1), Qr(:,:,1,1), &
                                Qc(:,:,2,1), Qr(:,:,2,1))
      CASE DEFAULT
   END SELECT
   !
   ! Write out a restart file.
   CALL WriteRestart(nxG, nyG, nxG, nyG, 1_IWPF, 1_IWPF, 1_IWPF, 1_IWPF, &
                     rISize, cISize, kxG, kyG, rank, nadv, time, HDF5_OUTPUT, &
                     uC, uR, vC, vR, Qc(:,:,1,1), Qr(:,:,1,1), &
                     Qc(:,:,2,1), Qr(:,:,2,1))

   ! Enter the main time stepping loop.
   loopBool = .TRUE.
   exitThisStep = .FALSE.
   writeBool = .FALSE.
   printBool = .FALSE.
   tloop: DO WHILE (loopBool)
      ! Increment the Q array by 1 time step.
      !
      ! NOTE: The only arrays that are correct after this step are Qr and Qc.
      CALL IntegrateOneStep(nxG, nyG, rISize, cISize, kxG, kyG, &
                            nu, dt, nVar, nStg, cS, nS, &
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
         CALL WriteRestart(nxG, nyG, nxG, nyG, 1_IWPF, 1_IWPF, 1_IWPF, 1_IWPF, &
                           rISize, cISize, kxG, kyG, rank, nadv, time, &
                           HDF5_OUTPUT, uC, uR, vC, vR, Qc(:,:,1,1), &
                           Qr(:,:,1,1), Qc(:,:,2,1), Qr(:,:,2,1))
      ELSE
         IF (time + dt >= writeTime) THEN
            writeBool = .TRUE.
         END IF
      END IF

      ! Check whether or not we will print information to the user.
      IF (printBool .OR. exitThisStep) THEN
         printTime = printTime + printPeriod
         printBool = .FALSE.
         WRITE(OUTPUT_UNIT,500) 'Simulation step number: ', nadv, &
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

