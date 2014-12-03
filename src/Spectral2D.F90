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
!!
!! This code solves the two-dimensional incompressible Navier-Stokes equations
!! in streamfunction-vorticity formulation with a pseudo-spectral scheme. The
!! equations are solved in Fourier space, with nonlinear terms being evaluated
!! pseudo-spectrally. Time integration is performed with Shu's TVD RK3 scheme.
!!
!! The DFTs required for the psuedo-spectral scheme are conducted using the
!! MPI-parallelized version of FFTW. All transforms are in-place. Currently,
!! a final data transposition is used in FFTW to realign the i memory direction
!! with the x coordinate direction. This could be relaxed to increase speed.
PROGRAM Spectral2D_p

   ! Required modules.
   USE ISO_C_BINDING
   USE ISO_FORTRAN_ENV,ONLY: OUTPUT_UNIT
   USE MPI
   USE Parameters_m,ONLY: IWPF, RWPF, IWPC, RWPC, CWPC, PI
   USE Alloc_m,ONLY: GetAllocSize, Alloc, Dealloc
   USE Analysis_m,ONLY: ComputeSpectrum
   USE Spectral_m,ONLY: GetGridSize, SetupMask, SetupTruncation, SetupPlans, &
                        SpectralFinalize, TransformR2C, TransformC2R, &
                        ComputeVelocity, NO_DEALIASING, DEALIAS_3_2_RULE, &
                        DEALIAS_2_3_RULE
   USE TimeIntegration_m,ONLY: TimeIntegrationSetup, IntegrateOneStep, &
                               ComputeTimeStep, RK3_TVD_SHU, &
                               TimeIntegrationFinalize
   USE SetIC_m,ONLY: SetTaylorGreen, SetTurbulentVelocity, &
                     TAYLOR_GREEN_VORTEX, SCHUMANN_VELOCITY_FIELD
   USE IO_m,ONLY: FileName, WriteGridHDF5, WriteRestart, FILE_NAME_LENGTH, &
                  HDF5_OUTPUT

   IMPLICIT NONE

   ! FFTW MPI procedure declarations.
   INCLUDE 'fftw3-mpi.f03'

   ! Methods by which to determine the time step.
   !
   !> Constant time stepping.
   INTEGER(KIND=IWPF),PARAMETER :: CONSTANT_TIME_STEP = 1_IWPF
   !> Time step determined by the CFL criterion.
   INTEGER(KIND=IWPF),PARAMETER :: DYNAMIC_TIME_STEP = 2_IWPF
   !
   ! Methods to determine when the simulation ends.
   !
   !> Fixed number of time step.
   INTEGER(KIND=IWPF),PARAMETER :: FIXED_STEP_COUNT = 1_IWPF
   !> Fixed end time for the simulation.
   INTEGER(KIND=IWPF),PARAMETER :: DESIRED_END_TIME = 2_IWPF

   ! User inputs to the simulation.
   !
   !> Number of grid points in the x direction.
   INTEGER(KIND=IWPF),PARAMETER :: nx = 64_IWPF
   !> Number of grid points in the y direction.
   INTEGER(KIND=IWPF),PARAMETER :: ny = 64_IWPF
   !
   ! Physical parameters and spatial scheme control.
   !
   !> Physical viscosity.
   REAL(KIND=RWPC),PARAMETER :: nu = 0.001_RWPC
   !> Dealiasing technique.
   INTEGER(KIND=IWPF),PARAMETER :: dealias = DEALIAS_3_2_RULE
   !
   ! Time integration and simulation ending method.
   !
   !> Time integration scheme.
   INTEGER(KIND=IWPF),PARAMETER :: timeScheme = RK3_TVD_SHU
   !> Desired simulation ending criterion.
   INTEGER(KIND=IWPF),PARAMETER :: endMethod = FIXED_STEP_COUNT
   !
   ! Parameters required for constant time stepping.
   !
   !> Number of steps for the simulation.
   INTEGER(KIND=IWPF),PARAMETER :: stepEnd = 10000000_IWPF
   !> Simulation time step.
   REAL(KIND=RWPC) :: dt = 1.0e-6_RWPC
   !
   ! Parameters required for dynamic time stepping.
   !
   !> End time for the simulation.
   REAL(KIND=RWPC) :: tEnd = 1.0_RWPC
   !
   ! Parameters to determine when data will be written.
   !
   !> Number of steps after which to write a data file.
   INTEGER(KIND=IWPF),PARAMETER :: writePeriod = 10000_IWPF
   !> Number of steps after which to print info to the user.
   INTEGER(KIND=IWPF),PARAMETER :: printPeriod = 1000_IWPF
   !
   ! Parameters to control initial conditions.
   !
   !> Desired initial conditions.
   INTEGER(KIND=IWPF),PARAMETER :: ics = SCHUMANN_VELOCITY_FIELD
   !
   ! Additional parameters required for turbulent initialization.
   !
   !> Desired peak wavenumber in the spectrum.
   REAL(KIND=RWPF),PARAMETER :: k0 = 4.0_RWPF
   !> RMS velocity intensity.
   REAL(KIND=RWPF),PARAMETER :: urms = 1.0_RWPF

   ! MPI related variables.
   !
   !> Process ID.
   INTEGER(KIND=IWPF) :: rank
   !> Number of processes in MPI_COMM_WORLD.
   INTEGER(KIND=IWPF) :: nprc
   !> Send buffers for determining max/min of vorticity and psi.
   REAL(KIND=RWPF),DIMENSION(2) :: maxSendBuff, minSendBuff
   !> Receive buffers for determining max/min of vorticity and psi.
   REAL(KIND=RWPF),DIMENSION(2) :: maxRecvBuff, minRecvBuff
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
   !> Maximum wavenumber magnitude allowed before truncation.
   REAL(KIND=RWPF) :: kMax
   !> Array of wavenumbers in the x direction.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kxG
   !> Array of wavenumbers in the y direction.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kyG
   !
   ! Local grid and wavenumber variables.
   !
   !> Number of grid points in the x direction for this process.
   INTEGER(KIND=IWPF) :: nxP
   !> Number of j columns given to this process by FFTW.
   INTEGER(KIND=IWPF) :: nyP
   !> Local j offset based on FFTW decomposition.
   INTEGER(KIND=IWPF) :: localJ
   !> Starting j index for this process.
   INTEGER(KIND=IWPF) :: j1
   !> Ending j index for this process.
   INTEGER(KIND=IWPF) :: j2
   !> Size of data arrays as required by FFTW (includes MPI padding).
   INTEGER(KIND=IWPF) :: allocLocal
   !> Size of the i dimension for real data arrays.
   INTEGER(KIND=IWPF) :: rISize
   !> Size of the i dimension for complex data arrays.
   INTEGER(KIND=IWPF) :: cISize
   !> Lowest x wavenumber on this process.
   INTEGER(KIND=IWPF) :: kxLP
   !> Maximum x wavenumber on this process.
   INTEGER(KIND=IWPF) :: kxMP
   !> Lowest y wavenumber on this process.
   INTEGER(KIND=IWPF) :: kyLP
   !> Maximum y wavenumber on this process.
   INTEGER(KIND=IWPF) :: kyMP
   !> Array of wavenumbers in the x direction for this process.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kxP
   !> Array of wavenumbers in the y direction for this process.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kyP
   !
   ! Main working arrays for the simulation.
   !
   !> Data array for the vorticity.
   TYPE(C_PTR) :: w
   !> Real cast of the vorticity.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: wR
   !> Complex cast of the vorticity.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: wC
   !> Data array for the streamfunction.
   TYPE(C_PTR) :: psi
   !> Real cast of the streamfunction.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: psiR
   !> Complex cast of the streamfunction.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: psiC
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
   !> Current simluation step.
   INTEGER(KIND=IWPF) :: nadv = 0_IWPF
   !
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
   !> Number of steps after which a data file will be written.
   INTEGER(KIND=IWPF) :: writeStep = writePeriod
   !> Number of steps after which information will be printed to the user.
   INTEGER(KIND=IWPF) :: printStep = printPeriod
   !> Buffer for output file names.
   CHARACTER(LEN=FILE_NAME_LENGTH) :: fname

   ! Initialize MPI.
   CALL MPI_INIT(ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprc, ierr)

   ! Initialize the FFTW MPI interface.
   CALL FFTW_MPI_INIT()

   ! Determine the working number of grid points based on the dealiasing
   ! technique being used in the simulation. Also calculated min/max wavenums.
   CALL GetGridSize(nx, ny, dealias, nxG, nyG, &
                    kxLG, kyLG, kxMG, kyMG, &
                    kxLU, kyLU, kxMU, kyMU, &
                    kMax, rISize, cISize, kxG, kyG)

   ! Get the allocation sizes as determined by the MPI FFTW routines.
   CALL GetAllocSize(nxG, nyG, kxG, kyG, kxMU, cISize, nxP, nyP, j1, j2, &
                     localJ, allocLocal, kxLP, kyLP, kxMP, kyMP, kxP, kyP)

   ! Set up masked wavenumbers for differentiation.
   CALL SetupMask(cISize, nyP, kxP, kyP, kxMG, kyMG)

   ! Set up wavenumber truncation bounds.
   CALL SetupTruncation(cISize, nyP, kMax, kxP, kyP)

   ! Write out the grid.
   CALL FileName('GRID', 'h5', fname)
   CALL WriteGridHDF5(MPI_COMM_WORLD, rank, fname, nxG, nyG, nxP, nyP, &
                      1_IWPF, nxP, j1, j2)

   ! Initialize the time integration module.
   CALL TimeIntegrationSetup(timeScheme, rISize, cISize, nyP, allocLocal)

   ! Allocate memory for the simulation.
   CALL Alloc(allocLocal, w)
   CALL Alloc(allocLocal, psi)
   CALL Alloc(allocLocal, u)
   CALL Alloc(allocLocal, v)
   !
   ! Create real/complex casts of arrays for use in physical/spectral space.
   CALL C_F_POINTER(w, wR, [rISize,nyP])
   CALL C_F_POINTER(w, wC, [cISize,nyP])
   wC(:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(psi, psiR, [rISize,nyP])
   CALL C_F_POINTER(psi, psiC, [cISize,nyP])
   psiC(:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(u, uR, [rISize,nyP])
   CALL C_F_POINTER(u, uC, [cISize,nyP])
   uC(:,:) = (0.0_RWPC, 0.0_RWPC)
   CALL C_F_POINTER(v, vR, [rISize,nyP])
   CALL C_F_POINTER(v, vC, [cISize,nyP])
   vC(:,:) = (0.0_RWPC, 0.0_RWPC)

   ! Initialize the FFT module, which creates the FFTW plans.
   CALL SetupPlans(nxG, nyG, nyP, rISize, cISize, wR, wC)

   ! Set the initial conditions.
   SELECT CASE (ics)
      CASE (TAYLOR_GREEN_VORTEX)
         CALL SetTaylorGreen(nxG, nyG, nxP, nyP, j1, j2, rISize, cISize, &
                             wC, wR, psiC, psiR)
      CASE (SCHUMANN_VELOCITY_FIELD)
         CALL SetTurbulentVelocity(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                                   kxP, kyP, kxLG, kyLG, kxMG, kyMG, k0, urms, &
                                   uC, vC, wC, psiC)
      CASE DEFAULT
   END SELECT
   !
   ! Write out a restart file.
   CALL WriteRestart(nxG, nyG, nxP, nyP, 1_IWPF, nxP, j1, j2, &
                     rISize, cISize, kxP, kyP, rank, nadv, time, nu, &
                     HDF5_OUTPUT, uC, uR, vC, vR, wC, wR, psiC, psiR)
   !
   ! Also write out the initial energy spectrum.
   CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                        uC, uR, vC, vR, PsiC, .FALSE.)
   CALL ComputeSpectrum(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                        kxP, kyP, kxLG, kyLG, kxMG, kyMG, .TRUE., &
                        0, uC, vC)

   ! Enter the main time stepping loop.
   loopBool = .TRUE.
   exitThisStep = .FALSE.
   writeBool = .FALSE.
   printBool = .FALSE.
   tloop: DO WHILE (loopBool)
      ! Increment the vorticity and streamfunction arrays by one time step.
      !
      ! NOTE: The only arrays that are correct after this step are wC and psiC.
      CALL IntegrateOneStep(nxG, nyG, nxP, nyP, rISize, cISize, kxP, kyP, &
                            nu, dt, uC, uR, vC, vR, wC, wR, psiC, psiR)

      ! Increment the time and step counters.
      time = time + dt
      nadv = nadv + 1_IWPF

      ! Compute the next time step.
      CALL ComputeTimeStep()

      ! Check whether or not we will print information to the user.
      IF (printBool .OR. exitThisStep) THEN
         printStep = printStep + printPeriod
         printBool = .FALSE.
         !
         ! Each process calculates the min/max of wC and psiC.
         maxSendBuff(1) = MAXVAL(ABS(wC(1:cISize,1:nyP)))
         maxSendBuff(2) = MAXVAL(ABS(psiC(1:cISize,1:nyP)))
         minSendBuff(1) = MINVAL(ABS(wC(1:cISize,1:nyP)))
         minSendBuff(2) = MINVAL(ABS(psiC(1:cISize,1:nyP)))
         !
         ! Reduce the values to the root process.
         CALL MPI_REDUCE(maxSendBuff, maxRecvBuff, 2, MPI_DOUBLE, MPI_MAX, &
                         0, MPI_COMM_WORLD, ierr)
         CALL MPI_REDUCE(minSendBuff, minRecvBuff, 2, MPI_DOUBLE, MPI_MIN, &
                         0, MPI_COMM_WORLD, ierr)
         !
         ! The root process prints the info to the screen for monitoring.
         IF (rank == 0_IWPF) THEN
            WRITE(OUTPUT_UNIT,500) 'Simulation step number: ', nadv, &
                                   '; Simulation time: ', time, &
                                   '; Max W: ', maxRecvBuff(1), &
                                   '; Min W: ', minRecvBuff(1), &
                                   '; Max Psi: ', maxRecvBuff(2), &
                                   '; Min Psi: ', minRecvBuff(2)
            500 FORMAT (A,I8.8,A,ES15.8,A,ES15.8,A,ES15.8,A,ES15.8,A,ES15.8)
         END IF
      ELSE
         IF (nadv + 1_IWPF == printStep) THEN
            printBool = .TRUE.
         END IF
      END IF

      ! Check whether or not we will write out data to file.
      IF (writeBool .OR. exitThisStep) THEN
         writeStep = writeStep + writePeriod
         writeBool = .FALSE.
         !
         ! Write out the restart file in HDF5 format.
         CALL WriteRestart(nxG, nyG, nxP, nyP, 1_IWPF, nxP, j1, j2, &
                           rISize, cISize, kxP, kyP, rank, nadv, time, nu, &
                           HDF5_OUTPUT, uC, uR, vC, vR, wC, wR, psiC, psiR)
         !
         ! Also write out the energy spectrum when a restart file is saved.
         CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                              uC, uR, vC, vR, PsiC, .FALSE.)
         CALL ComputeSpectrum(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                              kxP, kyP, kxLG, kyLG, kxMG, kyMG, .TRUE., &
                              nadv, uC, vC)
      ELSE
         IF (nadv + 1_IWPF == writeStep) THEN
            writeBool = .TRUE.
         END IF
      END IF

      ! Evaluate the exit conditions.
      IF (exitThisStep) THEN
         EXIT tloop
      END IF
      !
      ! Check when we should exit.
      SELECT CASE (endMethod)
         CASE (FIXED_STEP_COUNT)
            IF (nadv + 1_IWPF == stepEnd) THEN
               exitThisStep = .TRUE.
            END IF
         CASE (DESIRED_END_TIME)
            IF (time + dt >= tEnd) THEN
               exitThisStep = .TRUE.
               dt = tend - time
            END IF
      END SELECT
   END DO tloop

   ! Finalize the program.
   CALL Dealloc(w)
   CALL Dealloc(psi)
   CALL Dealloc(u)
   CALL Dealloc(v)
   wR => NULL()
   wC => NULL()
   psiR => NULL()
   psiC => NULL()
   uR => NULL()
   uC => NULL()
   vR => NULL()
   vC => NULL()
   DEALLOCATE(kxG)
   DEALLOCATE(kyG)
   DEALLOCATE(kxP)
   DEALLOCATE(kyP)
   CALL TimeIntegrationFinalize()
   CALL SpectralFinalize()
   CALL MPI_FINALIZE(ierr)

END PROGRAM Spectral2D_p

