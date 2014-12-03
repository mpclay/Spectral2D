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
!> @file SetIC.F90
!> @author Matthew Clay
!> @brief Module to set the initial conditions for the simulation.
!!
!! Currently we provide initialization routines for:
!!
!!    1. Taylor-Green vortex.
!!    2. Velocity field fitting a desired energy spectrum.
MODULE SetIC_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPF, IWPC, RWPF, RWPC, CWPC

   IMPLICIT NONE

   ! Available initial conditions.
   !
   !> Taylor-Green vortex initialization.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: TAYLOR_GREEN_VORTEX = 1_IWPF
   !> Turbulent velocity field initialization with Schumann spectrum.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: SCHUMANN_VELOCITY_FIELD = 2_IWPF

   ! Module procedures.
   PUBLIC :: SetTaylorGreen, SetTurbulentVelocity
   PRIVATE :: FourierVelocity, Schumann

CONTAINS

   !> Procedure to set the initial condition for the Taylor-Green vortex.
   !!
   !! The Taylor-Green vortex is an exact solution to the Navier-Stokes
   !! equations in a periodic domain. The flow functions of interest satisfy:
   !!
   !!    1. Psi(x,y,t) = sin(x)sin(y)F(t)
   !!    2. w(x,y,t) = 2sin(x)sin(y)F(t)
   !!    3. u(x,y,t) = sin(x)cos(y)F(t)
   !!    4. v(x,y,t) = -cos(x)sin(y)F(t),
   !!
   !! where F(t)=exp(-2*nu*t), and nu is the kinematic viscosity. Since this is
   !! an initial condition, F(0)=1 is used.
   !!
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction for this process.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] j1 Starting j index for this process.
   !> @param[in] j2 Ending j index for this process.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] wC Complex cast of vorticity array.
   !> @param[in] wR Real cast of vorticity array.
   !> @param[in] psiC Complex cast of streamfunction array.
   !> @param[in] psiR Real cast of streamfunction array.
   SUBROUTINE SetTaylorGreen(nxG, nyG, nxP, nyP, j1, j2, rISize, cISize, &
                             wC, wR, psiC, psiR)
      ! Required modules.
      USE Parameters_m,ONLY: PI
      USE Spectral_m,ONLY: TransformR2C
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nxP, nyP, j1, j2
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: wC, psiC
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: wR, psiR
      ! Local variables.
      ! Looping indices for physical space.
      INTEGER(KIND=IWPF) :: i, j
      ! Index to access the j rank of working arrays.
      INTEGER(KIND=IWPF) :: jI

      ! Loop over the physical grid points and fill in psi and vorticity.
      DO j = j1, j2
         jI = j - j1 + 1_IWPF
         DO i = 1, nxP
            psiR(i,jI) = SIN(2.0_RWPC*PI*REAL(i, RWPC)/REAL(nxG, RWPC))* &
                         SIN(2.0_RWPC*PI*REAL(j, RWPC)/REAL(nyG, RWPC))
            wR(i,jI) = 2.0_RWPC*SIN(2.0_RWPC*PI*REAL(i, RWPC)/REAL(nxG, RWPC))*&
                                SIN(2.0_RWPC*PI*REAL(j, RWPC)/REAL(nyG, RWPC))
         END DO
      END DO
      !
      ! Transform the signal to spectral space.
      CALL TransformR2C(nxG, nyG, nyP, rISize, cISize, psiR, psiC)
      CALL TransformR2C(nxG, nyG, nyP, rISize, cISize, wR, wC)
   END SUBROUTINE SetTaylorGreen

   !> Set the initial conditions of the simulation to have a velocity field
   !! fitting a desired energy spectrum.
   !!
   !> @param[in] rank MPI rank of this process.
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction on this process.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] rISize Size of real data in the i direction.
   !> @param[in] cISize Size of complex data in the i direction.
   !> @param[in] kxP Array of wavenumbers in the x direction on this process.
   !> @param[in] kyP Array of wavenumbers in the y direction on this process.
   !> @param[in] kxLG Lowest x wavenumber in the total simulation.
   !> @param[in] kyLG Lowest y wavenumber in the total simulation.
   !> @param[in] kxMG Maximum x wavenumber in the total simulation.
   !> @param[in] kyMG Maximum y wavenumber in the total simulation.
   !> @param[in] k0 Desired peak wavenumber in the spectrum.
   !> @param[in] urms Desired value of urms for the flow.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[in,out] wC Complex cast of vorticity array.
   !> @param[in,out] psiC Complex cast of streamfunction array.
   SUBROUTINE SetTurbulentVelocity(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                                   kxP, kyP, kxLG, kyLG, kxMG, kyMG, k0, urms, &
                                   uC, vC, wC, psiC)
      ! Required modules.
      USE Spectral_m,ONLY: ComputeVorticity, ComputePsi
      USE Analysis_m,ONLY: ComputeSpectrum
      USE Random_m,ONLY: InitRandomSeed
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rank, nxG, nyG, nxP, nyP, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLG, kyLG, kxMG, kyMG
      REAL(KIND=RWPF),INTENT(IN) :: k0, urms
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC, wC, psiC

      ! Initialize the random seed.
      CALL InitRandomSeed(rank)

      ! Initialize the Fourier velocity field.
      CALL FourierVelocity(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                           kxP, kyP, kxLG, kyLG, kxMG, kyMG, k0, urms, uC, vC)

      ! Calculate the initial energy spectrum.
      CALL ComputeSpectrum(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                           kxP, kyP, kxLG, kyLG, kxMG, kyMG, .TRUE., 0, uC, vC)

      ! Calculate the vorticity from the velocity field.
      CALL ComputeVorticity(nxP, nyP, cISize, kxP, kyP, uC, vC, wC)

      ! Calculate the streamfunction from the vorticity.
      CALL ComputePsi(cISize, nyP, kxP, kyP, wC, psiC)
   END SUBROUTINE SetTurbulentVelocity

   !> Subroutine to set a velocity field fitting a desired energy spectrum.
   !!
   !! This subroutine uses a modified form of Rogallo's method to construct a
   !! solenoidal complex velocity field.
   !!
   !> @param[in] rank MPI rank of this process.
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction on this process.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] rISize Size of real data in the i direction.
   !> @param[in] cISize Size of complex data in the i direction.
   !> @param[in] kxP Array of wavenumbers in the x direction on this process.
   !> @param[in] kyP Array of wavenumbers in the y direction on this process.
   !> @param[in] kxLG Lowest x wavenumber in the total simulation.
   !> @param[in] kyLG Lowest y wavenumber in the total simulation.
   !> @param[in] kxMG Maximum x wavenumber in the total simulation.
   !> @param[in] kyMG Maximum y wavenumber in the total simulation.
   !> @param[in] k0 Desired peak wavenumber in the spectrum.
   !> @param[in] urms Desired value of urms for the flow.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   SUBROUTINE FourierVelocity(rank, nxG, nyG, nxP, nyP, rISize, cISize, &
                              kxP, kyP, kxLG, kyLG, kxMG, kyMG, k0, urms, &
                              uC, vC)
      ! Required modules.
      USE MPI
      USE Parameters_m,ONLY: PI, II
      USE Analysis_m,ONLY: BinarySearch
      USE Spectral_m,ONLY: Truncate
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rank, nxG, nyG, nxP, nyP, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLG, kyLG, kxMG, kyMG
      REAL(KIND=RWPF),INTENT(IN) :: k0, urms
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC
      ! Local variables.
      ! Width of the bins.
      REAL(KIND=RWPF),PARAMETER :: dk = 1.0_RWPC
      ! Integer versions of the kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: kxI, kyI
      ! Real versions of the kx and ky wavenumbers.
      REAL(KIND=RWPF) :: kxR, kyR
      ! Magnitude of the wavenumber vector.
      REAL(KIND=RWPF) :: kMag
      ! Twice the energy content of a Fourier mode.
      REAL(KIND=RWPF) :: E
      ! Random phase given to each mode.
      REAL(KIND=RWPF) :: theta
      ! Complex coefficient for a given velocity mode.
      COMPLEX(KIND=CWPC) :: a
      ! Maximum wavenumber in the simulation.
      REAL(KIND=RWPF) :: kMax
      ! Maximum integer wavenumber bin we will use.
      INTEGER(KIND=IWPF) :: kIntMax
      ! Number of bins being used.
      INTEGER(KIND=IWPF) :: numBin
      ! Index for a mode in the bin array.
      INTEGER(KIND=IWPF) :: indx
      ! Array to count how many modes contribute to a bin.
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: binCntL
      ! Array to reduce the bin count to the root process.
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: binCntG
      ! Number of bin edges.
      INTEGER(KIND=IWPF) :: numEdges
      ! Array of bin edges for the binning routine.
      REAL(KIND=RWPF),DIMENSION(:),ALLOCATABLE :: binEdges
      ! Number to increase bin count by for a given mode.
      INTEGER(KIND=IWPF) :: cnt
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! Loop over all modes and determine their velocity magnitudes.
      DO j = 1, nyP
         !
         ! The wavenumber in the y direction.
         kyI = kyP(j)
         kyR = REAL(kyP(j), RWPF)
         DO i = 1, cISize
            !
            ! The wavenumber in the x direction.
            kxI = kxP(i)
            kxR = REAL(kxP(i), RWPF)
            !
            ! The magnitude of the wavenumber vector.
            kMag = SQRT(kxR**2 + kyR**2)
            !
            ! Calculate the desired value of the energy spectrum for this mode.
            E = Schumann(kMag, k0, urms)
            !
            ! Calculate the random phase given to this mode.
            CALL RANDOM_NUMBER(theta)
            theta = 2.0_RWPF*PI*theta
            !
            ! Set the complex coefficient for the velocity field.
            a = SQRT(REAL(E*dk, RWPC))*EXP(II*REAL(theta, RWPC))
            !
            ! Set the complex velocity components.
            IF ((kxI == 0_IWPF) .AND. (kyI == 0_IWPF)) THEN
               uC(i,j) = (0.0_RWPC, 0.0_RWPC)
               vC(i,j) = (0.0_RWPC, 0.0_RWPC)
            ELSE
               uC(i,j) = -1.0_RWPC*a*REAL(kyR, RWPC)/REAL(kMag, RWPC)
               vC(i,j) = a*REAL(kxR, RWPC)/REAL(kMag, RWPC)
            END IF
         END DO
      END DO
      !
      ! At this point we still have to normalize the velocity modes by the
      ! number that contributes to each bin in the energy spectrum.
      !
      ! The maximum wavenumber in the simulation.
      kMax = SQRT(REAL(kxMG*kxMG + kyMG*kyMG, RWPF))
      !
      ! The maximum integer wavenumber bin center we will use.
      kIntMax = INT(kMax + 1.0_RWPC, IWPF)
      !
      ! Number of bins being used.
      numBin = kIntMax + 1_IWPF
      !
      ! Allocate memory for the bin counters.
      ALLOCATE(BinCntL(numBin))
      BinCntL(:) = 0_IWPF
      ALLOCATE(BinCntG(numBin))
      BinCntG(:) = 0_IWPF
      !
      ! Number of bin edges there will be.
      numEdges = numBin + 1_IWPF
      !
      ! Fill in the wavenumber bin edges, which are used for the binning code.
      ALLOCATE(binEdges(numEdges))
      binEdges(1) = 0.0_RWPF
      DO i = 1, kIntMax
         binEdges(i+1) = 0.5_RWPF + REAL(i-1, RWPF)*1.0_RWPF
      END DO
      binEdges(numEdges) = REAL(kIntMax, RWPF) + 0.1_RWPF
      !
      ! Loop over all wavenumbers and count how many are in each bin.
      DO j = 1, nyP
         !
         ! The wavenumber in the y direction.
         kyI = kyP(j)
         kyR = REAL(kyP(j), RWPF)
         DO i = 1, cISize
            !
            ! The wavenumber in the x direction.
            kxI = kxP(i)
            kxR = REAL(kxP(i), RWPF)
            !
            ! The magnitude of the wavenumber vector.
            kMag = SQRT(kxR**2 + kyR**2)
            !
            ! Determine which bin the mode is in.
            CALL BinarySearch(numEdges, binEdges, kMag, indx)
            !
            ! When determining how many go in each bin, we must take into
            ! account what data is stored by FFTW. Because of Hermitian
            ! symmetry, when kx is between (1,Nx/2-1), we double the mode's
            ! contribution to take into account the modes not stored in memory.
            IF ((kxI == 0_IWPF) .OR. (kxI == kxMg)) THEN
               cnt = 1_IWPF
            ELSE
               cnt = 2_IWPF
            END IF
            !
            ! Increase the bin counter for this bin.
            binCntL(indx) = binCntL(indx) + cnt
         END DO
      END DO

      ! Reduce the bin counters to all processes.
      CALL MPI_ALLREDUCE(binCntL, binCntG, numBin, MPI_INTEGER, MPI_SUM, &
                         MPI_COMM_WORLD, ierr)

      ! Each process then goes through and normalizes their Fourier modes by
      ! the square root of the number of modes in each bin.
      DO j = 1, nyP
         !
         ! The wavenumber in the y direction.
         kyI = kyP(j)
         kyR = REAL(kyP(j), RWPF)
         DO i = 1, cISize
            !
            ! The wavenumber in the x direction.
            kxI = kxP(i)
            kxR = REAL(kxP(i), RWPF)
            !
            ! The magnitude of the wavenumber vector.
            kMag = SQRT(kxR**2 + kyR**2)
            !
            ! Determine which bin the mode is in.
            CALL BinarySearch(numEdges, binEdges, kMag, indx)
            !
            ! Normalize the Fourier mode.
            uC(i,j) = uC(i,j)/SQRT(REAL(binCntG(indx), RWPC))
            vC(i,j) = vC(i,j)/SQRT(REAL(binCntG(indx), RWPC))
         END DO
      END DO

      ! Truncate undesired Fourier modes.
      CALL Truncate(cISize, nyP, uC)
      CALL Truncate(cISize, nyP, vC)

      ! Free temporary memory.
      DEALLOCATE(BinCntL)
      DEALLOCATE(binEdges)
      DEALLOCATE(BinCntG)
   END SUBROUTINE FourierVelocity

   !> Energy spectrum attributed to Schumann.
   !!
   !! Need a citation.
   !!
   !> @param[in] k Wavenumber magnitude.
   !> @param[in] k0 Desired wavenumber peak for the spectrum.
   !> @param[in] urms Desired RMS velocity.
   !> @param[out] E Energy of the mode with wavenumber magnitude k.
   REAL(KIND=RWPF) FUNCTION Schumann(k, k0, urms) RESULT(E)
      IMPLICIT NONE
      ! Calling arguments.
      REAL(KIND=RWPF),INTENT(IN) :: k, k0, urms

      ! Calculate the energy following the schumann spectrum.
      E = 1.5_RWPF*urms**2*k/k0*EXP(-k/k0)
   END FUNCTION Schumann

END MODULE SetIC_m

