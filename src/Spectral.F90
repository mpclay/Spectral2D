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
!> @brief Module to perform spatial operations in spectral space.
!!
!! This module provides many routines to calculate quantities in physical and
!! spectral space, e.g., calculating velocities from the streamfunction, the
!! streamfunction from vorticity, etc. Aside from these operations, the module
!! also calculates the RHS of the governing equations for the time integrator.
!! All linear spatial derivatives are evaluated spectrally. The nonlinear terms
!! in the governing equations are evaluated with a pseudo-spectral approach.
MODULE Spectral_m

   ! Required modules.
   USE ISO_C_BINDING
   USE Parameters_m,ONLY: IWPF, RWPF, IWPC, RWPC, CWPC

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

   !> Structure used for wavenumber truncation.
   !!
   !! We store the minimum and maximum i indices to truncate for each j row.
   TYPE,PRIVATE :: Truncation_t
      !> Lower i index for truncation.
      INTEGER(KIND=IWPF) :: i1
      !> Upper i index for truncation.
      INTEGER(KIND=IWPF) :: i2
   END TYPE
   !
   !> Dealiasing technique set by the calling code.
   INTEGER(KIND=IWPF),PRIVATE :: dealias
   !> Data objects to store information about truncation regions.
   TYPE(Truncation_t),DIMENSION(:),ALLOCATABLE,PRIVATE :: trunc
   !
   ! Variables for masking modes during differentiation.
   !
   !> Whether or not the maximum x wavenumber is on this process for masking.
   LOGICAL,PRIVATE :: xMaskBool
   !> Whether or not the maximum y wavenumber is on this process for masking.
   LOGICAL,PRIVATE :: yMaskBool
   !> Index for the maximum x wavenumber if it is on this process.
   INTEGER(KIND=IWPF),PRIVATE :: xMaskInd
   !> Index for the maximum y wavenumber if it is on this process.
   INTEGER(KIND=IWPF),PRIVATE :: yMaskInd

   ! Variables for working with the FFTW library.
   !
   !> FFT plan for computing the forward (R to C) DFT.
   TYPE(C_PTR),PRIVATE :: r2cPlan
   !> FFT plan for computing the reverse (C to R) DFT.
   TYPE(C_PTR),PRIVATE :: c2rPlan

   ! Module procedures.
   PUBLIC :: GetGridSize
   PUBLIC :: SetupMask, SetupPlans, SetupTruncation, SpectralFinalize
   PUBLIC :: ComputeRHS, ComputePsi, ComputeVelocity, ComputeVorticity
   PUBLIC :: TransformR2C, TransformC2R
   PUBLIC :: DuDx, DuDy
   PUBLIC :: Truncate

CONTAINS

   !> Adjust the working grid size based on the dealiasing technique.
   !!
   !! In this routine we also set module variables that are needed when working
   !! in spectral space. These include the min/max wavenumbers being used in
   !! each direction, and the necessary variables for the truncation scheme
   !! (if one is being used).
   !!
   !> @param[in] nx Number of grid points in the x direction.
   !> @param[in] ny Number of grid points in the y direction.
   !> @param[in] dealias_ Desired dealiasing technique.
   !> @param[out] nxG Working number of grid points in the x direction.
   !> @param[out] nyG Working number of grid points in the y direction.
   !> @param[out] kxLG Absolute smallest x wavenumber.
   !> @param[out] kyLG Absolute smallest y wavenumber.
   !> @param[out] kxMG Absolute largest x wavenumber.
   !> @param[out] kyMG Absolute largest y wavenumber.
   !> @param[out] kxLU Smallest useful x wavenumber.
   !> @param[out] kyLU Smallest useful y wavenumber.
   !> @param[out] kxMU Largest useful x wavenumber.
   !> @param[out] kyMU Largest useful y wavenumber.
   !> @param[out] kMax Maximum wavenumber magnitude used for truncation.
   !> @param[out] rISize Dimension of i extent for real data arrays.
   !> @param[out] cISize Dimension of i extent for complex data arrays.
   !> @param[in,out] kxG Global list of wavenumbers in the x direction.
   !> @param[in,out] kyG Global list of wavenumbers in the y direction.
   SUBROUTINE GetGridSize(nx, ny, dealias_, nxG, nyG, &
                          kxLG, kyLG, kxMG, kyMG, &
                          kxLU, kyLU, kxMU, kyMU, &
                          kMax, rISize, cISize, kxG, kyG)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nx, ny, dealias_
      INTEGER(KIND=IWPF),INTENT(OUT) :: nxG, nyG
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLG, kyLG, kxMG, kyMG
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLU, kyLU, kxMU, kyMU
      INTEGER(KIND=IWPF),INTENT(OUT) :: rISize, cISize
      REAL(KIND=RWPF),INTENT(OUT) :: kMax
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: kxG, kyG
      ! Local variables.
      ! Looping index for wavenumbers in the x and y directions.
      INTEGER(KIND=IWPF) :: i
      ! Used to figure out wavenumber in the y direction.
      INTEGER(KIND=IWPF) :: kytmp

      ! Store the desired dealiasing technique.
      dealias = dealias_
      !
      ! Adjust the working grid size as necessary.
      SELECT CASE (dealias)
         CASE (NO_DEALIASING)
            ! Don't do anything to the grid size.
            nxG = nx
            nyG = ny
            !
            ! Set the min/max wavenumbers on the working grid size.
            kxLG = 0_IWPF
            kyLG = -nyG/2_IWPF + 1_IWPF
            kxMG = nxG/2_IWPF
            kyMG = nyG/2_IWPF
            !
            ! Set the min/max wavenumbers actually being used.
            kxLU = 0_IWPF
            kyLU = -ny/2_IWPF + 1_IWPF
            kxMU = nx/2_IWPF
            kyMU = ny/2_IWPF
            !
            ! Maximum wavenumber magnitude allowed by truncation.
            kMax = REAL(kxMU, RWPF)
         CASE (DEALIAS_3_2_RULE)
            ! The grid is multiplied by a factor of 3/2 in each direction.
            nxG = 3_IWPF*nx/2_IWPF
            nyG = 3_IWPF*ny/2_IWPF
            !
            ! The global min/max wavenumbers are controlled by nxG/nyG.
            kxLG = 0_IWPF
            kyLG = -nyG/2_IWPF + 1_IWPF
            kxMG = nxG/2_IWPF
            kyMG = nyG/2_IWPF
            !
            ! For the 3/2 rule, the max/min useful wavenumbers are determined
            ! by the initial grid size.
            kxLU = 0_IWPF
            kyLU = -ny/2_IWPF + 1_IWPF
            kxMU = nx/2_IWPF
            kyMU = ny/2_IWPF
            !
            ! Maximum wavenumber magnitude allowed by truncation.
            kMax = REAL(kxMU, RWPF)
         CASE (DEALIAS_2_3_RULE)
         CASE (DEALIAS_GRID_SHIFT)
         CASE DEFAULT
            ! Default will be no dealiasing, but with a warning.
            nxG = nx
            nyG = ny
      END SELECT
      !
      ! Size of i dimension for in-place transforms.
      rISize = 2_IWPF*(nxG/2_IWPF + 1_IWPF)
      cISize = nxG/2_IWPF + 1_IWPF
      !
      ! Allocate and fill in the global wavenumber arrays. These MUST correspond
      ! to the way the Fourier coefficients are stored in the working arrays in
      ! the FFTW routines!
      !
      !     1. The x direction is the first direction to get FFTd. Because all
      !        of the data is real, only nxG/2+1 complex numbers are used to
      !        store the Fourier coefficients due to symmetry.
      !     2. When the FFT is taken in the y direction, all numbers are complex
      !        (from the x transform), so we have full storage of nyG complex
      !        Fourier coefficients.
      !
      !  - X wavenumbers.
      ALLOCATE(kxG(cISize))
      DO i = 1, cISize
         kxG(i) = i - 1_IWPF
      END DO
      !
      ! - Y wavenumbers.
      ALLOCATE(kyG(nyG))
      DO i = 1, nyG
         kytmp = i - 1_IWPF
         IF (kytmp > nyG/2_IWPF) THEN
            kyG(i) = kytmp - nyG
         ELSE
            kyG(i) = kytmp
         END IF
      END DO
   END SUBROUTINE GetGridSize

   !> Subroutine to set up the masked wavenumbers for differentiation.
   !!
   !> @param[in] cISize Number of wavenumbers in i for complex arrays.
   !> @param[in] nyP Number of grid points in y for this process.
   !> @param[in] kxP Array of x wavenumbers for this process.
   !> @param[in] kyP Array of y wavenumbers for this process.
   !> @param[in] kxMG Maximum total wavenumber in x.
   !> @param[in] kyMG Maximum total wavenumber in y.
   SUBROUTINE SetupMask(cISize, nyP, kxP, kyP, kxMG, kyMG)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyP, kxMG, kyMG
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j

      ! Check x wavenumbers to see if kxMG is present.
      xMaskBool = .FALSE.
      xMaskInd = -1_IWPF
      iLoop: DO i = 1, cISize
         IF (kxP(i) == kxMG) THEN
            xMaskBool = .TRUE.
            xMaskInd = i
            EXIT iLoop
         END IF
      END DO iLoop
      !
      ! Check y wavenumbers to see if kyMG is present.
      yMaskBool = .FALSE.
      yMaskInd = -1_IWPF
      jLoop: DO j = 1, nyP
         IF (kyP(j) == kyMG) THEN
            yMaskBool = .TRUE.
            yMaskInd = j
            EXIT jLoop
         END IF
      END DO jLoop
   END SUBROUTINE SetupMask

   !> Set up the truncation scheme being used.
   !!
   !! For each y row of data on a given process, we identify the starting and
   !! ending wavenumbers that will be truncated.
   !!
   !> @param[in] cISize Size of i dimension for complex arrays.
   !> @param[in] Number of points in the y direction for this process.
   !> @param[in] kMax Maximum wavenumber magnitude allowed by truncation.
   !> @param[in] kxP Array of x wavenumbers for this process.
   !> @param[in] kyP Array of y wavenumbers for this process.
   SUBROUTINE SetupTruncation(cISize, nyP, kMax, kxP, kyP)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISIze, nyP
      REAL(KIND=RWPF),INTENT(IN) :: kMax
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      ! Local variables.
      ! Magnitude of the local wavenumber.
      REAL(KIND=RWPF) :: kMag
      ! Wavenumbers in the x and y direction.
      REAL(KIND=RWPF) :: kx, ky
      ! Small epsilon to make sure we do not miss modes with kMag = kMax.
      REAL(KIND=RWPF),PARAMETER :: eps = 1.0E-10_RWPF
      ! Looping index for the kx and ky directions.
      INTEGER(KIND=IWPF) :: i, j

      ! Allocate memory for the truncation objects.
      ALLOCATE(trunc(nyP))
      !
      ! For each j row, set the starting and ending wavenumbers for truncation.
      jLoop: DO j = 1, nyP
         !
         ! Wavenumber in the y direction.
         ky = REAL(kyP(j), RWPF)
         iLoop: DO i = 1, cISize
            !
            ! Wavenumber in the x direction.
            kx = REAL(kxP(i), RWPF)
            !
            ! Overall wavenumber magnitude.
            kMag = SQRT(kx**2 + ky**2)
            !
            ! Check to see if this wavenumber is beyond the threshold.
            IF (kMag > kMax + eps) THEN
               trunc(j)%i1 = i
               trunc(j)%i2 = cISize
               EXIT iLoop
            END IF
         END DO iLoop
      END DO jLoop
   END SUBROUTINE SetupTruncation

   !> Initialize the spectral module.
   !!
   !! In this routine we create the FFTW MPI plans that are used to transform
   !! signals to/from spectral space.
   !!
   !! NOTE: When actually making this parallel, we need to figure out what
   !! y dimension the FFTW routine expects (local or global).
   !!
   !> @param[in] nxG Global size of the arrays in the x direction.
   !> @param[in] nyG Global size of the arrays in the y direction.
   !> @param[in] nyP Number of points in the y direction for this process.
   !> @param[in] rISize Size of memory in i direction for Qr.
   !> @param[in] cISize Size of memory in j direction for Qc.
   !> @param[in] dataR Real cast of data array.
   !> @param[in] dataC Complex cast of data array.
   SUBROUTINE SetupPlans(nxG, nyG, nyP, rISize, cISize, dataR, dataC)
      ! Required modules.
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nyP, rISize, cISize
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: dataR
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: dataC

      ! The plan for the forward transform (R to C). Reverse ordering for C.
      r2cPlan = FFTW_MPI_PLAN_DFT_R2C_2D(INT(nyG, C_INTPTR_T), &
                                         INT(nxG, C_INTPTR_T), &
                                         dataR, dataC, &
                                         MPI_COMM_WORLD, FFTW_MEASURE)
      !
      ! The plan for the reverse transform (C to R). Reverse ordering for C.
      c2rPlan = FFTW_MPI_PLAN_DFT_C2R_2D(INT(nyG, C_INTPTR_T), &
                                         INT(nxG, C_INTPTR_T), &
                                         dataC, dataR, &
                                         MPI_COMM_WORLD, FFTW_MEASURE)
      !
      ! Zero out the array in case FFTW wrote garbage values to it.
      dataC(:,:) = (0.0_RWPC, 0.0_RWPC)
   END SUBROUTINE SetupPlans

   !> Subroutine to compute the RHS of the ODE system.
   !!
   !! Here we calculate the RHS of the ODE system for the Fourier vorticity.
   !!
   !!    dw/dt = -G - nu*k**2*w
   !!
   !! where w is the vorticity, G is the nonlinear term, nu is the viscosity,
   !! and psi is the streamfunction. In physical space, the nonlinear term is
   !!
   !!    G = d(w*u_j)/dx_j.
   !!
   !! This term is evaluated pseudo-spectrally in the following way:
   !!
   !!    1. Invert w and u_j from Fourier space to physical space.
   !!    2. Perform the product w*u_j in physical space.
   !!    3. Take the forward transform of w*u_j to go to spectral space.
   !!    4. Perform differentiation in spectral space.
   !!
   !! After the vorticity has been updated by the time integrator, the stream-
   !! function is determined by the Poisson equation (not in this routine):
   !!
   !!    psi = omega/k**2
   !!
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction on this process.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] rISize Size of real data in the i direction.
   !> @param[in] cISize Size of complex data in the i direction.
   !> @param[in] kxP Array of wavenumbers in the x direction on this process.
   !> @param[in] kyP Array of wavenumbers in the y direction on this process.
   !> @param[in] nu The physical viscosity.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] uR Real cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[in,out] vR Real cast of v velocity array.
   !> @param[in,out] wC Complex cast of vorticity array.
   !> @param[in,out] wR Real cast of vorticity array.
   !> @param[in] psiC Complex cast of streamfunction array.
   !> @param[in] psiR Real cast of streamfunction array.
   !> @param[out] dW Increment for the vorticity array.
   SUBROUTINE ComputeRHS(nxG, nyG, nxP, nyP, rISize, cISize, kxP, kyP, nu, &
                         uC, uR, vC, vR, wC, wR, psiC, psiR, dW)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nxP, nyP, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      REAL(KIND=RWPC),INTENT(IN) :: nu
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC, wC
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: uR, vR, wR
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(IN) :: psiC
      REAL(KIND=CWPC),DIMENSION(rISize,nyP),INTENT(IN) :: psiR
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(OUT) :: dW
      ! Local variables.
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: i, j
      ! Wavenumber in the x direction.
      REAL(KIND=RWPC) :: kx
      ! Wavenumber in the y direction.
      REAL(KIND=RWPC) :: ky

      ! Steps to form the nonlinear product term.
      !
      ! 1. Differentiate the streamfunction to get u and v velocity.
      !
      !     A. u = dPsi/dy (physical) = i*k2*Psi (spectral)
      !     B. v = -dPsi/dx (physical) = -i*k1*Psi (spectral)
      !
      CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                           uC, uR, vC, vR, psiC, .FALSE.)
      !
      ! 3. Invert the vorticity and velocity components to physical space.
      !
      !     NOTE: Proper truncation has already been taken care of for this
      !     inversion. At the end of each RK stage w and psi are explicitly
      !     truncated, so u and v are truncated as well.
      !
      CALL TransformC2R(rISize, cISize, nyP, uC, uR)
      CALL TransformC2R(rISize, cISize, nyP, vC, vR)
      CALL TransformC2R(rISize, cISize, nyP, wC, wR)
      !
      ! 4. Form the nonlinear product terms w*u and w*v in physical space.
      !
      !    NOTE: These scratch results for the nonlinear terms are stored in the
      !    u and v arrays themselves.
      !
      DO j = 1, nyP
         DO i = 1, nxP
            uR(i,j) = uR(i,j)*wR(i,j)
            vR(i,j) = vR(i,j)*wR(i,j)
         END DO
      END DO
      !
      ! 5. Transform the nonlinear terms and vorticity back to spectral space.
      !
      CALL TransformR2C(nxG, nyG, nyP, rISize, cISize, uR, uC)
      CALL TransformR2C(nxG, nyG, nyP, rISize, cISize, vR, vC)
      CALL TransformR2C(nxG, nyG, nyP, rISize, cISize, wR, wC)
      !
      ! 6. Truncate the high wavenumber spectral content from the signal.
      !
      CALL Truncate(cISize, nyP, uC)
      CALL Truncate(cISize, nyP, vC)
      !
      ! 7. Perform differentation of the nonlinear term in spectral space.
      !
      !    NOTE: After this step, the nonlinear terms are ready to go.
      !
      CALL DuDx(cISize, nyP, kxP, kyP, uC, uC)
      CALL DuDy(cISize, nyP, kxP, kyP, vC, vC)

      ! Form the RHS of the ODE system in spectral space.
      !
      ! 1. Loop over all wavennumbers for the vorticity and increment dW.
      !
      DO j = 1, nyP
         ! The wavenumber in the y direction.
         ky = REAL(kyP(j), RWPC)
         DO i = 1, cISize
            ! The wavenumber in the x direction.
            kx = REAL(kxP(i), RWPC)
            !
            ! Increment dW.
            dW(i,j) = -uC(i,j) - vC(i,j) - nu*(kx**2 + ky**2)*wC(i,j)
         END DO
      END DO
   END SUBROUTINE ComputeRHS

   !> Routine to solve the Poisson equation for the stream function.
   !!
   !! The stream function is related to the vorticity by:
   !!
   !!    d^2(Psi)/dx^2 + d^2(Psi)/dy^2 = -w
   !!
   !! In spectral space, this is
   !!
   !!    psi = w/k**2
   !!
   !! In this routine we solve the poisson equation in spectral space.
   !!
   !! NOTE: Currently this routine is broken into 4 chunks in order to handle
   !! the singularity at (kx=0,ky=0) in a safe manner. This should be done in
   !! three chunks instead, where all kx for ky >= 1 are used, then all kx>=1
   !! for ky=0, then kx=0 and ky=0. This could prove to be faster as well,
   !! since the largest set of calculations will all be over contiguous memory.
   !!
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] kxP Array of wavenumbers for the x direction.
   !> @param[in] kyP Array of wavenumbers for the y direction.
   !> @param[in] Qw Current state vector for vorticity.
   !> @param[out] Qp Current state vector for the streamfunction.
   SUBROUTINE ComputePsi(cISize, nyP, kxP, kyP, Qw, Qp)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyP
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(IN) :: Qw
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(OUT) :: Qp
      ! Local variables.
      ! Wavenumber in the x direction.
      REAL(KIND=RWPC) :: kx
      ! Wavenumber in the y direction.
      REAL(KIND=RWPC) :: ky
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: i, j

      ! Loop over all useful wavenumbers and compute the streamfunction. To
      ! avoid division by zero, we proceed in four main steps.
      !
      ! 1. Loop over all wavenumbers that have neither kx=0 or ky=0.
      DO j = 2, nyP
         ! The wavenumber in the y direction.
         ky = REAL(kyP(j), RWPC)
         DO i = 2, cISize
            ! The wavenumber in the x direction.
            kx = REAL(kxP(i), RWPC)
            !
            ! Determine the streamfunction.
            Qp(i,j) = Qw(i,j)/(kx**2 + ky**2)
         END DO
      END DO
      !
      ! 2. Loop over all nonzero x wavenumbers that have ky=0.
      DO i = 2, cISize
         ! The wavenumber in the x direction.
         kx = REAL(kxP(i), RWPC)
         !
         ! Determine the streamfunction.
         Qp(i,1) = Qw(i,1)/kx**2
      END DO
      !
      ! 3. Loop over all nonzero y wavenumbers that have kx=0.
      DO j = 2, nyP
         ! The wavenumber in the y direction.
         ky = REAL(kyP(j), RWPC)
         !
         ! Determine the streamfunction.
         Qp(1,j) = Qw(1,j)/ky**2
      END DO
      !
      ! 4. Explicitly set the zero mode to zero, since this is indeterminate.
      Qp(1,1) = (0.0_RWPC, 0.0_RWPC)
   END SUBROUTINE ComputePsi

   !> Subroutine to calculate the velocity from the streamfunction.
   !!
   !> @param[in] nxP Number of grid points in the x direction for this process.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] kxP Array of wavenumbers in the x direction for this process.
   !> @param[in] kyP Array of wavenumbers in the y direction for this process.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] uR Real cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[in,out] vR Real cast of v velocity array.
   !> @param[in] PsiC Complex cast of psi array.
   !> @param[in] trans Whether or not to transform velocity from C to R.
   SUBROUTINE ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                              uC, uR, vC, vR, PsiC, trans)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxP, nyP, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: uR, vR
      COMPLEX(CWPC),DIMENSION(cISize,nyP),INTENT(IN) :: PsiC
      LOGICAL,INTENT(IN) :: trans

      ! 1. Getting the u velocity.
      !
      ! Zero out the working array.
      uC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Differentiate the streamfunction to get the u velocity.
      CALL DuDy(cISize, nyP, kxP, kyP, PsiC, uC)

      ! 2. Getting the v velocity.
      !
      ! Zero out the working array.
      vC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Differentiate the streamfunction to get the v velocity.
      CALL DuDx(cISize, nyP, kxP, kyP, PsiC, vC)
      vC(:,:) = -1.0_RWPC*vC(:,:)

      ! If requested, transform the velocities back to physical space.
      IF (trans) THEN
         CALL TransformC2R(rISize, cISize, nyP, uC, uR)
         CALL TransformC2R(rISize, cISize, nyP, vC, vR)
      END IF
   END SUBROUTINE ComputeVelocity

   !> Compute the vorticity in spectral space from the velocity.
   !!
   !! NOTE: The uC and vC arrays are destroyed by this operation.
   !!
   !> @param[in] nxP Number of grid points in the x direction for this process.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] kxP Array of wavenumbers for the x direction for this process.
   !> @param[in] kyP Array of wavenumbers for the y direction for this process.
   !> @param[in] uC Complex u velocity array.
   !> @param[in] vC Complex v velocity array.
   !> @param[in,out] wC Complex vorticity array.
   SUBROUTINE ComputeVorticity(nxP, nyP, cISize, kxP, kyP, uC, vC, wC)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxP, nyP, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC, wC
      ! Local variables.
      ! Looping indices for the kx and ky directions.
      INTEGER(KIND=IWPF) :: i, j

      ! Differentiate the v velocity with respect to x.
      CALL DuDx(cISize, nyP, kxP, kyP, vC, vC)
      !
      ! Differentiate the u velocity with respect to y.
      CALL DuDy(cISize, nyP, kxP, kyP, uC, uC)
      !
      ! Form the complex vorticity.
      DO j = 1, nyP
         DO i = 1, cISize
            wC(i,j) = vC(i,j) - uC(i,j)
         END DO
      END DO
   END SUBROUTINE ComputeVorticity

   !> Routine to differentiate a signal in the x direction.
   !!
   !! In this routine we multiply a signal by i*kx, where kx is the wavenumber
   !! in the x direction. If the maximum wavenumber in the x direction is
   !! present on this process, that wavenumber will be explicitly masked in
   !! the differentiation process.
   !!
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] nyP Number of y grid points in the simulation for this process.
   !> @param[in] kxP Array of wavenumbers for the x direction for this process. 
   !> @param[in] kyP Array of wavenumbers for the y direction for this process.
   !> @param[in] sigIn Signal being differentiated.
   !> @param[out] sigOut Differentiated signal
   SUBROUTINE DuDx(cISize, nyP, kxP, kyP, sigIn, sigOut)
      ! Required modules.
      USE Parameters_m,ONLY: II
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyP
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(IN) :: sigIn
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(OUT) :: sigOut
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Wavenumber in the x direction.
      REAL(KIND=RWPC) :: kx

      ! Loop over memory in contiguous fashion and differentiate.
      DO j = 1, nyP
         DO i = 1, cISize
            ! Wavenumber in the x direction.
            kx = REAL(kxP(i), RWPC)
            !
            ! Differentiate the signal in x.
            sigOut(i,j) = II*kx*sigIn(i,j)
         END DO
      END DO

      ! Zero out undesired modes if the max wavenumber in x is present.
      IF (xMaskBool) THEN
         sigOut(xMaskInd,:) = (0.0_RWPC, 0.0_RWPC)
      END IF
   END SUBROUTINE DuDx

   !> Routine to differentiate a signal in the y direction.
   !!
   !! In this routine we multiply a signal by i*ky, where ky is the wavenumber
   !! in the y direction. If the maximum wavenumber in the y direction is
   !! present on this process, that wavenumber will be explicitly masked in
   !! the differentiation process.
   !!
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] nyP Number of y grid points in the simulation for this process.
   !> @param[in] kxP Array of wavenumbers for the x direction for this process.
   !> @param[in] kyP Array of wavenumbers for the y direction for this process.
   !> @param[in] sigIn Signal being differentiated.
   !> @param[out] sigOut Differentiated signal
   SUBROUTINE DuDy(cISize, nyP, kxP, kyP, sigIn, sigOut)
      ! Required modules.
      USE Parameters_m,ONLY: II
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyP
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(IN) :: sigIn
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(OUT) :: sigOut
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Wavenumber in the y direction.
      REAL(KIND=RWPC) :: ky

      ! Loop over memory in contiguous fashion and differentiate.
      DO j = 1, nyP
         ! The wavenumber in the y direction.
         ky = REAL(kyP(j), RWPC)
         DO i = 1, cISize
            ! Differentiate the signal in y.
            sigOut(i,j) = II*ky*sigIn(i,j)
         END DO
      END DO

      ! Zero out undesired modes if the max wavenumber in y is present.
      IF (yMaskBool) THEN
         sigOut(:,yMaskInd) = (0.0_RWPC, 0.0_RWPC)
      END IF
   END SUBROUTINE DuDy

   !> Subroutine to transform a signal from physical to spectral space.
   !!
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] rISize Size of the i dimension in real data arrays.
   !> @param[in] cISize Size of the i dimension in complex data arrays.
   !> @param[in,out] rData Real data to be transformed.
   !> @param[in,out] cData Complex data array to accept transform.
   SUBROUTINE TransformR2C(nxG, nyG, nyP, rISize, cISize, rData, cData)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nyP, rISize, cISize
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: rData
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: cData
      ! Local variables.
      ! Normalization factor for the forward FFT.
      REAL(KIND=RWPC) :: norm

      ! Perform the forward transform.
      CALL FFTW_MPI_EXECUTE_DFT_R2C(r2cPlan, rData, cData)
      norm = REAL(nxG*nyG, RWPC)
      cData(:,:) = cData(:,:)/norm
   END SUBROUTINE TransformR2C

   !> Subroutine to transform a signal from spectral to physical space.
   !!
   !> @param[in] rISize Size of the i dimension in real data arrays.
   !> @param[in] cISize Size of the i dimension in complex data arrays.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in,out] cData Complex data array to accept transform.
   !> @param[in,out] rData Real data to be transformed.
   SUBROUTINE TransformC2R(rISize, cISize, nyP, cData, rData)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize, nyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: cData
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: rData

      ! Perform the forward transform.
      CALL FFTW_MPI_EXECUTE_DFT_C2R(c2rPlan, cData, rData)
   END SUBROUTINE TransformC2R

   !> Subroutine to truncate undesired Fourier modes.
   !!
   !> @param[in] cISize Size of the i dimension of complex arrays.
   !> @param[in] nyP Number of y points on this process.
   !> @param[in] dataC Complex data array to truncate.
   SUBROUTINE Truncate(cISize, nyP, dataC)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyP
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: dataC
      ! Local variables.
      ! Looping index for truncation objects.
      INTEGER(KIND=IWPF) :: j

      ! Loop over each truncation object and zero out the complex numbers for
      ! the i range it specifies.
      DO j = 1, nyP
         dataC(trunc(j)%i1:trunc(j)%i2,j) = (0.0_RWPC, 0.0_RWPC)
      END DO
   END SUBROUTINE Truncate

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

