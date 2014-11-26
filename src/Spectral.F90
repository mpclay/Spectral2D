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

   ! Information pertaining to dealiasing.
   !
   !> Structure used for wavenumber truncation.
   !!
   !! We store the minimum and maximum indices in each direction for each block
   !! of contiguous memory to be trunctated in wavenumber space. Each mode will
   !! have a Fourier coefficient of zero in the region defined by this object.
   TYPE,PRIVATE :: Truncation_t
      !> Lower index in the x direction.
      INTEGER(KIND=IWPF) :: kxL
      !> Upper index in the x direction.
      INTEGER(KIND=IWPF) :: kxU
      !> Lower index in the y direction.
      INTEGER(KIND=IWPF) :: kyL
      !> Upper index in the y direction.
      INTEGER(KIND=IWPF) :: kyU
   END TYPE
   !
   !> Dealiasing technique set by the calling code.
   INTEGER(KIND=IWPF),PRIVATE :: dealias
   !> Number of contiguous memory chunks in spectral space to be truncated.
   INTEGER(KIND=IWPF),PRIVATE :: numTrunc
   !> Data objects to store information about truncation regions.
   TYPE(Truncation_t),DIMENSION(:),ALLOCATABLE,PRIVATE :: trunc

   ! Variables for working with the FFTW library.
   !
   !> FFT plan for computing the forward (R to C) DFT.
   TYPE(C_PTR),PRIVATE :: r2cPlan
   !> FFT plan for computing the reverse (C to R) DFT.
   TYPE(C_PTR),PRIVATE :: c2rPlan

   ! Module procedures.
   PUBLIC :: GetGridSize, SpectralSetup, SpectralFinalize
   PUBLIC :: ComputeRHS, ComputePsi, ComputeVelocity, ComputeVorticity
   PUBLIC :: TransformR2C, TransformC2R
   PUBLIC :: DuDx, DuDy

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
   !> @param[out] rISize Dimension of i extent for real data arrays.
   !> @param[out] cISize Dimension of i extent for complex data arrays.
   !> @param[in,out] kxG Global list of wavenumbers in the x direction.
   !> @param[in,out] kyG Global list of wavenumbers in the y direction.
   SUBROUTINE GetGridSize(nx, ny, dealias_, nxG, nyG, &
                          kxLG, kyLG, kxMG, kyMG, &
                          kxLU, kyLU, kxMU, kyMU, &
                          rISize, cISize, kxG, kyG)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nx, ny, dealias_
      INTEGER(KIND=IWPF),INTENT(OUT) :: nxG, nyG
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLG, kyLG, kxMG, kyMG
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLU, kyLU, kxMU, kyMU
      INTEGER(KIND=IWPF),INTENT(OUT) :: rISize, cISize
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
            ! Set the min/max wavenumber actually being used.
            kxLU = 0_IWPF
            kyLU = -ny/2_IWPF + 1_IWPF
            kxMU = nx/2_IWPF
            kyMU = ny/2_IWPF
            !
            ! Do not allocate any truncation objects.
            numTrunc = 0_IWPF
         CASE (DEALIAS_3_2_RULE)
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
      !        of the data is real, only nxG/2+1 complex numbers are required
      !        to store the Fourier coefficients due to symmetry.
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
   !> @param[in] rISize Size of memory in i direction for Qr.
   !> @param[in] cISize Size of memory in j direction for Qc.
   !> @param[in] nVar Number of variables in the Q array.
   !> @param[in] nStr Number of time storage locations in the Q array.
   !> @param[in] Qr Real cast of data array.
   !> @param[in] Qc Complex cast of data array.
   SUBROUTINE SpectralSetup(nxG, nyG, rISize, cISize, nVar, nStr, Qr, Qc)
      ! Required modules.
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, rISize, cISize, nVar, nStr
      REAL(KIND=RWPC),DIMENSION(rISize,nyG,nVar,nStr),INTENT(INOUT) :: Qr
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG,nVar,nStr),INTENT(INOUT) :: Qc

      ! The plan for the forward transform (R to C). Reverse ordering for C.
      r2cPlan = FFTW_MPI_PLAN_DFT_R2C_2D(INT(nyG, C_INTPTR_T), &
                                         INT(nxG, C_INTPTR_T), &
                                         Qr(:,:,1,1), Qc(:,:,1,1), &
                                         MPI_COMM_WORLD, FFTW_ESTIMATE)
      !
      ! The plan for the reverse transform (C to R). Reverse ordering for C.
      c2rPlan = FFTW_MPI_PLAN_DFT_C2R_2D(INT(nyG, C_INTPTR_T), &
                                         INT(nxG, C_INTPTR_T), &
                                         Qc(:,:,1,1), Qr(:,:,1,1), &
                                         MPI_COMM_WORLD, FFTW_ESTIMATE)
   END SUBROUTINE SpectralSetup

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
   !! After the vorticity has been updated by the time integrator, the stream
   !! function is determined by the Poisson equation (not in this routine):
   !!
   !!    psi = omega/k**2
   !!
   !> @param[in] nxG Number of grid points in the x direction.
   !> @param[in] nyG Number of grid points in the y direction.
   !> @param[in] rISize Size of real data in the i direction.
   !> @param[in] cISize Size of complex data in the i direction.
   !> @param[in] kxG Array of wavenumbers in the x direction.
   !> @param[in] kyG Array of wavenumbers in the y direction.
   !> @param[in] nu The physical viscosity.
   !> @param[in] nVar Number of variables in Q.
   !> @param[in] uC Complex cast of u velocity array.
   !> @param[in] uR Real cast of u velocity array.
   !> @param[in] vC Complex cast of v velocity array.
   !> @param[in] vR Real cast of v velocity array.
   !> @param[in] Qc Complex cast of Q array at this time step.
   !> @param[in] Qr Real cast of Q array at this time step.
   !> @param[in] dQ Increment for the Q array.
   SUBROUTINE ComputeRHS(nxG, nyG, rISize, cISize, &
                         kxG, kyG, nu, nVar, &
                         uC, uR, vC, vR, Qc, Qr, dQ)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, rISize, cISize, nVar
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      REAL(KIND=RWPC),INTENT(IN) :: nu
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyG),INTENT(INOUT) :: uR, vR
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG,2),INTENT(INOUT) :: Qc
      REAL(KIND=RWPC),DIMENSION(rISize,nyG,2),INTENT(INOUT) :: Qr
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(OUT) :: dQ
      ! Local variables.
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: i, j
      ! Wavenumber in the x direction.
      REAL(KIND=RWPC) :: kx
      ! Wavenumber in the y direction.
      REAL(KIND=RWPC) :: ky
      ! Wavenumbers being masked for differentiation in the x direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kxMask
      ! Wavenumbers being masked for differentiation in the y direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kyMask

      ! Steps to form the nonlinear product term.
      !
      ! 1. Differentiate the streamfunction to get u and v velocity.
      !
      !     A. u = dPsi/dy (physical) = i*k2*Psi (spectral)
      !     B. v = -dPsi/dx (physical) = -i*k1*Psi (spectral)
      !
      CALL ComputeVelocity(nxG, nyG, rISize, cISize, kxG, kyG, &
                           uC, uR, vC, vR, Qc(:,:,2), .FALSE.)
      !
      ! 2. Explicitly truncate higher frequency modes before inversion.
      !
      !    NOTE: This is already done for u and v, since Q is explicitly
      !    truncated at the end of each RK stage.
      !
      ! 3. Invert the vorticity and velocity components to physical space.
      !
      CALL TransformC2R(rISize, cISize, nyG, uC, uR)
      CALL TransformC2R(rISize, cISize, nyG, vC, vR)
      CALL TransformC2R(rISize, cISize, nyG, Qc(:,:,1), Qr(:,:,1))
      !
      ! 4. Form the nonlinear product terms w*u and w*v in physical space.
      !
      !    NOTE: These scratch results for the nonlinear terms are stored in the
      !    u and v arrays themselves.
      !
      DO j = 1, nyG
         DO i = 1, nxG
            uR(i,j) = uR(i,j)*Qr(i,j,1)
            vR(i,j) = vR(i,j)*Qr(i,j,1)
         END DO
      END DO
      !
      ! 5. Transform the nonlinear terms and vorticity back to spectral space.
      !
      CALL TransformR2C(nxG, nyG, rISize, cISize, uR, uC)
      CALL TransformR2C(nxG, nyG, rISize, cISize, vR, vC)
      CALL TransformR2C(nxG, nyG, rISize, cISize, Qr(:,:,1), Qc(:,:,1))
      !
      ! 6. Truncate the high wavenumber spectral content from the signal.
      !
      ! 7. Perform differentation of the nonlinear term in spectral space.
      !
      !    NOTE: After this step, the nonlinear terms are ready to go.
      !
      kxMask(1) = nxG/2_IWPF + 1_IWPF
      kyMask(1) = nyG/2_IWPF + 1_IWPF
      CALL DuDx(cISize, nyG, kxG, kyG, 1_IWPF, kxMask, uC, uC)
      CALL DuDy(cISize, nyG, kxG, kyG, 1_IWPF, kyMask, vC, vC)

      ! Form the RHS of the ODE system in spectral space.
      !
      ! 1. Loop over all wavennumbers for the vorticity and increment dQ.
      !
      DO j = 1, nyG
         ! The wavenumber in the y direction.
         ky = REAL(kyG(j), RWPC)
         DO i = 1, cISize
            ! The wavenumber in the x direction.
            kx = REAL(kxG(i), RWPC)
            !
            ! Increment dQ.
            dQ(i,j) = -uC(i,j) - vC(i,j) - nu*(kx**2 + ky**2)*Qc(i,j,1)
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
   !> @param[in] nyG Total number of y grid points in the simulation.
   !> @param[in] kxG Array of wavenumbers for the x direction.
   !> @param[in] kyG Array of wavenumbers for the y direction.
   !> @param[in] Qw Current state vector for vorticity.
   !> @param[out] Qp Current state vector for the streamfunction.
   SUBROUTINE ComputePsi(cISize, nyG, kxG, kyG, Qw, Qp)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyG
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(IN) :: Qw
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(OUT) :: Qp
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
      DO j = 2, nyG
         ! The wavenumber in the y direction.
         ky = REAL(kyG(j), RWPC)
         DO i = 2, cISize
            ! The wavenumber in the x direction.
            kx = REAL(kxG(i), RWPC)
            !
            ! Determine the streamfunction.
            Qp(i,j) = Qw(i,j)/(kx**2 + ky**2)
         END DO
      END DO
      !
      ! 2. Loop over all nonzero x wavenumbers that have ky=0.
      DO i = 2, cISize
         ! The wavenumber in the x direction.
         kx = REAL(kxG(i), RWPC)
         !
         ! Determine the streamfunction.
         Qp(i,1) = Qw(i,1)/kx**2
      END DO
      !
      ! 3. Loop over all nonzero y wavenumbers that have kx=0.
      DO j = 2, nyG
         ! The wavenumber in the y direction.
         ky = REAL(kyG(j), RWPC)
         !
         ! Determine the streamfunction.
         Qp(1,j) = Qw(1,j)/ky**2
      END DO
      !
      ! 4. Explicitly set the zero mode to zero, since this is indeterminant.
      Qp(1,1) = (0.0_RWPC, 0.0_RWPC)
   END SUBROUTINE ComputePsi

   !> Subroutine to calculate the velocity from the streamfunction.
   !!
   !> @param[in] nxG Number of grid points in the x direction.
   !> @param[in] nyG Number of grid points in the y direction.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] kxG Array of wavenumbers in the x direction.
   !> @param[in] kyG Array of wavenumbers in the y direction.
   !> @param[in] uC Complex cast of u velocity array.
   !> @param[in] uR Real cast of u velocity array.
   !> @param[in] vC Complex cast of v velocity array.
   !> @param[in] vR Real cast of v velocity array.
   !> @param[in] PsiC Complex cast of psi array.
   !> @param[in] trans Whether or not to transform velocity from C to R.
   SUBROUTINE ComputeVelocity(nxG, nyG, rISize, cISize, kxG, kyG, &
                              uC, uR, vC, vR, PsiC, trans)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyG),INTENT(INOUT) :: uR, vR
      COMPLEX(CWPC),DIMENSION(cISize,nyG),INTENT(INOUT) :: PsiC
      LOGICAL,INTENT(IN) :: trans
      ! Local variables.
      ! Masked wavenumber for differentiation in the y direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kyMask
      ! Masked wavenumber for differentiation in the x direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kxMask

      ! 1. Getting the u velocity.
      !
      ! Zero out the working array.
      uC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Differentiate the streamfunction to get the u velocity.
      kyMask(1) = nyG/2_IWPF + 1_IWPF
      CALL DuDy(cISize, nyG, kxG, kyG, 1_IWPF, kyMask, PsiC, uC)

      ! 2. Getting the v velocity.
      !
      ! Zero out the working array.
      vC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Differentiate the streamfunction to get the v velocity.
      kxMask(1) = nxG/2_IWPF + 1_IWPF
      CALL DuDx(cISize, nyG, kxG, kyG, 1_IWPF, kxMask, PsiC, vC)
      vC(:,:) = -1.0_RWPC*vC(:,:)

      ! If requested, transform the velocities back to physical space.
      IF (trans) THEN
         CALL TransformC2R(rISize, cISize, nyG, uC, uR)
         CALL TransformC2R(rISize, cISize, nyG, vC, vR)
      END IF
   END SUBROUTINE ComputeVelocity

   !> Compute the vorticity in spectral space from the velocity.
   !!
   !! NOTE: The uC and vC arrays are destroyed by this operation.
   !!
   !> @param[in] nxG Number of grid points in the x direction.
   !> @param[in] nyG Number of grid points in the y direction.
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] kxG Array of wavenumbers for the x direction.
   !> @param[in] kyG Array of wavenumbers for the y direction.
   !> @param[in] uC Complex u velocity array.
   !> @param[in] vC Complex v velocity array.
   !> @param[in,out] wC Complex vorticity array.
   SUBROUTINE ComputeVorticity(nxG, nyG, cISize, kxG, kyG, uC, vC, wC)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(INOUT) :: uC, vC, wC
      ! Local variables.
      ! Masked wavenumbers in the x direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kxMask
      ! Masked wavenumbers in the y direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kyMask
      ! Looping indices for the kx and ky directions.
      INTEGER(KIND=IWPF) :: i, j

      ! Differentiate the v velocity with respect to x.
      kxMask(1) = nxG/2_IWPF + 1_IWPF
      CALL DuDx(cISize, nyG, kxG, kyG, 1_IWPF, kxMask, vC, vC)
      !
      ! Differentiate the u velocity with respect to y.
      kyMask(1) = nyG/2_IWPF + 1_IWPF
      CALL DuDy(cISize, nyG, kxG, kyG, 1_IWPF, kyMask, uC, uC)
      !
      ! Form the complex vorticity.
      DO j = 1, nyG
         DO i = 1, cISize
            wC(i,j) = vC(i,j) - uC(i,j)
         END DO
      END DO
   END SUBROUTINE ComputeVorticity

   !> Routine to differentiate a signal in the x direction.
   !!
   !! In this routine we multiply a signal by i*kx, where kx is the wavenumber
   !! in the x direction. Wavenumbers in the mask array will have their values
   !! explicitly set to zero after differentiation. These will correspond to all
   !! i indices for each j index in the mask array.
   !!
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] nyG Total number of y grid points in the simulation.
   !> @param[in] kxG Array of wavenumbers for the x direction.
   !> @param[in] kyG Array of wavenumbers for the y direction.
   !> @param[in] numMask Number of y wavenumbers being masked.
   !> @param[in] kxMask Array of wavenumbers being masked.
   !> @param[in] sigIn Signal being differentiated.
   !> @param[out] sigOut Differentiated signal
   SUBROUTINE DuDx(cISize, nyG, kxG, kyG, numMask, kxMask, sigIn, sigOut)
      ! Required modules.
      USE Parameters_m,ONLY: II
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyG, numMask
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      INTEGER(KIND=IWPF),DIMENSION(numMask),INTENT(IN) :: kxMask
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(IN) :: sigIn
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(OUT) :: sigOut
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Wavenumber in the x direction.
      REAL(KIND=RWPC) :: kx
      ! Which wavenumber index is being masked.
      INTEGER(KIND=IWPF) :: maskInd

      ! Loop over memory in contiguous fashion and differentiate.
      DO j = 1, nyG
         DO i = 1, cISize
            ! Wavenumber in the x direction.
            kx = REAL(kxG(i), RWPC)
            !
            ! Differentiate the signal in x.
            sigOut(i,j) = II*kx*sigIn(i,j)
         END DO
      END DO

      ! Zero out undesired modes.
      DO i = 1, numMask
         ! The wavenumber being masked.
         maskInd = kxMask(i)
         !
         ! Zero out all modes with this kx wavenumber.
         sigOut(maskInd,:) = (0.0_RWPC, 0.0_RWPC)
      END DO
   END SUBROUTINE DuDx

   !> Routine to differentiate a signal in the y direction.
   !!
   !! In this routine we multiply a signal by i*ky, where ky is the wavenumber
   !! in the y direction. Wavenumbers in the mask array will have their values
   !! explicitly set to zero after differentiation. These will correspond to all
   !! i indices for each j index in the mask array.
   !!
   !> @param[in] cISize Size of the i dimension for complex arrays.
   !> @param[in] nyG Total number of y grid points in the simulation.
   !> @param[in] kxG Array of wavenumbers for the x direction.
   !> @param[in] kyG Array of wavenumbers for the y direction.
   !> @param[in] numMask Number of y wavenumbers being masked.
   !> @param[in] kyMask Array of wavenumbers being masked.
   !> @param[in] sigIn Signal being differentiated.
   !> @param[out] sigOut Differentiated signal
   SUBROUTINE DuDy(cISize, nyG, kxG, kyG, numMask, kyMask, sigIn, sigOut)
      ! Required modules.
      USE Parameters_m,ONLY: II
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: cISize, nyG, numMask
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      INTEGER(KIND=IWPF),DIMENSION(numMask),INTENT(IN) :: kyMask
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(IN) :: sigIn
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(OUT) :: sigOut
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Wavenumber in the y direction.
      REAL(KIND=RWPC) :: ky
      ! Which wavenumber index is being masked.
      INTEGER(KIND=IWPF) :: maskInd

      ! Loop over memory in contiguous fashion and differentiate.
      DO j = 1, nyG
         ! The wavenumber in the y direction.
         ky = REAL(kyG(j), RWPC)
         DO i = 1, cISize
            ! Differentiate the signal in y.
            sigOut(i,j) = II*ky*sigIn(i,j)
         END DO
      END DO

      ! Zero out undesired modes.
      DO i = 1, numMask
         ! The wavenumber being masked.
         maskInd = kyMask(i)
         !
         ! Zero out all modes with this ky wavenumber.
         sigOut(:,maskInd) = (0.0_RWPC, 0.0_RWPC)
      END DO
   END SUBROUTINE DuDy

   !> Subroutine to transform a signal from physical to spectral space.
   !!
   !! This is primarily intended for other modules that need to do transforms.
   !!
   !> @param[in] nxW Number of grid points in the x direction.
   !> @param[in] nyW Number of grid points in the y direction.
   !> @param[in] rISize Size of the i dimension in real data arrays.
   !> @param[in] cISize Size of the i dimension in complex data arrays.
   !> @param[in,out] rData Real data to be transformed.
   !> @param[in,out] cData Complex data array to accept transform.
   SUBROUTINE TransformR2C(nxW, nyW, rISize, cISize, rData, cData)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxW, nyW, rISize, cISize
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(INOUT) :: rData
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyW),INTENT(INOUT) :: cData
      ! Local variables.
      ! Normalization factor for the forward FFT.
      REAL(KIND=RWPC) :: norm

      ! Perform the forward transform.
      CALL FFTW_MPI_EXECUTE_DFT_R2C(r2cPlan, rData, cData)
      norm = REAL(nxW*nyW, RWPC)
      cData(:,:) = cData(:,:)/norm
   END SUBROUTINE TransformR2C

   !> Subroutine to transform a signal from spectral to physical space.
   !!
   !! This is primarily intended for other modules that need to do transforms.
   !!
   !> @param[in] rISize Size of the i dimension in real data arrays.
   !> @param[in] cISize Size of the i dimension in complex data arrays.
   !> @param[in] nyW Number of grid points in the j direction.
   !> @param[in,out] cData Complex data array to accept transform.
   !> @param[in,out] rData Real data to be transformed.
   SUBROUTINE TransformC2R(rISize, cISize, nyW, cData, rData)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize, nyW
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyW),INTENT(INOUT) :: cData
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(INOUT) :: rData

      ! Perform the forward transform.
      CALL FFTW_MPI_EXECUTE_DFT_C2R(c2rPlan, cData, rData)
   END SUBROUTINE TransformC2R

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

