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
   PUBLIC :: ComputeRHS, ComputePsi
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
   !> @param[out] nxW Working number of grid points in the x direction.
   !> @param[out] nyW Working number of grid points in the y direction.
   !> @param[out] kxLT Absolute smallest x wavenumber (truncated).
   !> @param[out] kyLT Absolute smallest y wavenumber (truncated).
   !> @param[out] kxMT Absolute largest x wavenumber (truncated).
   !> @param[out] kyMT Absolute largest y wavenumber (truncated).
   !> @param[out] kxLU Smallest useful x wavenumber.
   !> @param[out] kyLU Smallest useful y wavenumber.
   !> @param[out] kxMU Largest useful x wavenumber.
   !> @param[out] kyMU Largest useful y wavenumber.
   SUBROUTINE GetGridSize(nx, ny, dealias_, nxW, nyW, &
                          kxLT, kyLT, kxMT, kyMT, &
                          kxLU, kyLU, kxMU, kyMU)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nx, ny, dealias_
      INTEGER(KIND=IWPF),INTENT(OUT) :: nxW, nyW
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLT, kyLT, kxMT, kyMT
      INTEGER(KIND=IWPF),INTENT(OUT) :: kxLU, kyLU, kxMU, kyMU

      ! Store the desired dealiasing technique.
      dealias = dealias_
      !
      ! Adjust the working grid size as necessary.
      SELECT CASE (dealias)
         CASE (NO_DEALIASING)
            ! Don't do anything to the grid size.
            nxW = nx
            nyW = ny
            !
            ! Set the min/max wavenumbers on the working grid size.
            kxLT = 0_IWPF
            kyLT = -nyW/2_IWPF + 1_IWPF
            kxMT = nxW/2_IWPF
            kyMT = nyW/2_IWPF
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
            nxW = nx
            nyW = ny
      END SELECT
   END SUBROUTINE GetGridSize

   !> Initialize the spectral module.
   !!
   !! In this routine we create the FFTW MPI plans that are used to transform
   !! signals to/from spectral space.
   !!
   !> @param[in] nxW Working size of the arrays in the x direction.
   !> @param[in] nyW Working size of the arrays in the y direction.
   !> @param[in] rISize Size of memory in i direction for Qr.
   !> @param[in] cISize Size of memory in j direction for Qc.
   !> @param[in] nVar Number of variables in the Q array.
   !> @param[in] nStr Number of time storage locations in the Q array.
   !> @param[in] Qr Real cast of data array.
   !> @param[in] Qc Complex cast of data array.
   SUBROUTINE SpectralSetup(nxW, nyW, rISize, cISize, nVar, nStr, Qr, Qc)
      ! Required modules.
      USE MPI,ONLY: MPI_COMM_WORLD
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxW, nyW, rISize, cISize, nVar, nStr
      REAL(KIND=RWPC),DIMENSION(rISize,nyW,nVar,nStr),INTENT(INOUT) :: Qr
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyW,nVar,nStr),INTENT(INOUT) :: Qc

      ! The plan for the forward transform (R to C). Reverse ordering for C.
      r2cPlan = FFTW_MPI_PLAN_DFT_R2C_2D(INT(nyW, C_INTPTR_T), &
                                         INT(nxW, C_INTPTR_T), &
                                         Qr(:,:,1,1), Qc(:,:,1,1), &
                                         MPI_COMM_WORLD, FFTW_ESTIMATE)
      !
      ! The plan for the reverse transform (C to R). Reverse ordering for C.
      c2rPlan = FFTW_MPI_PLAN_DFT_C2R_2D(INT(nyW, C_INTPTR_T), &
                                         INT(nxW, C_INTPTR_T), &
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
   !> @param[in] rISize Size of real data in the i direction (with padding).
   !> @param[in] nxW Number of useful grid points in the i direction for
   !!            physical space arrays.
   !> @param[in] nyW Size of data arrays in the y direction.
   !> @param[in] kxLT Absolute smallest x wavenumber (truncated).
   !> @param[in] kyLT Absolute smallest y wavenumber (truncated).
   !> @param[in] kxMT Absolute largest x wavenumber (truncated).
   !> @param[in] kyMT Absolute largest y wavenumber (truncated).
   !> @param[in] kxLU Smallest useful x wavenumber.
   !> @param[in] kyLU Smallest useful y wavenumber.
   !> @param[in] kxMU Largest useful x wavenumber.
   !> @param[in] kyMU Largest useful y wavenumber.
   !> @param[in] Q Current state vector.
   !> @param[in] nu The physical viscosity.
   !> @param[in,out] uC Working array for the u velocity complex cast.
   !> @param[in,out] vC Working array for the v velocity complex cast.
   !> @param[in,out] wC Scratch array for the vorticity complex cast.
   !> @param[in,out] uR Working array for the u velocity complex cast.
   !> @param[in,out] vR Working array for the v velocity complex cast.
   !> @param[in,out] wR Scratch array for the vorticity real cast.
   !> @param[out] dQ Increment for w and psi in spectral space.
   SUBROUTINE ComputeRHS(rISize, nxW, nyW, &
                         kxLT, kyLT, kxMT, kyMT, &
                         kxLU, kyLU, kxMU, kyMU, &
                         Q, nu, uC, vC, wC, uR, vR, wR, dQ)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, nxW, nyW
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLT, kyLT, kxMT, kyMT
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLU, kyLU, kxMU, kyMU
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT,2),INTENT(INOUT) :: Q
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(INOUT) :: uC
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(INOUT) :: vC
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(INOUT) :: wC
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(INOUT) :: uR, vR, wR
      COMPLEX(KIND=CWPC),DIMENSION(kxLU:kxMU,kyLU:kxMU),INTENT(OUT) :: dQ
      REAL(KIND=RWPC),INTENT(IN) :: nu
      ! Local variables.
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: kx, ky
      ! Wavenumbers being masked for differentiation in the x direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kxMask
      ! Wavenumbers being masked for differentiation in the y direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kyMask
      ! Looping indices for physical space.
      INTEGER(KIND=IWPF) :: i, j
      ! Normalization factor for the FFT.
      REAL(KIND=RWPC) :: norm

      ! Steps to form the nonlinear product term.
      !
      ! 1. Differentiate the streamfunction to get u and v velocity.
      !
      !     A. u = dPsi/dy (physical) = i*k2*Psi (spectral)
      !     B. v = -dPsi/dx (physical) = -i*k1*Psi (spectral)
      !
      uC(:,:) = (0.0_RWPC, 0.0_RWPC)
      vC(:,:) = (0.0_RWPC, 0.0_RWPC)
      kxMask(1) = kxMU
      kyMask(1) = kyMU
      CALL DuDy(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kyMask, Q(:,:,2), uC)
      CALL DuDx(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kxMask, Q(:,:,2), vC)
      vC(kxLU:kxMU,kyLU:kyMU) = -1.0_RWPC*vC(kxLU:kxMU,kyLU:kyMU)
      !
      ! 2. Explicitly truncate higher frequency modes before inversion.
      !
      !    NOTE: This is already done for u and v, since Q is explicitly
      !    truncated at the end of each RK stage.
      !
      !    We accomplish vorticity truncation by only copying part of the array.
      !
      wC(:,:) = (0.0_RWPC, 0.0_RWPC)
      wC(kxLU:kxMU,kyLU:kyMU) = Q(kxLU:kxMU,kyLU:kyMU,1)
      !
      ! 3. Invert the vorticity and velocity components to physical space.
      !
      CALL FFTW_MPI_EXECUTE_DFT_C2R(c2rPlan, uC, uR)
      CALL FFTW_MPI_EXECUTE_DFT_C2R(c2rPlan, vC, vR)
      CALL FFTW_MPI_EXECUTE_DFT_C2R(c2rPlan, wC, wR)
      !
      ! 4. Form the nonlinear product terms w*u and w*v in physical space.
      !
      DO j = 1, nyW
         DO i = 1, nxW
            uR(i,j) = uR(i,j)*wR(i,j)
            vR(i,j) = vR(i,j)*wR(i,j)
         END DO
      END DO
      !
      ! 5. Transform the nonlinear terms back to spectral space.
      !
      CALL FFTW_MPI_EXECUTE_DFT_R2C(r2cPlan, uR, uC)
      CALL FFTW_MPI_EXECUTE_DFT_R2C(r2cPlan, vR, vC)
      norm = REAL(nxW*nyW, RWPC)
      uC(:,:) = uC(:,:)/norm
      vC(:,:) = vC(:,:)/norm
      !
      ! 6. Perform differentation of the nonlinear term in spectral space.
      !
      !    NOTE: After this step, the nonlinear terms are ready to go.
      !
      CALL DuDx(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kxMask, uC, uC)
      CALL DuDy(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kyMask, vC, vC)

      ! Form the RHS of the ODE system in spectral space.
      !
      ! NOTE: Here we only loop over the wavenumbers being used, not the total
      ! number of wavenumbers resolved by the grid.
      !
      ! 1. Loop over all wavennumbers for the vorticity and increment dQ.
      !
      DO ky = kyLU, kyMU
         DO kx = kxLU, kxMU
            dQ(kx,ky) = -uC(kx,ky) - vC(kx,ky) - nu*(kx**2 + ky**2)*Q(kx,ky,1)
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
   !> @param[in] kxLT Absolute smallest x wavenumber (truncated).
   !> @param[in] kyLT Absolute smallest y wavenumber (truncated).
   !> @param[in] kxMT Absolute largest x wavenumber (truncated).
   !> @param[in] kyMT Absolute largest y wavenumber (truncated).
   !> @param[in] kxLU Smallest useful x wavenumber.
   !> @param[in] kyLU Smallest useful y wavenumber.
   !> @param[in] kxMU Largest useful x wavenumber.
   !> @param[in] kyMU Largest useful y wavenumber.
   !> @param[in] Qw Current state vector for vorticity.
   !> @param[out] Qp Current state vector for the streamfunction.
   SUBROUTINE ComputePsi(kxLT, kyLT, kxMT, kyMT, &
                         kxLU, kyLU, kxMU, kyMU, &
                         Qw, Qp)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLT, kyLT, kxMT, kyMT
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLU, kyLU, kxMU, kyMU
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(IN) :: Qw
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(OUT) :: Qp
      ! Local variables.
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: kx, ky

      ! Loop over all useful wavenumbers and compute the streamfunction.
      !
      ! NOTE: revisit this in the future...there is divide by zero going on.
      DO ky = kyLU, kyMU
         DO kx = kxLU, kxMU
            Qp(kx,ky) = Qw(kx,ky)/(REAL(kx, RWPC)**2 + REAL(ky, RWPC)**2)
         END DO
      END DO
      !
      ! Explicitly set the zero mode to zero, since this is indeterminant.
      Qp(0,0) = (0.0_RWPC, 0.0_RWPC)
   END SUBROUTINE ComputePsi

   !> Routine to differentiate a signal in the x direction.
   !!
   !! In this routine we multiply a signal by i*kx, where kx is the wavenumber
   !! in the x direction. Wavenumbers in the mask array will have their values
   !! explicitly set to zero after differentiation. These will correspond to all
   !! j indices for each i index in the mask array.
   !!
   !> @param[in] kxL Lowest wavenumber index in the data arrays.
   !> @param[in] kxM Maximum wavenumber index in the data arrays.
   !> @param[in] kyL Lowest wavenumber index in the data arrays.
   !> @param[in] kyM Maximum wavenumber index in the data arrays.
   !> @param[in] kxStart Starting wavenumber in x for differentiation.
   !> @param[in] kxEnd Ending wavenumber in x for differentiation.
   !> @param[in] kyStart Starting wavenumber in y for differentiation.
   !> @param[in] kyEnd Ending wavenumber in y for differentiation.
   !> @param[in] numMask Number of x wavenumbers being masked.
   !> @param[in] kxMask Array of wavenumbers being masked.
   !> @param[in] sigIn Signal being differentiated.
   !> @param[out] sigOut Differentiated signal
   SUBROUTINE DuDx(kxL, kxM, kyL, kyM, kxStart, kxEnd, kyStart, kyEnd, &
                   numMask, kxMask, sigIn, sigOut)
      ! Required modules.
      USE Parameters_m,ONLY: II
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: kxL, kxM, kyL, kyM, numMask
      INTEGER(KIND=IWPF),INTENT(IN) :: kxStart, kxEnd, kyStart, kyEnd
      INTEGER(KIND=IWPF),DIMENSION(numMask),INTENT(IN) :: kxMask
      COMPLEX(KIND=CWPC),DIMENSION(kxL:kxM,kyL:kyM),INTENT(IN) :: sigIn
      COMPLEX(KIND=CWPC),DIMENSION(kxL:kxM,kyL:kyM),INTENT(OUT) :: sigOut
      ! Local variables.
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: kx, ky
      ! Looping index for masked wavenumbers.
      INTEGER(KIND=IWPF) :: i

      ! Loop over memory in contiguous fashion and differentiate.
      DO ky = kyStart, kyEnd
         DO kx = kxStart, kxEnd
            ! Differentiate the signal in x.
            sigOut(kx,ky) = II*REAL(kx, RWPC)*sigIn(kx,ky)
         END DO
      END DO

      ! Zero out undesired modes.
      DO i = 1, numMask
         ! The wavenumber being masked.
         kx = kxMask(i)
         !
         ! Zero out all modes with this kx wavenumber.
         sigOut(kx,kyStart:kyEnd) = (0.0_RWPC, 0.0_RWPC)
      END DO
   END SUBROUTINE DuDx

   !> Routine to differentiate a signal in the y direction.
   !!
   !! In this routine we multiply a signal by i*ky, where ky is the wavenumber
   !! in the y direction. Wavenumbers in the mask array will have their values
   !! explicitly set to zero after differentiation. These will correspond to all
   !! i indices for each j index in the mask array.
   !!
   !> @param[in] kxL Lowest wavenumber index in the data arrays.
   !> @param[in] kxM Maximum wavenumber index in the data arrays.
   !> @param[in] kyL Lowest wavenumber index in the data arrays.
   !> @param[in] kyM Maximum wavenumber index in the data arrays.
   !> @param[in] kxStart Starting wavenumber in x for differentiation.
   !> @param[in] kxEnd Ending wavenumber in x for differentiation.
   !> @param[in] kyStart Starting wavenumber in y for differentiation.
   !> @param[in] kyEnd Ending wavenumber in y for differentiation.
   !> @param[in] numMask Number of y wavenumbers being masked.
   !> @param[in] kyMask Array of wavenumbers being masked.
   !> @param[in] sigIn Signal being differentiated.
   !> @param[out] sigOut Differentiated signal
   SUBROUTINE DuDy(kxL, kxM, kyL, kyM, kxStart, kxEnd, kyStart, kyEnd, &
                   numMask, kyMask, sigIn, sigOut)
      ! Required modules.
      USE Parameters_m,ONLY: II
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: kxL, kxM, kyL, kyM, numMask
      INTEGER(KIND=IWPF),INTENT(IN) :: kxStart, kxEnd, kyStart, kyEnd
      INTEGER(KIND=IWPF),DIMENSION(numMask),INTENT(IN) :: kyMask
      COMPLEX(KIND=CWPC),DIMENSION(kxL:kxM,kyL:kyM),INTENT(IN) :: sigIn
      COMPLEX(KIND=CWPC),DIMENSION(kxL:kxM,kyL:kyM),INTENT(OUT) :: sigOut
      ! Local variables.
      ! Looping indices for kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: kx, ky
      ! Looping index for masked wavenumbers.
      INTEGER(KIND=IWPF) :: i

      ! Loop over memory in contiguous fashion and differentiate.
      DO ky = kyStart, kyEnd
         DO kx = kxStart, kxEnd
            ! Differentiate the signal in x.
            sigOut(kx,ky) = II*REAL(ky, RWPC)*sigIn(kx,ky)
         END DO
      END DO

      ! Zero out undesired modes.
      DO i = 1, numMask
         ! The wavenumber being masked.
         ky = kyMask(i)
         !
         ! Zero out all modes with this kx wavenumber.
         sigOut(kxStart:kxEnd,ky) = (0.0_RWPC, 0.0_RWPC)
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

