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
MODULE SetIC_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPF, IWPC, RWPF, RWPC, CWPC

   IMPLICIT NONE

   ! Available initial conditions.
   !
   !> Cosine shear layer moving in the x direction.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: COSINE_SHEAR_X = 1_IWPF

   ! Module procedures.
   PUBLIC :: SetICCosineShearX

CONTAINS

   !> Procedure to set the initial condition to u(x,y) = cos(y).
   !!
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] nxW Number of grid points in the x direction.
   !> @param[in] nyW Number of grid points in the y direction.
   !> @param[in] kxLT Absolute smallest x wavenumber (truncated).
   !> @param[in] kyLT Absolute smallest y wavenumber (truncated).
   !> @param[in] kxMT Absolute largest x wavenumber (truncated).
   !> @param[in] kyMT Absolute largest y wavenumber (truncated).
   !> @param[in] kxLU Smallest useful x wavenumber.
   !> @param[in] kyLU Smallest useful y wavenumber.
   !> @param[in] kxMU Largest useful x wavenumber.
   !> @param[in] kyMU Largest useful y wavenumber.
   !> @param[in] nVar Number of variables in Q.
   !> @param[in] nStg Number of stages in Q for time stepping.
   !> @param[in] uC Complex cast of u velocity array.
   !> @param[in] uR Real cast of u velocity array.
   !> @param[in] vC Complex cast of v velocity array.
   !> @param[in] vR Real cast of v velocity array.
   !> @param[in] Qc Complex cast of Q array.
   !> @param[in] Qr Real cast of Q array.
   SUBROUTINE SetICCosineShearX(rISize, cISize, nxW, nyW, &
                                kxLT, kyLT, kxMT, kyMT, &
                                kxLU, kyLU, kxMU, kyMU, &
                                nVar, nStg, &
                                uC, uR, vC, vR, Qc, Qr)
      ! Required modules.
      USE Parameters_m,ONLY: PI
      USE Spectral_m,ONLY: ComputeRHS, ComputePsi, TransformR2C, DuDy
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize, nxW, nyW
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLT, kyLT, kxMT, kyMT
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLU, kyLU, kxMU, kyMU
      INTEGER(KIND=IWPF),INTENT(IN) :: nVar, nStg
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(INOUT) :: uR, vR
      COMPLEX(CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT,nVar,nStg),INTENT(INOUT) :: Qc
      REAL(KIND=RWPC),DIMENSION(rISize,nyW,nVar,nStg),INTENT(INOUT) :: Qr
      ! Local variables.
      ! Looping indices for physical space.
      INTEGER(KIND=IWPF) :: i, j
      ! Looping indices for spectral space.
      INTEGER(KIND=IWPF) :: kx, ky
      ! Masked wavenumber for differentiation.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kyMask

      ! Loop over the physical grid points and fill in the velocity profile.
      DO j = 1, nyW
         DO i = 1, nxW
            uR(i,j) = COS(REAL(j, RWPC)*2.0_RWPC*REAL(PI, RWPC)/REAL(nyW, RWPC))
         END DO
      END DO
      !
      ! Transform the signal to spectral space.
      CALL TransformR2C(nxW, nyW, rISize, cISize, uR, uC)
      !
      ! Truncate the signal to the useful wavenubers.
      !
      ! Differentiate the signal to form the vorticity.
      kyMask(1) = kyMU
      CALL DuDy(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kyMask, uC, Qc(:,:,1,1))
      !
      ! Solve for the streamfunction based on the vorticity.
      CALL ComputePsi(kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,1,1), Qc(:,:,2,1))
   END SUBROUTINE SetICCosineShearX

END MODULE SetIC_m

