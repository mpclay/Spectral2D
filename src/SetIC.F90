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
   !> @param[in] wC Complex cast of vorticity array.
   !> @param[in] wR Real cast of vorticity array.
   !> @param[in] psiC Complex cast of streamfunction array.
   !> @param[in] psiR Real cast of streamfunction array.
   SUBROUTINE SetICCosineShearX(nxG, nyG, rISize, cISize, kxG, kyG, &
                                uC, uR, vC, vR, wC, wR, psiC, psiR)
      ! Required modules.
      USE Parameters_m,ONLY: PI
      USE Spectral_m,ONLY: ComputeVorticity, ComputePsi, TransformR2C
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(INOUT) :: uC, vC, wC, psiC
      REAL(KIND=RWPC),DIMENSION(rISize,nyG),INTENT(INOUT) :: uR, vR, wR, psiR
      ! Local variables.
      ! Looping indices for physical space.
      INTEGER(KIND=IWPF) :: i, j

      ! Loop over the physical grid points and fill in the velocity profile.
      DO j = 1, nyG
         DO i = 1, nxG
            uR(i,j) = COS(REAL(j, RWPC)*2.0_RWPC*REAL(PI, RWPC)/REAL(nyG, RWPC))
         END DO
      END DO
      !
      ! Transform the signal to spectral space.
      CALL TransformR2C(nxG, nyG, rISize, cISize, uR, uC)
      CALL TransformR2C(nxG, nyG, rISize, cISize, vR, vC)
      !
      ! Truncate the signal to the useful wavenubers.
      !
      ! Differentiate the signal to form the vorticity.
      CALL ComputeVorticity(nxG, nyG, cISize, kxG, kyG, uC, vC, wC)
      !
      ! Solve for the streamfunction based on the vorticity.
      CALL ComputePsi(cISize, nyG, kxG, kyG, wC, psiC)
   END SUBROUTINE SetICCosineShearX

END MODULE SetIC_m

