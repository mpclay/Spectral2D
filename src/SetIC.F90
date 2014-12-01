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
MODULE SetIC_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPF, IWPC, RWPF, RWPC, CWPC

   IMPLICIT NONE

   ! Available initial conditions.
   !
   !> Taylor-Green vortex initialization.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: TAYLOR_GREEN_VORTEX = 1_IWPF

   ! Module procedures.
   PUBLIC :: SetTaylorGreen

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

END MODULE SetIC_m

