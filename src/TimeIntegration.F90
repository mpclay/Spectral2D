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
!> @file TimeIntegration.F90
!> @author Matthew Clay
!> @brief Time integration schemes.
MODULE TimeIntegration_m

   ! Required modules.
   USE ISO_C_BINDING,ONLY: C_PTR
   USE Parameters_m,ONLY: IWPC, IWPF, RWPC, RWPF, CWPC

   IMPLICIT NONE

   ! Available time integrators.
   !
   !> TVD RK3 of Shu.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: RK3_TVD_SHU = 1_IWPF

   ! Local variables.
   !
   !> Which time integrator is being used.
   INTEGER(KIND=IWPF),PRIVATE :: scheme

   ! Scratch arrays for time integration.
   !
   !> Increment of Q for the useful wavenumbers in the simulation.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: dQ

   ! Module procedures.
   PUBLIC :: TimeIntegrationSetup, IntegrateOneStep, ComputeTimeStep
   PRIVATE :: RK3TVDTimeUpdate

CONTAINS

   !> Subroutine to set up the time integration module.
   !!
   !> @param[in] scheme_ Which time integration scheme to use.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] nyG Global number of grid points in the y direction.
   !> @param[out] nStg Number of stages needed for time integration scheme.
   !> @param[out] cS Location in Q for this time stage.
   !> @param[out] nS Location in Q for next time stage.
   SUBROUTINE TimeIntegrationSetup(scheme_, rISize, cISize, nyG, nStg, cS, nS)
      ! Required modules.
      USE ISO_C_BINDING,ONLY: C_F_POINTER
      USE Alloc_m,ONLY: Alloc
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: scheme_, rISize, cISize, nyG
      INTEGER(KIND=IWPF),INTENT(OUT) :: nStg, cS, nS

      ! Set which scheme we are going to use.
      SELECT CASE (scheme_)
         CASE (RK3_TVD_SHU)
            ! This is a low-storage scheme. Only two storage locations needed.
            scheme = RK3_TVD_SHU
            nStg = 2_IWPF
            cS = 1_IWPF
            nS = 2_IWPF
         CASE DEFAULT
            ! Add warning message.
      END SELECT

      ! Allocate memory for scratch arrays.
      CALL Alloc(cISize, nyG, 1_IWPF, 1_IWPF, dQ)
   END SUBROUTINE TimeIntegrationSetup

   !> Subroutine to advance one step in time.
   !!
   !> @param[in] nxG Number of grid points in the x direction.
   !> @param[in] nyG Number of grid points in the y direction.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] kxG Array of wavenumbers for the x direction.
   !> @param[in] kyG Array of wavenumbers for the y direction.
   !> @param[in] nu Physical viscosity.
   !> @param[in] dt Time step.
   !> @param[in] nVar Number of variables in Q.
   !> @param[in] nStg Number of stages in Q for time stepping.
   !> @param[in,out] cS Current stage. Stage in Q that is for current step.
   !> @param[in,out] nS Next stage. Stage in Q to hold updated step.
   !> @param[in] uC Complex cast of u velocity array.
   !> @param[in] uR Real cast of u velocity array.
   !> @param[in] vC Complex cast of v velocity array.
   !> @param[in] vR Real cast of v velocity array.
   !> @param[in] Qc Complex cast of Q array.
   !> @param[in] Qr Real cast of Q array.
   SUBROUTINE IntegrateOneStep(nxG, nyG, rISize, cISize, kxG, kyG, &
                               nu, dt, nVar, nStg, cS, nS, &
                               uC, uR, vC, vR, Qc, Qr)
      ! Required modules.
      USE Spectral_m,ONLY: ComputeRHS, ComputePsi
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyG),INTENT(IN) :: kyG
      REAL(KIND=RWPC),INTENT(IN) :: nu, dt
      INTEGER(KIND=IWPF),INTENT(IN) :: nVar, nStg
      INTEGER(KIND=IWPF),INTENT(INOUT) :: cS, nS
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyG),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyG),INTENT(INOUT) :: uR, vR
      COMPLEX(CWPC),DIMENSION(cISize,nyG,nVar,nStg),INTENT(INOUT) :: Qc
      REAL(KIND=RWPC),DIMENSION(rISize,nyG,nVar,nStg),INTENT(INOUT) :: Qr
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j

      ! Zero out everything in the next stage holder and dQ.
      Qc(:,:,:,nS) = (0.0_RWPC, 0.0_RWPC)
      dQ(:,:) = (0.0_RWPC, 0.0_RWPC)

      ! Calculate the RHS with the data at the start of the time step.
      CALL ComputeRHS(nxG, nyG, rISize, cISize, &
                      kxG, kyG, nu, nVar, &
                      uC, uR, vC, vR, Qc(:,:,:,cS), Qr(:,:,:,cS), dQ)
      !
      ! Determine omega at the first intermediate RK stage.
      DO j = 1, nyG
         DO i = 1, cISize
            Qc(i,j,1,nS) = Qc(i,j,1,cS) + dt*dQ(i,j)
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(cISize, nyG, kxG, kyG, Qc(:,:,1,nS), Qc(:,:,2,nS))

      ! Calculate the RHS using the first intermediate stage.
      dQ(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL ComputeRHS(nxG, nyG, rISize, cISize, &
                      kxG, kyG, nu, nVar, &
                      uC, uR, vC, vR, Qc(:,:,:,nS), Qr(:,:,:,nS), dQ)
      !
      ! Determine the vorticity at the second intermediate RK stage.
      DO j = 1, nyG
         DO i = 1, cISize
            Qc(i,j,1,nS) = 0.75_RWPC*Qc(i,j,1,cS) + &
                           0.25_RWPC*Qc(i,j,1,nS) + &
                           0.25_RWPC*dt*dq(i,j)
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(cISize, nyG, kxG, kyG, Qc(:,:,1,nS), Qc(:,:,2,nS))

      ! Calculate the RHS using the second intermediate stage.
      dQ(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL ComputeRHS(nxG, nyG, rISize, cISize, &
                      kxG, kyG, nu, nVar, &
                      uC, uR, vC, vR, Qc(:,:,:,nS), Qr(:,:,:,nS), dQ)
      !
      ! Determine the vorticity at the next time step.
      DO j = 1, nyG
         DO i = 1, cISize
            Qc(i,j,1,cS) = Qc(i,j,1,cS)/3.0_RWPC + &
                           2.0_RWPC*Qc(i,j,1,nS)/3.0_RWPC + &
                           dt*2.0_RWPC*dQ(i,j)/3.0_RWPC
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(cISize, nyG, kxG, kyG, Qc(:,:,1,cS), Qc(:,:,2,cS))

      ! For this low-storage scheme we do not change cS from 1.
      cS = 1_IWPF
      nS = 2_IWPF
   END SUBROUTINE IntegrateOneStep

   !> Subroutine to calculate the time step.
   SUBROUTINE ComputeTimeStep()
      IMPLICIT NONE
   END SUBROUTINE ComputeTimeStep

   !> Update one step of the RK3 TVD scheme of Shu.
   SUBROUTINE RK3TVDTimeUpdate()
      IMPLICIT NONE
   END SUBROUTINE RK3TVDTimeUpdate

END MODULE TimeIntegration_m

