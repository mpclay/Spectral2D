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
   !> Scratch array for the vorticity when computing the nonlinear term.
   TYPE(C_PTR),PRIVATE :: w
   !> Real cast of w for working in physical space.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: wR
   !> Complex cast of w for working in spectral space.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: wC

   ! Module procedures.
   PUBLIC :: TimeIntegrationSetup, IntegrateOneStep, ComputeTimeStep
   PRIVATE :: RK3TVDTimeUpdate

CONTAINS

   !> Subroutine to set up the time integration module.
   !!
   !> @param[in] scheme_ Which time integration scheme to use.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] nyW Working number of grid points in the y direction.
   !> @param[in] kxLU Smallest useful x wavenumber.
   !> @param[in] kyLU Smallest useful y wavenumber.
   !> @param[in] kxMU Largest useful x wavenumber.
   !> @param[in] kyMU Largest useful y wavenumber.
   !> @param[out] nStg Number of stages needed for time integration scheme.
   !> @param[out] cS Location in Q for this time stage.
   !> @param[out] nS Location in Q for next time stage.
   SUBROUTINE TimeIntegrationSetup(scheme_, rISize, cISize, nyW, &
                                   kxLU, kyLU, kxMU, kyMU, nStg, cS, nS)
      ! Required modules.
      USE ISO_C_BINDING,ONLY: C_F_POINTER
      USE Alloc_m,ONLY: Alloc
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: scheme_, rISize, cISize, nyw
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLU, kyLU, kxMU, kyMU
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
      CALL Alloc(cISize, nyW, kxLU, kyLU, dQ)
      CALL Alloc(cISize, nyW, w)
      CALL C_F_POINTER(w, wR, [rISize,nyW])
      CALL C_F_POINTER(w, wC, [cISize,nyW])
      WRITE(*,*) 'SHAPE(dQ) = ', SHAPE(dQ)
      WRITE(*,*) 'SHAPE(wR) = ', SHAPE(wR)
      WRITE(*,*) 'SHAPE(wC) = ', SHAPE(wC)
   END SUBROUTINE TimeIntegrationSetup

   !> Subroutine to advance one step in time.
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
   !> @param[in] nu Physical viscosity.
   !> @param[in] dt Time step.
   !> @param[in] nVar Number of variables in Q.
   !> @param[in] nStg Number of stages in Q for time stepping.
   !> @param[out] cStg Current stage. Stage in Q that is for current step.
   !> @param[out] nStg Next stage. Stage in Q to hold updated step.
   !> @param[in] uC Complex cast of u velocity array.
   !> @param[in] uR Real cast of u velocity array.
   !> @param[in] vC Complex cast of v velocity array.
   !> @param[in] vR Real cast of v velocity array.
   !> @param[in] Qc Complex cast of Q array.
   !> @param[in] Qr Real cast of Q array.
   SUBROUTINE IntegrateOneStep(rISize, cISize, nxW, nyW, &
                               kxLT, kyLT, kxMT, kyMT, &
                               kxLU, kyLU, kxMU, kyMU, &
                               nu, dt, &
                               nVar, nStg, cS, nS, &
                               uC, uR, vC, vR, Qc, Qr)
      ! Required modules.
      USE Spectral_m,ONLY: ComputeRHS, ComputePsi
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize, nxW, nyW
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLT, kyLT, kxMT, kyMT
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLU, kyLU, kxMU, kyMU
      REAL(KIND=RWPC),INTENT(IN) :: nu, dt
      INTEGER(KIND=IWPF),INTENT(IN) :: nVar, nStg
      INTEGER(KIND=IWPF),INTENT(INOUT) :: cS, nS
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(INOUT) :: uR, vR
      COMPLEX(CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT,nVar,nStg),INTENT(INOUT) :: Qc
      REAL(KIND=RWPC),DIMENSION(rISize,nyW,nVar,nStg),INTENT(INOUT) :: Qr
      ! Local variables.
      ! Looping indices for wavenumbers.
      INTEGER(KIND=IWPF) :: kx, ky

      ! Zero out everything in the next stage holder and dQ.
      Qc(:,:,:,nS) = (0.0_RWPC, 0.0_RWPC)
      dQ(:,:) = (0.0_RWPC, 0.0_RWPC)

      ! Calculate the RHS with the data at the start of the time step.
      CALL ComputeRHS(rISize, nxW, nyW, &
                      kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,:,cS), nu, &
                      uC, vC, wC, uR, vR, wR, dQ)
      !
      ! Determine omega at the first intermediate RK stage.
      DO ky = kyLU, kyMU
         DO kx = kxLU, kxMU
            Qc(kx,ky,1,nS) = Qc(kx,ky,1,cS) + dt*dQ(kx,ky)
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,1,nS), Qc(:,:,2,nS))

      ! Calculate the RHS using the first intermediate stage.
      dQ(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL ComputeRHS(rISize, nxW, nyW, &
                      kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,:,nS), nu, &
                      uC, vC, wC, uR, vR, wR, dQ)
      !
      ! Determine the vorticity at the second intermediate RK stage.
      DO ky = kyLU, kyMU
         DO kx = kxLU, kxMU
            Qc(kx,ky,1,nS) = 0.75_RWPC*Qc(kx,ky,1,cS) + &
                             0.25_RWPC*Qc(kx,ky,1,nS) + &
                             0.25_RWPC*dt*dq(kx,ky)
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,1,nS), Qc(:,:,2,nS))

      ! Calculate the RHS using the second intermediate stage.
      dQ(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL ComputeRHS(rISize, nxW, nyW, &
                      kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,:,nS), nu, &
                      uC, vC, wC, uR, vR, wR, dQ)
      !
      ! Determine the vorticity at the next time step.
      DO ky = kyLU, kyMU
         DO kx = kxLU, kxMU
            Qc(kx,ky,1,cS) = Qc(kx,ky,1,cS)/3.0_RWPC + &
                             2.0_RWPC*Qc(kx,ky,1,nS)/3.0_RWPC + &
                             dt*2.0_RWPC*dQ(kx,ky)/3.0_RWPC
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(kxLT, kyLT, kxMT, kyMT, &
                      kxLU, kyLU, kxMU, kyMU, &
                      Qc(:,:,1,cS), Qc(:,:,2,cS))

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

