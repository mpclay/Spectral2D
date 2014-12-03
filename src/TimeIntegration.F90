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
!> @brief Time integration module.
!!
!! Right now this module only implements the TVD RK3 scheme of Shu. This is a
!! low-storage scheme, which is able to work with the data array allocation in
!! the main code. By that we mean that it is not straightforward to use a Q
!! array for multiple variables and time stages with the allocation procedure
!! of MPI FFTW. Because the MPI FFTW routines require larger arrays than the
!! actual (cISize X nyP) number of complex variables, and we cannot cast that
!! array in a manner that can be easily used in Fortran, we are stuck with
!! individual arrays for each variable. In the future, we can probably achieve
!! the desired effect by casting arrays like (1:allocLocal,1:nVar,1:nStg) by
!! sacrificing the desired (1:cISize,1:nyP,1:nVar,1:nStg) shape.
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
   !> Increment of vorticity.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),ALLOCATABLE,PRIVATE :: dW
   !> Scratch data array for the vorticity.
   TYPE(C_PTR) :: wS
   !> Real cast of the scratch vorticity array.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: wSR
   !> Complex cast of the scratch vorticity array.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: wSC
   !> Scratch data array for the streamfunction.
   TYPE(C_PTR) :: psiS
   !> Real cast of the scratch streamfunction array.
   REAL(KIND=RWPC),DIMENSION(:,:),POINTER :: psiSR
   !> Complex cast of the scratch streamfunction array.
   COMPLEX(KIND=CWPC),DIMENSION(:,:),POINTER :: psiSC

   ! Module procedures.
   PUBLIC :: TimeIntegrationSetup, TimeIntegrationFinalize
   PUBLIC :: IntegrateOneStep, ComputeTimeStep

CONTAINS

   !> Subroutine to set up the time integration module.
   !!
   !> @param[in] scheme_ Which time integration scheme to use.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] allocLocal Size of data arrays required by FFTW MPI.
   SUBROUTINE TimeIntegrationSetup(scheme_, rISize, cISize, nyP, allocLocal)
      ! Required modules.
      USE ISO_C_BINDING,ONLY: C_F_POINTER
      USE Alloc_m,ONLY: Alloc
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: scheme_, rISize, cISize, nyP, allocLocal

      ! Set which scheme we are going to use.
      SELECT CASE (scheme_)
         CASE (RK3_TVD_SHU)
            ! This is a low-storage scheme. Only two storage locations needed.
            scheme = RK3_TVD_SHU
         CASE DEFAULT
            ! Add warning message.
      END SELECT

      ! Allocate memory for scratch arrays.
      CALL Alloc(cISize, nyP, 1_IWPF, 1_IWPF, dW)
      CALL Alloc(allocLocal, wS)
      CALL Alloc(allocLocal, psiS)
      !
      ! Cast useable portions of the arrays for fortran use.
      CALL C_F_POINTER(wS, wSR, [rISize,nyP])
      CALL C_F_POINTER(wS, wSC, [cISize,nyP])
      wSC(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL C_F_POINTER(psiS, psiSR, [rISize,nyP])
      CALL C_F_POINTER(psiS, psiSC, [cISize,nyP])
      psiSC(:,:) = (0.0_RWPC, 0.0_RWPC)
   END SUBROUTINE TimeIntegrationSetup

   !> Subroutine to advance one step in time.
   !!
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of x grid points for this process.
   !> @param[in] nyP Number of y grid points for this process.
   !> @param[in] rISize Size of i dimension for real numbers.
   !> @param[in] cISize Size of i dimension for complex numbers.
   !> @param[in] kxP Array of wavenumbers for the x direction for this process.
   !> @param[in] kyP Array of wavenumbers for the y direction for this process.
   !> @param[in] nu Physical viscosity.
   !> @param[in] dt Time step.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] uR Real cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[in,out] vR Real cast of v velocity array.
   !> @param[in,out] wC Complex cast of vorticity array.
   !> @param[in,out] wR Real cast of vorticity array.
   !> @param[in,out] psiC Complex cast of streamfunction array.
   !> @param[in,out] psiR Real cast of streamfunction array.
   SUBROUTINE IntegrateOneStep(nxG, nyG, nxP, nyP, rISize, cISize, kxP, kyP, &
                               nu, dt, uC, uR, vC, vR, wC, wR, psiC, psiR)
      ! Required modules.
      USE Spectral_m,ONLY: ComputeRHS, ComputePsi, Truncate
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nxP, nyP, rISize, cISize
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      REAL(KIND=RWPC),INTENT(IN) :: nu, dt
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC, wC, psiC
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: uR, vR, wR, psiR
      ! Local variables.
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j

      ! Zero out everything in the next stage holder and dQ.
      wSC(:,:) = (0.0_RWPC, 0.0_RWPC)
      psiSC(:,:) = (0.0_RWPC, 0.0_RWPC)
      dW(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Calculate the RHS with the data at the start of the time step.
      CALL ComputeRHS(nxG, nyG, nxP, nyP, rISize, cISize, kxP, kyP, nu, &
                      uC, uR, vC, vR, wC, wR, psiC, psiR, dW)
      !
      ! Determine omega at the first intermediate RK stage.
      DO j = 1, nyP
         DO i = 1, cISize
            wSC(i,j) = wC(i,j) + dt*dW(i,j)
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(cISize, nyP, kxP, kyP, wSC, psiSC)
      !
      ! Truncate undesired high-wavenumber modes.
      CALL Truncate(cISize, nyP, wSC)
      CALL Truncate(cISize, nyP, psiSC)

      ! Calculate the RHS using the first intermediate stage.
      dW(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL ComputeRHS(nxG, nyG, nxP, nyP, rISize, cISize, kxP, kyP, nu, &
                      uC, uR, vC, vR, wSC, wSR, psiSC, psiSR, dW)
      !
      ! Determine the vorticity at the second intermediate RK stage.
      DO j = 1, nyP
         DO i = 1, cISize
            wSC(i,j) = 0.75_RWPC*wC(i,j) + &
                       0.25_RWPC*wSC(i,j) + &
                       0.25_RWPC*dt*dW(i,j)
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(cISize, nyP, kxP, kyP, wSC, psiSC)
      !
      ! Truncate undesired high-wavenumber modes.
      CALL Truncate(cISize, nyP, wSC)
      CALL Truncate(cISize, nyP, psiSC)

      ! Calculate the RHS using the second intermediate stage.
      dW(:,:) = (0.0_RWPC, 0.0_RWPC)
      CALL ComputeRHS(nxG, nyG, nxP, nyP, rISize, cISize, kxP, kyP, nu, &
                      uC, uR, vC, vR, wSC, wSR, psiSC, psiSR, dW)
      !
      ! Determine the vorticity at the next time step.
      DO j = 1, nyP
         DO i = 1, cISize
            wC(i,j) = wC(i,j)/3.0_RWPC + &
                      2.0_RWPC*wSC(i,j)/3.0_RWPC + &
                      dt*2.0_RWPC*dW(i,j)/3.0_RWPC
         END DO
      END DO
      !
      ! Now with vorticity updated, solve for the streamfunction.
      CALL ComputePsi(cISize, nyP, kxP, kyP, wC, psiC)
      !
      ! Truncate undesired high-wavenumber modes.
      CALL Truncate(cISize, nyP, wC)
      CALL Truncate(cISize, nyP, psiC)
   END SUBROUTINE IntegrateOneStep

   !> Subroutine to calculate the time step.
   SUBROUTINE ComputeTimeStep()
      IMPLICIT NONE
   END SUBROUTINE ComputeTimeStep

   !> Routine to finalize the time integration module.
   SUBROUTINE TimeIntegrationFinalize()
      ! Required modules.
      USE ISO_C_BINDING
      USE Alloc_m,ONLY: Dealloc
      IMPLICIT NONE
      INCLUDE 'fftw3-mpi.f03'

      ! Free all working arrays.
      DEALLOCATE(dW)
      CALL Dealloc(wS)
      CALL Dealloc(psiS)
      wSR => NULL()
      wSC => NULL()
      psiSR => NULL()
      psiSC => NULL()
   END SUBROUTINE TimeIntegrationFinalize

END MODULE TimeIntegration_m

