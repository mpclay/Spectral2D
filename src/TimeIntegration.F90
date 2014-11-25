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
   USE Parameters_m,ONLY: IWPC, IWPF, RWPC, RWPF, CWPC

   IMPLICIT NONE

   ! Available time integrators.
   !
   !> TVD RK3 of Shu.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: RK3_TVD_SHU = 1_IWPF

   ! Module procedures.
   PUBLIC :: TimeIntegrationSetup, IntegrateOneStep, ComputeTimeStep
   PRIVATE :: RK3TVDTimeUpdate

CONTAINS

   !> Subroutine to set up the time integration module.
   SUBROUTINE TimeIntegrationSetup()
      IMPLICIT NONE
   END SUBROUTINE TimeIntegrationSetup

   !> Subroutine to advance one step in time.
   SUBROUTINE IntegrateOneStep()
      IMPLICIT NONE
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

