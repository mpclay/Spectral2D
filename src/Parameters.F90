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
!> @file Parameters.F90
!> @author Matthew Clay
!> @brief Defines working precision and universal constants.
MODULE Parameters_m

   ! Required modules.
   USE ISO_FORTRAN_ENV,ONLY: REAL64, INT32
   USE ISO_C_BINDING,ONLY: C_DOUBLE, C_DOUBLE_COMPLEX, C_INTPTR_T

   IMPLICIT NONE

   ! Variable precision.
   !
   !> Working precision for reals with C.
   INTEGER,PARAMETER,PUBLIC :: RWPC = C_DOUBLE
   !> Working precision for reals with Fortran.
   INTEGER,PARAMETER,PUBLIC :: RWPF = REAL64
   !> Working precision for complex with C.
   INTEGER,PARAMETER,PUBLIC :: CWPC = C_DOUBLE_COMPLEX
   !> Working precision for integers with C.
   INTEGER,PARAMETER,PUBLIC :: IWPC = C_INTPTR_T
   !> Working precision for integer with Fortran.
   INTEGER,PARAMETER,PUBLIC :: IWPF = INT32

   ! Fundamental constants.
   !
   !> Pi.
   REAL(KIND=RWPF),PARAMETER,PUBLIC :: PI = ACOS(-1.0_RWPF)

END MODULE Parameters_m

