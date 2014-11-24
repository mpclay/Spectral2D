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
!> @file Messages.F90
!> @author Matthew Clay
!> @brief Error and warning messages.
MODULE Messages_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPF

   IMPLICIT NONE

   ! Module procedures.
   PUBLIC :: Error, ErrorCheck

CONTAINS

   !> Error handling and messages.
   !!
   !! Use this to pass in an error message to be printed to the user. The code
   !! can be halted if the user desires.
   !!
   !> @param[in] comm MPI communicator.
   !> @param[in] rank MPI rank.
   !> @param[in] msg Message to be printed.
   !> @param[in,optional] halt Whether or not to halt the code.
   SUBROUTINE Error(comm, rank, msg, halt)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: OUTPUT_UNIT
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: comm, rank
      CHARACTER(LEN=*),INTENT(IN) :: msg
      LOGICAL,OPTIONAL,INTENT(IN) :: halt
      ! Local variables.
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr
      ! Logic to halt the code.
      LOGICAL :: halt_

      ! Handle optional arguments.
      IF (PRESENT(halt)) THEN
         halt_ = halt
      ELSE
         halt_ = .FALSE.
      END IF

      ! Print out the message to the user.
      IF (rank == 0) THEN
         WRITE(OUTPUT_UNIT,10,ADVANCE="NO") "SPECTRAL ERROR: ", TRIM(msg), "."
         IF (halt_) THEN
            WRITE(OUTPUT_UNIT,20) " HALTING."
         ELSE
            WRITE(OUTPUT_UNIT,20) " Continuing execution."
         END IF
      END IF
      10 FORMAT (A,A,A)
      20 FORMAT (A)

      ! If necessary, halt execution.
      IF (halt_) THEN
         CALL MPI_FINALIZE(comm, ierr)
         CALL EXIT(1)
      END IF
   END SUBROUTINE Error

   !> Check for errors.
   !!
   !! In the following we refer to 'the call' and 'the function' as the function
   !! the user called somewhere else that returned an error token.
   !!
   !> @param[in] comm MPI communicator.
   !> @param[in] rank MPI rank.
   !> @param[in] str Some information about the call being made.
   !> @param[in] ierr The error value returned from the function call.
   !> @param[in,optional] error Expected value from the function when no error
   !! has occured. Default value is 0.
   SUBROUTINE ErrorCheck(comm, rank, str, ierr, error)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: OUTPUT_UNIT
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: str
      INTEGER(KIND=IWPF),INTENT(IN) :: comm, rank, ierr
      INTEGER(KIND=IWPF),OPTIONAL,INTENT(IN) :: error
      ! Local variables
      ! Error handling passed to this function.
      INTEGER(KIND=IWPF) :: error_
      ! Local error handling.
      INTEGER(KIND=IWPF) :: err

      ! Handling optional arguments.
      IF (PRESENT(error)) THEN
         error_ = error
      ELSE
         error_ = 0_IWPF
      END IF

      ! An error occured if ierr is not equal to error.
      IF (ierr /= error_) THEN
         IF (rank == 0_IWPF) THEN
            WRITE(OUTPUT_UNIT,10) "SPECTRAL ERROR: an error occured for ", &
                                  TRIM(str), ". HALTING."
         END IF
         10 FORMAT (A,A,A)
         CALL MPI_FINALIZE(comm, err)
         CALL EXIT(1)
      END IF
   END SUBROUTINE ErrorCheck

END MODULE Messages_m

