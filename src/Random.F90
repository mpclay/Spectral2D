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
!> @file Random.F90
!> @author Matthew Clay
!> @brief Module to handle random numbers, etc. for the code.
MODULE Random_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPF

   IMPLICIT NONE

   !> Size of the seed (compiler dependent).
   INTEGER(KIND=IWPF),PRIVATE :: seedSize
   !> Random seed.
   INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE,PRIVATE :: seed

CONTAINS

   !> Procedure to initialize the random seed.
   !!
   !> @param[in] rank Rank of this process.
   SUBROUTINE InitRandomSeed(rank)
      ! Required modules.
      USE MPI
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rank
      ! Local variables.
      ! IO unit for reading random data to fill the seed.
      INTEGER(KIND=IWPF) :: un
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! The size of the seed is compiler dependent.
      CALL RANDOM_SEED(SIZE=seedSize)

      ! Allocate memory for the seed.
      IF (ALLOCATED(seed)) THEN
         DEALLOCATE(seed)
      END IF
      ALLOCATE(seed(seedSize), STAT=ierr)

      ! The root process reads random data to fill the seed and then broadcasts
      ! it to the other processes.
      IF (rank == 0) THEN
         OPEN(NEWUNIT=un, FILE="/dev/urandom", ACCESS="STREAM", &
              FORM="UNFORMATTED", ACTION="READ", STATUS="OLD")
         READ(UNIT=un) seed
         CLOSE(UNIT=un)
      END IF

      ! Make sure all processes have the same random seed.
      CALL MPI_BCAST(seed, seedSize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      CALL RANDOM_SEED(PUT=seed)
   END SUBROUTINE InitRandomSeed

END MODULE Random_m

