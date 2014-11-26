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
!> @file IO.F90
!> @author Matthew Clay
!> @brief File to handle IO for the program.
MODULE IO_m

   ! Required modules.
   USE ISO_FORTRAN_ENV
   USE HDF5
   USE Parameters_m,ONLY: IWPF, IWPC, RWPF, RWPC, CWPC

   IMPLICIT NONE

   !> File name length for output files.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: FILE_NAME_LENGTH = 256_IWPF

   !> Interface for the file name function.
   INTERFACE FileName
      MODULE PROCEDURE FileNameWithoutNum
      MODULE PROCEDURE FileNameWithNum
   END INTERFACE FileName

   ! Module procedures.
   PUBLIC :: WriteGridPlot3d, WriteRestartPlot3d
   PRIVATE :: FileNameWithoutNum, FileNameWithNum

CONTAINS

   !> Routine to form a file name with a root name and file suffix.
   !!
   !> @param[in] root Root name for the file.
   !> @param[in] suffix Suffix for the file.
   !> @param[out] fname File name: root.suffix
   SUBROUTINE FileNameWithoutNum(root, suffix, fname)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: root, suffix
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fname

      ! Form the output file name.
      WRITE(fname,10) TRIM(root), '.', TRIM(suffix)
      10 FORMAT (A,A,A)
   END SUBROUTINE FileNameWithoutNum

   !> Routine to form a file name with a root name, number, and file suffix.
   !!
   !> @param[in] root Root name for the file.
   !> @param[in] num Number of the file.
   !> @param[in] suffix Suffix for the file.
   !> @param[out] fname File name: root_num.suffix
   SUBROUTINE FileNameWithNum(root, num, suffix, fname)
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: root, suffix
      INTEGER(KIND=IWPF),INTENT(IN) :: num
      CHARACTER(LEN=FILE_NAME_LENGTH),INTENT(OUT) :: fname

      ! Form the output file name.
      WRITE(fname,10) TRIM(root), '_', num, '.', TRIM(suffix)
      10 FORMAT (A,A,I8.8,A,A)
   END SUBROUTINE FileNameWithNum

   !> Routine to write the grid in Plot3d format.
   !!
   !> @param[in] nxW Number of grid points in the x direction.
   !> @param[in] nyW Number of grid points in the y direction.
   SUBROUTINE WriteGridPlot3d(nxW, nyW)
      ! Required modules.
      USE Parameters_m,ONLY: PI
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxW, nyW
      ! Local variables.
      ! Looping indices for x and y directions.
      INTEGER(KIND=IWPF) :: i, j

      ! Open the file for reading.
      OPEN(UNIT=100,FILE='grid.xyz',STATUS='REPLACE',FORM='FORMATTED')
      !
      ! Write out the number of domains.
      WRITE(100,150) 1_IWPF
      150 FORMAT (I6.6)
      !
      ! Write out the grid dimensions.
      WRITE(100,200) nxW, nyW
      200 FORMAT (I6.6,1X,I6.6)
      !
      ! Write out the grid.
      WRITE(100,300) ((2.0_RWPF*PI*REAL(i, RWPF)/REAL(nxW, RWPF), &
                       i=1,nxW), j=1,nyW)
      WRITE(100,300) ((2.0_RWPF*PI*REAL(j, RWPF)/REAL(nyW, RWPF), &
                       i=1,nxW), j=1,nyW)
      300 FORMAT (ES15.8)
      !
      ! Close the grid file.
      CLOSE(UNIT=100)
   END SUBROUTINE WriteGridPlot3d

   !> Routine to write a solution file in Plot3d format.
   !!
   !! We use the uC and uR arrays as scratch for transforming omega and then the
   !! streamfunction. After omega and the streamfunction are written, we
   !! calculate u and v, transform them, and write them out to file.
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
   !> @param[in] cS Current stage in Q to use.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] uR Real cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[in,out] vR Real cast of v velocity array.
   !> @param[in] Qc Complex cast of Q array.
   !> @param[in] Qr Real cast of Q array.
   !> @param[in] fname Name for the output file.
   SUBROUTINE WriteRestartPlot3d(rISize, cISize, nxW, nyW, &
                                 kxLT, kyLT, kxMT, kyMT, &
                                 kxLU, kyLU, kxMU, kyMU, &
                                 nVar, nStg, cS, &
                                 uC, uR, vC, vR, Qc, Qr, fname)
      ! Required modules.
      USE Spectral_m,ONLY: TransformC2R, DuDx, DuDy
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize, nxW, nyW
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLT, kyLT, kxMT, kyMT
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLU, kyLU, kxMU, kyMU
      INTEGER(KIND=IWPF),INTENT(IN) :: nVar, nStg, cS
      COMPLEX(KIND=CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(INOUT) :: uR, vR
      COMPLEX(CWPC),DIMENSION(kxLT:kxMT,kyLT:kyMT,nVar,nStg),INTENT(IN) :: Qc
      REAL(KIND=RWPC),DIMENSION(rISize,nyW,nVar,nStg),INTENT(IN) :: Qr
      CHARACTER(LEN=*),INTENT(IN) :: fname
      ! Local variables.
      ! Looping indices for physical space.
      INTEGER(KIND=IWPF) :: i, j
      ! Looping indices for spectral space.
      INTEGER(KIND=IWPF) :: kx, ky
      ! Masked wavenumber for differentiation in the y direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kyMask
      ! Masked wavenumber for differentiation in the x direction.
      INTEGER(KIND=IWPF),DIMENSION(1) :: kxMask

      ! Open the output file for writing and fill in the header.
      OPEN(UNIT=100,FILE=TRIM(fname),STATUS='REPLACE',FORM='FORMATTED')
      WRITE(100,150) 1_IWPF
      WRITE(100,200) nxW, nyW
      150 FORMAT (I6.6)
      200 FORMAT (I6.6,1X,I6.6)

      ! 1. Write out the vorticity to file.
      !
      ! Zero out the scratch array.
      uC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Fill in the useful data for transforming vorticity.
      uC(kxLU:kxMU,kyLU:kyMU) = Qc(kxLU:kxMU,kyLU:kyMU,1,cS)
      !
      ! Transform the vorticity to physical space.
      CALL TransformC2R(rISize, cISize, nyW, uC, uR)
      !
      ! Write to file.
      WRITE(100,300) ((uR(i,j), i=1,nxW), j=1,nyW)
      300 FORMAT (ES15.8)

      ! 2. Write out the streamfunction to file.
      !
      ! Zero out the scratch array.
      uC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Fill in the useful data for transforming the streamfunction.
      uC(kxLU:kxMU,kyLU:kyMU) = Qc(kxLU:kxMU,kyLU:kyMU,2,cS)
      !
      ! Transform the streamfunction to physical space.
      CALL TransformC2R(rISize, cISize, nyW, uC, uR)
      !
      ! Write to file.
      WRITE(100,300) ((uR(i,j), i=1,nxW), j=1,nyW)

      ! 3. Write out the u velocity to file.
      !
      ! Zero out the working array.
      uC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Differentiate the streamfunction to get the u velocity.
      kyMask(1) = kyMU
      CALL DuDy(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kyMask, Qc(:,:,2,cS), uC)
      !
      ! Invert the velocity to physical space.
      CALL TransformC2R(rISize, cISize, nyW, uC, uR)
      !
      ! Write to file.
      WRITE(100,300) ((uR(i,j), i=1,nxW), j=1,nyW)

      ! 3. Write out the v velocity to file.
      !
      ! Zero out the working array.
      vC(:,:) = (0.0_RWPC, 0.0_RWPC)
      !
      ! Differentiate the streamfunction to get the v velocity.
      kxMask(1) = kxMU
      CALL DuDx(kxLT, kxMT, kyLT, kyMT, kxLU, kxMU, kyLU, kyMU, &
                1_IWPF, kxMask, Qc(:,:,2,cS), vC)
      !
      ! Invert the velocity to physical space.
      CALL TransformC2R(rISize, cISize, nyW, vC, vR)
      !
      ! Write to file.
      WRITE(100,300) ((vR(i,j), i=1,nxW), j=1,nyW)

      ! Close the restart file.
      CLOSE(UNIT=100)
   END SUBROUTINE WriteRestartPlot3d

END MODULE IO_m

