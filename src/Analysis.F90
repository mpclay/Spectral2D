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
!> @file Analysis.F90
!> @author Matthew Clay
!> @brief Module providing analysis routines.
MODULE Analysis_m

   ! Required modules.
   USE Parameters_m,ONLY: IWPC, IWPF, RWPC, RWPF, CWPC
   USE IO_m,ONLY: FILE_NAME_LENGTH

   IMPLICIT NONE

   !> Name of the time-series data analysis output file.
   CHARACTER(LEN=FILE_NAME_LENGTH),PRIVATE :: eulOut

   ! Module procedures.
   PUBLIC :: SetupAnalysis, Analysis, ComputeSpectrum, BinarySearch

CONTAINS

   !> Procedure to set up the data analysis module.
   !!
   !> @param[in] fRoot Desired root name for the Eulerian statistics file.
   SUBROUTINE SetupAnalysis(fRoot)
      ! Required modules.
      USE IO_m,ONLY: FileName
      IMPLICIT NONE
      ! Calling arguments.
      CHARACTER(LEN=*),INTENT(IN) :: fRoot

      ! Form the Eulerian file name.
      CALL FileName(fRoot, 'DAT', eulOut)

      ! Open the Eulerian statistics file and write the header.
      OPEN(UNIT=100,FILE=eulOut,STATUS='REPLACE',ACTION='WRITE', &
           FORM='FORMATTED')
      WRITE(100,200) '#', 'NADV', 'TIME', 'TKE', 'EPS', 'EnstDiss', &
                     'xSkew', 'ySkew'
      200 FORMAT (A,T3,A,T17,A,T33,A,T49,A,T65,A,T81,A,T97,A)
      !
      ! Close out the file.
      CLOSE(UNIT=100)
   END SUBROUTINE SetupAnalysis

   !> Procedure to perform data analysis during simulations.
   !!
   !! This routine will calculate the two-dimensional skewness factors for the
   !! x and y directions, the total turbulent kinetic energy, the enstrophy
   !! dissipation rate, and the energy/enstrophy dissipation spectra.
   !!
   !> @param[in] rank MPI rank of this process.
   !> @param[in] nadv Current simulation step.
   !> @param[in] time Current simulation time.
   !> @param[in] nu Physical viscosity.
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction on this process.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] rISize Size of real data in the i direction.
   !> @param[in] cISize Size of complex data in the i direction.
   !> @param[in] kxP Array of wavenumbers in the x direction on this process.
   !> @param[in] kyP Array of wavenumbers in the y direction on this process.
   !> @param[in] kxLG Lowest x wavenumber in the total simulation.
   !> @param[in] kyLG Lowest y wavenumber in the total simulation.
   !> @param[in] kxMG Maximum x wavenumber in the total simulation.
   !> @param[in] kyMG Maximum y wavenumber in the total simulation.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] uR Real cast of the u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[in,out] vR Real cast of the v velocity array.
   !> @param[in,out] wC Complex cast of the vorticity field.
   !> @param[in,out] wR Real cast of the vorticity field.
   !> @param[in,out] psiC Complex cast of the streamfunction.
   !> @param[in,out] psiR Real cast of the streamfunction.
   SUBROUTINE Analysis(rank, nadv, time, nu, nxG, nyG, nxP, nyP, &
                       rISize, cISize, kxP, kyP, kxLG, kyLG, kxMG, kyMG, &
                       uC, uR, vC, vR, wC, wR, psiC, psiR)
      ! Required modules.
      USE ISO_FORTRAN_ENV,ONLY: OUTPUT_UNIT
      USE MPI
      USE Spectral_m,ONLY: ComputeVelocity, DuDx, DuDy, &
                           TransformR2C, TransformC2R
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rank, nadv, nxG, nyG, nxP, nyP
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize
      REAL(KIND=RWPF),INTENT(IN) :: time, nu
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLG, kyLG, kxMG, kyMG
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC, wC, psiC
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: uR, vR, wR, psiR
      ! Local variables.
      ! Send/receive buffers for calculating skewness.
      REAL(KIND=RWPC),DIMENSION(3) :: sBuff, rBuff
      ! Skewness calculated in the x direction.
      REAL(KIND=RWPC) :: xSkew
      ! Skewness calculated in the y direction.
      REAL(KIND=RWPC) :: ySkew
      ! Normalization factor for averaging.
      REAL(KIND=RWPC) :: fact
      ! Total turbulent kinetic energy.
      REAL(KIND=RWPF) :: tke
      ! Mean dissipation rate.
      REAL(KIND=RWPF) :: eps
      ! Mean enstrophy dissipation rate.
      REAL(KIND=RWPF) :: enstDiss
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! Inform the user that Eulerian stats are being collected.
      IF (rank == 0_IWPF) THEN
         WRITE(OUTPUT_UNIT,100) 'Writing Eulerian statistics at t = ', time, &
                                ' and nadv = ', nadv
         100 FORMAT (A,ES15.8,A,I8.8)
      END IF

      ! Calculate the skewness factor for the x direction.
      !
      ! NOTE: while doing this, the v velocity array is used as a scratch array
      ! for the x derivative of vorticity.
      !
      ! 1. Calculate the velocity in Fourier space. NOTE: this is currently
      !    inefficient, since the v velocity is not being used here.
      CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                           uC, uR, vC, vR, PsiC, .FALSE.)
      !
      ! 2. Differentiate the u velocity and vorticity with respect to x. NOTE:
      !    after this step, dw/dx is stored in the v velocity array!
      CALL DuDx(cISize, nyP, kxP, kyP, uC, uC)
      CALL DuDx(cISize, nyP, kxP, kyP, wC, vC)
      !
      ! 3. Transform the signals to physical space.
      CALL TransformC2R(rISize, cISize, nyP, uC, uR)
      CALL TransformC2R(rISize, cISize, nyP, vC, vR)
      !
      ! 4. Form the summations required for the skewness with the data on this
      !    process.
      sBuff(:) = 0.0_RWPC
      DO j = 1, nyP
         DO i = 1, nxP
            sBuff(1) = sBuff(1) + uR(i,j)*vR(i,j)**2
            sBuff(2) = sBuff(2) + uR(i,j)**2
            sBuff(3) = sBuff(3) + vR(i,j)**2
         END DO
      END DO
      !
      ! 5. Reduce all summations to the root process.
      CALL MPI_REDUCE(sBuff, rBuff, 3, MPI_DOUBLE, MPI_SUM, 0, &
                      MPI_COMM_WORLD, ierr)
      !
      ! 6. The root process performs the final averaging for the skewness.
      IF (rank == 0_IWPF) THEN
         ! Average each of the summations.
         fact = REAL(nxG, RWPC)*REAL(nyG, RWPC)
         rBuff(1) = rBuff(1)/fact
         rBuff(2) = rBuff(2)/fact
         rBuff(3) = rBuff(3)/fact
         !
         ! Form the x direction skewness.
         xSkew = -2.0_RWPC*rBuff(1)/(SQRT(rBuff(2))*rBuff(3))
      END IF

      ! Calculate the skewness factor for the y direction.
      !
      ! NOTE: while doing this, the u velocity array is used as a scratch array
      ! for the y derivative of vorticity.
      !
      ! 1. Calculate the velocity in Fourier space. NOTE: this is currently
      !    inefficient, since the u velocity is not being used here.
      CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                           uC, uR, vC, vR, PsiC, .FALSE.)
      !
      ! 2. Differentiate the v velocity and vorticity with respect to y. NOTE:
      !    after this step, dw/dy is stored in the u velocity array!
      CALL DuDy(cISize, nyP, kxP, kyP, vC, vC)
      CALL DuDy(cISize, nyP, kxP, kyP, wC, uC)
      !
      ! 3. Transform the signals to physical space.
      CALL TransformC2R(rISize, cISize, nyP, vC, vR)
      CALL TransformC2R(rISize, cISize, nyP, uC, uR)
      !
      ! 4. Form the summations required for the skewness with the data on this
      !    process.
      sBuff(:) = 0.0_RWPC
      DO j = 1, nyP
         DO i = 1, nxP
            sBuff(1) = sBuff(1) + vR(i,j)*uR(i,j)**2
            sBuff(2) = sBuff(2) + vR(i,j)**2
            sBuff(3) = sBuff(3) + uR(i,j)**2
         END DO
      END DO
      !
      ! 6. Reduce all summations to the root process.
      CALL MPI_REDUCE(sBuff, rBuff, 3, MPI_DOUBLE, MPI_SUM, 0, &
                      MPI_COMM_WORLD, ierr)
      !
      ! 7. The root process performs the final averaging for the skewness.
      IF (rank == 0_IWPF) THEN
         ! Average each of the summations.
         fact = REAL(nxG, RWPC)*REAL(nyG, RWPC)
         rBuff(1) = rBuff(1)/fact
         rBuff(2) = rBuff(2)/fact
         rBuff(3) = rBuff(3)/fact
         !
         ! Form the y direction skewness.
         ySkew = -2.0_RWPC*rBuff(1)/(SQRT(rBuff(2))*rBuff(3))
      END IF

      ! Calculate the energy, dissipation, and enstrophy dissipation spectra.
      !
      ! NOTE: we also get the total TKE and enstrophy dissipation from this.
      CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxP, kyP, &
                           uC, uR, vC, vR, PsiC, .FALSE.)
      CALL ComputeSpectrum(rank, nu, nxG, nyG, nxP, nyP, rISize, cISize, &
                           kxP, kyP, kxLG, kyLG, kxMG, kyMG, .TRUE., &
                           nadv, uC, vC, tke, eps, enstDiss)

      ! Write out the data to file.
      OPEN(UNIT=100,FILE=eulOut,STATUS='OLD',ACTION='WRITE', &
           FORM='FORMATTED',ACCESS='APPEND')
      WRITE(100,150) nadv, time, tke, eps, enstDiss, xSkew, ySkew
      150 FORMAT (I12.10,6ES16.8)
      CLOSE(UNIT=100)
   END SUBROUTINE Analysis

   !> Procedure to compute the energy spectrum.
   !!
   !> @param[in] rank MPI rank of this process.
   !> @param[in] nu Physical viscosity of the fluid.
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction on this process.
   !> @param[in] nyP Number of grid points in the y direction on this process.
   !> @param[in] rISize Size of real data in the i direction.
   !> @param[in] cISize Size of complex data in the i direction.
   !> @param[in] kxP Array of wavenumbers in the x direction on this process.
   !> @param[in] kyP Array of wavenumbers in the y direction on this process.
   !> @param[in] kxLG Lowest x wavenumber in the total simulation.
   !> @param[in] kyLG Lowest y wavenumber in the total simulation.
   !> @param[in] kxMG Maximum x wavenumber in the total simulation.
   !> @param[in] kyMG Maximum y wavenumber in the total simulation.
   !> @param[in] output Whether or not to write the spectrum to file.
   !> @param[in] outNum Number for the output spectrum.
   !> @param[in,out] uC Complex cast of u velocity array.
   !> @param[in,out] vC Complex cast of v velocity array.
   !> @param[out] tke Total turbulent kinetic energy. Only for root process.
   !> @param[out] eps Dissipation rate. Only for root process.
   !> @param[out] enstDiss Enstrophy dissipation rate. Only for root process.
   SUBROUTINE ComputeSpectrum(rank, nu, nxG, nyG, nxP, nyP, rISize, cISize, &
                              kxP, kyP, kxLG, kyLG, kxMG, kyMG, output, &
                              outNum, uC, vC, tke, eps, enstDiss)
      ! Required modules.
      USE MPI
      USE IO_m,ONLY: FILE_NAME_LENGTH, FileName
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: rank, nxG, nyG, nxP, nyP, rISize, cISize
      REAL(KIND=RWPF),INTENT(IN) :: nu
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxP
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyP
      INTEGER(KIND=IWPF),INTENT(IN) :: kxLG, kyLG, kxMG, kyMG, outNum
      LOGICAL,INTENT(IN) :: output
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC
      REAL(KIND=RWPF),INTENT(OUT) :: tke, eps, enstDiss
      ! Local variables.
      ! Maximum wavenumber in the simulation.
      REAL(KIND=RWPF) :: kMax
      ! Maximum integer wavenumber bin we will use.
      INTEGER(KIND=IWPF) :: kIntMax
      ! Width of the bins.
      REAL(KIND=RWPF) :: dk
      ! Number of bins being used.
      INTEGER(KIND=IWPF) :: numBin
      ! Index for a mode in the bin array.
      INTEGER(KIND=IWPF) :: indx
      ! Array to bin kinetic energy locally.
      REAL(KIND=RWPF),DIMENSION(:),ALLOCATABLE :: EkBinL
      ! Array to reduce kinetic energy spectrum to the root process.
      REAL(KIND=RWPF),DIMENSION(:),ALLOCATABLE :: EkBinG
      ! Array to count how many modes contribute to a bin.
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: binCntL
      ! Array to reduce the bin count to the root process.
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: binCntG
      ! Number of bin edges.
      INTEGER(KIND=IWPF) :: numEdges
      ! Array of bin edges for the binning routine.
      REAL(KIND=RWPF),DIMENSION(:),ALLOCATABLE :: binEdges
      ! Bin centers, only used on root process for output.
      INTEGER(KIND=IWPF),DIMENSION(:),ALLOCATABLE :: kBin
      ! Integer versions of the kx and ky wavenumbers.
      INTEGER(KIND=IWPF) :: kxI, kyI
      ! Real versions of the kx and ky wavenumbers.
      REAL(KIND=RWPF) :: kxR, kyR
      ! Magnitude of the wavenumber vector.
      REAL(KIND=RWPF) :: kMag
      ! Twice the energy content of a Fourier mode.
      REAL(KIND=RWPF) :: Ehat
      ! Number to increase bin count by for a given mode.
      INTEGER(KIND=IWPF) :: cnt
      ! File name for the spectrum if it is being saved.
      CHARACTER(LEN=FILE_NAME_LENGTH) :: fname
      ! Looping indices.
      INTEGER(KIND=IWPF) :: i, j
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! The maximum wavenumber in the simulation.
      kMax = SQRT(REAL(kxMG*kxMG + kyMG*kyMG, RWPF))
      !
      ! The maximum integer wavenumber bin center we will use.
      kIntMax = INT(kMax + 1.0_RWPC, IWPF)
      !
      ! Width of the bins will just be 1 for now.
      dk = 1.0_RWPF
      !
      ! Number of bins being used.
      numBin = kIntMax + 1_IWPF
      !
      ! Allocate local and global storage for the kinetic energy.
      ALLOCATE(EkBinL(numBin))
      EkBinL(:) = 0.0_RWPF
      IF (rank == 0) THEN
         ALLOCATE(EkBinG(numBin))
         EkBinG(:) = 0.0_RWPF
      END IF
      !
      ! Allocate memory for the bin counters.
      ALLOCATE(BinCntL(numBin))
      BinCntL(:) = 0_IWPF
      IF (rank == 0) THEN
         ALLOCATE(BinCntG(numBin))
         BinCntG(:) = 0_IWPF
      END IF
      !
      ! Number of bin edges there will be.
      numEdges = numBin + 1_IWPF
      !
      ! Fill in the wavenumber bin edges, which are used for the binning code.
      ALLOCATE(binEdges(numEdges))
      binEdges(1) = 0.0_RWPF
      DO i = 1, kIntMax
         binEdges(i+1) = 0.5_RWPF + REAL(i-1, RWPF)*1.0_RWPF
      END DO
      binEdges(numEdges) = REAL(kIntMax, RWPF) + 0.1_RWPF

      ! Loop over all wavenumbers on this process and bin the kinetic energy.
      DO j = 1, nyP
         !
         ! The wavenumber in the y direction.
         kyI = kyP(j)
         kyR = REAL(kyP(j), RWPF)
         DO i = 1, cISize
            !
            ! The wavenumber in the x direction.
            kxI = kxP(i)
            kxR = REAL(kxP(i), RWPF)
            !
            ! The magnitude of the wavenumber vector.
            kMag = SQRT(kxR**2 + kyR**2)
            !
            ! Calculate the kinetic energy for this mode.
            Ehat = REAL(uC(i,j)*CONJG(uC(i,j)) + vC(i,j)*CONJG(vC(i,j)), RWPF)
            !
            ! Determine which bin the mode is in.
            CALL BinarySearch(numEdges, binEdges, kMag, indx)
            !
            ! When calculating the energy contribution of a mode, we must take
            ! into account what data is stored by FFTW. Because of Hermitian
            ! symmetry, when kx is between (1,Nx/2-1), we double the mode's
            ! contribution to take into account the modes not stored in memory.
            IF ((kxI == 0_IWPF) .OR. (kxI == kxMg)) THEN
               EkBinL(indx) = EkBinL(indx) + 0.5_RWPF*Ehat
               cnt = 1_IWPF
            ELSE
               EkBinL(indx) = EkBinL(indx) + Ehat
               cnt = 2_IWPF
            END IF
            !
            ! Increase the bin counter for this bin.
            binCntL(indx) = binCntL(indx) + cnt
         END DO
      END DO

      ! Reduce the binned kinetic energy and counters to root.
      CALL MPI_REDUCE(EkBinL, EkBinG, numBin, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(binCntL, binCntG, numBin, MPI_INTEGER, MPI_SUM, 0, &
                      MPI_COMM_WORLD, ierr)

      ! The root process can now go through and calculate the energy spectrum.
      IF (rank == 0) THEN
         ! Allocate memory for the bin centers.
         ALLOCATE(kBin(numBin))
         !
         ! Fill in the bin vector and normalize the data by the bin width.
         DO i = 1, numBin
            ! The bin centers are located at the integers.
            kBin(i) = i - 1_IWPF
            !
            ! The energy spectrum must be normalized by the bin width.
            EkBinG(i) = EkBinG(i)/dk
         END DO
         !
         ! Save the spectrum to file if instructed to do so.
         IF (output) THEN
            CALL FileName('SPECTRUM', outNum, 'DAT', fname)
            OPEN(UNIT=10,FILE=fname,STATUS='REPLACE',FORM='FORMATTED')
            WRITE(10,50) '#', 'k', 'E(k)', 'D(k)', 'DEnst(k)'
            DO i = 1, numBin
               WRITE(10,100) kBin(i), EkBinG(i), &
                             REAL(kBin(i), RWPF)**2*EkBinG(i), &
                             REAL(kBin(i), RWPF)**4*EkBinG(i)
            END DO
            50 FORMAT (A,T3,A,T11,A,T28,A,T45,A)
            100 FORMAT (I6.6,ES16.8,1X,ES16.8,1X,ES16.8)
            CLOSE(UNIT=10)
         END IF
         !
         ! Root process calculates the tke, dissipation rate, and the enstrophy
         ! dissipation rate.
         tke = 0.0_RWPF
         eps = 0.0_RWPF
         enstDiss = 0.0_RWPF
         DO i = 1, numBin
            tke = tke + EkBinG(i)*dk
            eps = eps + REAL(kBin(i), RWPF)**2*EkBinG(i)*dk
            enstDiss = enstDiss + REAL(kBin(i), RWPF)**4*EkBinG(i)*dk
         END DO
         eps = 2.0_RWPF*nu*eps
         enstDiss = 2.0_RWPF*nu*enstDiss
      END IF

      ! Free temporary memory.
      DEALLOCATE(EkBinL)
      DEALLOCATE(BinCntL)
      DEALLOCATE(binEdges)
      IF (rank == 0) THEN
         DEALLOCATE(EkBinG)
         DEALLOCATE(BinCntG)
         DEALLOCATE(kBin)
      END IF
   END SUBROUTINE ComputeSpectrum

   !> Binary search routine.
   !!
   !! This routine returns an index i such that f is in the interval
   !! [fVec(i),fVec(i+1)) (note the open bracket on the right). If the value is
   !! exactly equal to fVec(end), it is placed in the (end-1) bin, i.e., the
   !! open bracket is closed for the last bin in the interval.
   !!
   !! Note that we assume the fVec array is monotone increasing.
   !!
   !> @param[in] n Length of the vector fVec.
   !> @param[in] fVec Value of the function at the bin endpoints.
   !> @param[in] f Value of the function used for searching.
   !> @param[out] indx Index of f in fVec.
   SUBROUTINE BinarySearch(n, fVec, f, indx)
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: n
      REAL(KIND=RWPF),DIMENSION(n),INTENT(IN) :: fVec
      REAL(KIND=RWPF),INTENT(IN) :: f
      INTEGER(KIND=IWPF),INTENT(OUT) :: indx
      ! Local variables.
      ! Lower index in the reduced interval for the binary search algorithm.
      INTEGER(KIND=IWPF) :: iL
      ! Upper index in the reduced interval for the binary search algorithm.
      INTEGER(KIND=IWPF) :: iU
      ! Middle index of the interval [iL, iU].
      INTEGER(KIND=IWPF) :: iM

      ! Initialize iL and iU. Note that iL is initialized below the actual
      ! interval, and that iU is initialized to the exact length of the
      ! interval. With this iU initialization (and the following logic), the
      ! last bin will be a closed interval, as discussed above.
      iL = 0_IWPF
      iU = n

      ! Loop until we get the index for this f.
      DO WHILE ((iU - iL) > 1_IWPF)
         ! Calculate the middle index.
         iM = (iL + iU)/2_IWPF
         !
         ! If f is greater than fVec(iM), the new interval becomes [iM, iU],
         ! otherwise, the interval is [iL, iM].
         IF (f >= fVec(iM)) THEN
            iL = iM
         ELSE
            iU = iM
         END IF
      END DO

      ! Return the index to the user.
      indx = iL
   END SUBROUTINE BinarySearch

END MODULE Analysis_m

