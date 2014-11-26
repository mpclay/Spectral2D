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
!!
!! This module provides routines to write out restart files and grid files for
!! the simulation. If HDF5 output is used, the files will be written in
!! parallel.
!!
!! TODO:
!!
!!    1. Add Plot3d format, both ASCII and binary formatted.
MODULE IO_m

   ! Required modules.
   USE HDF5
   USE Parameters_m,ONLY: IWPF, IWPC, RWPF, RWPC, CWPC

   IMPLICIT NONE

   ! Supported output file types.
   !
   !> HDF5 file output.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: HDF5_OUTPUT = 1_IWPF

   !> File name length for output files.
   INTEGER(KIND=IWPF),PARAMETER,PUBLIC :: FILE_NAME_LENGTH = 256_IWPF

   !> Interface for the file name function.
   INTERFACE FileName
      MODULE PROCEDURE FileNameWithoutNum
      MODULE PROCEDURE FileNameWithNum
   END INTERFACE FileName

   !> Interface for writing attributes to HDF5 files.
   INTERFACE WriteAttribute
      MODULE PROCEDURE WriteAttributeInteger
   END INTERFACE WriteAttribute

   ! Module procedures.
   PUBLIC :: WriteRestart, WriteGridHDF5
   PRIVATE :: WriteRestartHDF5
   PRIVATE :: WritePHDF5Dataset, WriteAttributeInteger
   PRIVATE :: FileNameWithoutNum, FileNameWithNum

CONTAINS

   !> Write a restart file.
   !!
   !! In this routine we have to perform some transformations to get all of the
   !! data into physical space, which is the desired output format for data
   !! visualization.
   !!
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction for this process.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] i1 Starting i index for the data on this process.
   !> @param[in] i2 Ending i index for the data on this process.
   !> @param[in] j1 Starting j index for the data on this process.
   !> @param[in] j2 Ending j index for the data on this process.
   !> @param[in] rISize Size of the i dimension for real data arrays.
   !> @param[in] cISize Size of the i dimension for complex data arrays.
   !> @param[in] kxG Wavenumbers in the x direction.
   !> @param[in] kyG Wavenumbers in the y direction.
   !> @param[in] rank MPI process ID for this process.
   !> @param[in] outNum Number of the restart file.
   !> @param[in] outType Desired output file type.
   !> @param[in] uC Complex cast of u velocity data.
   !> @param[in] uR Real cast of u velocity data.
   !> @param[in] vC Complex cast of v velocity data.
   !> @param[in] vR Real cast of v velocity data.
   !> @param[in] wC Complex cast of vorticity.
   !> @param[in] wR Real cast of vorticity.
   !> @param[in] psiC Complex cast of the streamfunction.
   !> @param[in] psiR Real cast of the streamfunction.
   SUBROUTINE WriteRestart(nxG, nyG, nxP, nyP, i1, i2, j1, j2, &
                           rISize, cISize, kxG, kyG, rank, outNum, outType, &
                           uC, uR, vC, vR, wC, wR, psiC, psiR)
      ! Required modules.
      USE MPI,ONLY: MPI_COMM_WORLD
      USE Spectral_m,ONLY: ComputeVelocity, TransformC2R, TransformR2C
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nxP, nyP, i1, i2, j1, j2
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, cISize, rank, outNum, outType
      INTEGER(KIND=IWPF),DIMENSION(cISize),INTENT(IN) :: kxG
      INTEGER(KIND=IWPF),DIMENSION(nyP),INTENT(IN) :: kyG
      COMPLEX(KIND=CWPC),DIMENSION(cISize,nyP),INTENT(INOUT) :: uC, vC, wC, psiC
      REAL(KIND=RWPC),DIMENSION(rISize,nyP),INTENT(INOUT) :: uR, vR, wR, psiR
      ! Local variables.
      ! Output file name.
      CHARACTER(LEN=FILE_NAME_LENGTH) :: fname

      ! Get variables ready to write out the initial restart file.
      CALL ComputeVelocity(nxP, nyP, rISize, cISize, kxG, kyG, &
                           uC, uR, vC, vR, psiC, .TRUE.)
      CALL TransformC2R(rISize, cISize, nyG, wC, wR)
      CALL TransformC2R(rISize, cISize, nyG, psiC, psiR)
      !
      ! Write out the restart file.
      SELECT CASE (outType)
         CASE (HDF5_OUTPUT)
            CALL FileName('REST', outNum, 'h5', fname)
            CALL WriteRestartHDF5(MPI_COMM_WORLD, rank, fname, rISize, &
                                  nxG, nyG, nxP, nyP, i1, j1, &
                                  wR, psiR, uR, vR)
         CASE DEFAULT
      END SELECT
      !
      ! Transform back to spectral space.
      CALL TransformR2C(nxG, nyG, rISize, cISize, wR, wC)
      CALL TransformR2C(nxG, nyG, rISize, cISize, psiR, psiC)
   END SUBROUTINE WriteRestart

   !> Write an HDF5 grid file.
   !!
   !> @param[in] comm MPI communicator from the calling code.
   !> @param[in] rank MPI rank in the process listing.
   !> @param[in] fName Name of the grid file.
   !> @param[in] nxG Total number of grid points in the x direction.
   !> @param[in] nyG Total number of grid points in the y direction.
   !> @param[in] nxP Number of grid points in the x direction for this process.
   !> @param[in] nyP Number of grid points in the y direction for this process.
   !> @param[in] i1 Starting i index for this process.
   !> @param[in] i2 Ending i index for this process.
   !> @param[in] j1 Starting j index for this process.
   !> @param[in] j2 Ending j index for this process.
   SUBROUTINE WriteGridHDF5(comm, rank, fName, nxG, nyG, nxP, nyP, &
                            i1, i2, j1, j2)
      ! Required modules.
      USE HDF5
      USE MPI,ONLY: MPI_INFO_NULL
      USE Parameters_m,ONLY: PI
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: comm, rank
      CHARACTER(LEN=*),INTENT(IN) :: fName
      INTEGER(KIND=IWPF),INTENT(IN) :: nxG, nyG, nxP, nyP, i1, i2, j1, j2
      ! Local variables.
      ! HDF5 ID for the restart file.
      INTEGER(KIND=HID_T) :: file_id
      ! HDF5 ID for the header dataset.
      INTEGER(KIND=HID_T) :: dataset_id
      ! HDF5 ID for the FlowData group.
      INTEGER(KIND=HID_T) :: group_id
      ! HDF5 data space identifier.
      INTEGER(KIND=HID_T) :: dataspace_id
      ! Property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! Needed whem making the property list.
      INTEGER(KIND=IWPF) :: info
      ! Offset for the hyperslab of data taken from the restart file.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: offset
      ! Number of blocks to select from the dataspace for hyperslabs.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: counter
      ! Strides to take in the dataspace for the hyperslab.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: stride
      ! Size of the hyperslab block for this process.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: block
      ! Total size of the restart file.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dimsT
      ! X grid coordinates for this process.
      REAL(KIND=RWPC),DIMENSION(:,:),ALLOCATABLE :: x
      ! Y grid coordinates for this process.
      REAL(KIND=RWPC),DIMENSION(:,:),ALLOCATABLE :: y
      ! Looping indices for i and j directions.
      INTEGER(KIND=IWPF) :: i, j
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! Allocate memory for the grid points.
      ALLOCATE(x(i1:i2,j1:j2))
      ALLOCATE(y(i1:i2,j1:j2))
      !
      ! Fill in the grid points for this process.
      DO j = j1, j2
         DO i = i1, i2
            x(i,j) = REAL(i, RWPC)*2.0_RWPC*PI/REAL(nxG, RWPC)
            y(i,j) = REAL(j, RWPC)*2.0_RWPC*PI/REAL(nyG, RWPC)
         END DO
      END DO

      ! The offset for the data in the HDF5 file are based on the starting
      ! global indices for this process. We subtract 1 to get [0, 0, 0] for the
      ! first process, etc.
      offset = [INT(i1, HSIZE_T) - 1_HSIZE_T, &
                INT(j1, HSIZE_T) - 1_HSIZE_T]

      ! For each hyperslab, we will just select 1 block from the dataspace.
      counter = [1_HSIZE_T, 1_HSIZE_T]

      ! We aren't doing any striding with the data--each process has a
      ! contiguous chunk of the data.
      stride = [1_HSIZE_T, 1_HSIZE_T]

      ! Size of the data for this process.
      block = [INT(nxP, HSIZE_T), INT(nyP, HSIZE_T)]

      ! Dimensions of the complete restart file.
      dimsT = [INT(nxG, HSIZE_T), INT(nyG, HSIZE_T)]

      ! Initialize the Fortran-HDF5 interface.
      CALL H5OPEN_F(ierr)

      ! Set up file access property list with parallel I/O access.
      info = MPI_INFO_NULL
      CALL H5PCREATE_F(H5P_FILE_ACCESS_F, plist_id, ierr)
      CALL H5PSET_FAPL_MPIO_F(plist_id, comm, info, ierr)

      ! Open the file for writing.
      CALL H5FCREATE_F(fName, H5F_ACC_TRUNC_F, file_id, ierr, &
                       access_prp=plist_id)

      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)

      ! Create the 'Header' dataspace.
      CALL H5SCREATE_F(H5S_NULL_F, dataspace_id, ierr)

      ! Create the 'Header' dataset.
      CALL H5DCREATE_F(file_id, 'Header', H5T_NATIVE_INTEGER, dataspace_id, &
                       dataset_id, ierr)

      ! Write header attributes.
      CALL WriteAttribute(dataset_id, H5T_NATIVE_INTEGER, 'nx', nxG)
      CALL WriteAttribute(dataset_id, H5T_NATIVE_INTEGER, 'ny', nyG)

      ! Close the 'Header' dataset
      CALL H5DCLOSE_F(dataset_id, ierr)

      ! Close the 'Header' dataspace
      CALL H5SCLOSE_F(dataspace_id, ierr)

      ! Create the 'FlowData' data group
      CALL H5GCREATE_F(file_id, 'Domain_00001', group_id, ierr)

      ! Write the data to the file using collective writeout.
      CALL WritePHDF5Dataset(nxP, nyP, group_id, 'x', x(i1:i2,j1:j2), &
                             offset, counter, stride, block, dimsT)
      CALL WritePHDF5Dataset(nxP, nyP, group_id, 'y', y(i1:i2,j1:j2), &
                             offset, counter, stride, block, dimsT)

      ! Close the 'FlowData' data group.
      CALL H5GCLOSE_F(group_id, ierr)

      ! Close the file.
      CALL H5FCLOSE_F(file_id, ierr)

      ! Close the Fortran-HDF5 interface.
      CALL H5CLOSE_F(ierr)

      ! Free temporary memory.
      DEALLOCATE(x)
      DEALLOCATE(y)
   END SUBROUTINE WriteGridHDF5

   !> Write an HDF5 restart file.
   !!
   !> @param[in] comm MPI communicator from the calling code.
   !> @param[in] rank MPI rank in the main process list.
   !> @param[in] fName Name of the restart file.
   !> @param[in] rISize Size of real data arrays in the i direction.
   !> @param[in] nxW Number of cells in the x direction for the simulation.
   !> @param[in] nyW Number of cells in the y direction for the simulation.
   !> @param[in] nxP Number of cells in the x direction for this process.
   !> @param[in] nyP Number of cells in the y direction for this process.
   !> @param[in] i1 Starting i index for this process.
   !> @param[in] j1 Starting j index for this process.
   !> @param[in] w Vorticity in the domain.
   !> @param[in] psi Streamfunction in the domain.
   !> @param[in] u Velocity in the x direction.
   !> @param[in] v Velocity in the y direction.
   SUBROUTINE WriteRestartHDF5(comm, rank, fName, rISize, nxW, nyW, &
                               nxP, nyP, i1, j1, w, psi, u, v)
      ! Required modules.
      USE HDF5
      USE MPI,ONLY: MPI_INFO_NULL
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: comm, rank
      CHARACTER(LEN=*),INTENT(IN) :: fName
      INTEGER(KIND=IWPF),INTENT(IN) :: rISize, nxW, nyW, nxP, nyP, i1, j1
      REAL(KIND=RWPC),DIMENSION(rISize,nyW),INTENT(IN) :: w, psi, u, v
      ! Local variables.
      ! HDF5 ID for the restart file.
      INTEGER(KIND=HID_T) :: file_id
      ! HDF5 ID for the header dataset.
      INTEGER(KIND=HID_T) :: dataset_id
      ! HDF5 ID for the FlowData group.
      INTEGER(KIND=HID_T) :: group_id
      ! HDF5 data space identifier.
      INTEGER(KIND=HID_T) :: dataspace_id
      ! Property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! Needed whem making the property list.
      INTEGER(KIND=IWPF) :: info
      ! Offset for the hyperslab of data taken from the restart file.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: offset
      ! Number of blocks to select from the dataspace for hyperslabs.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: counter
      ! Strides to take in the dataspace for the hyperslab.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: stride
      ! Size of the hyperslab block for this process.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: block
      ! Total size of the restart file.
      INTEGER(KIND=HSIZE_T),DIMENSION(2) :: dimsT
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! The offset for the data in the HDF5 file are based on the starting
      ! global indices for this process. We subtract 1 to get [0, 0, 0] for the
      ! first process, etc.
      offset = [INT(i1, HSIZE_T) - 1_HSIZE_T, &
                INT(j1, HSIZE_T) - 1_HSIZE_T]

      ! For each hyperslab, we will just select 1 block from the dataspace.
      counter = [1_HSIZE_T, 1_HSIZE_T]

      ! We aren't doing any striding with the data--each process has a
      ! contiguous chunk of the data.
      stride = [1_HSIZE_T, 1_HSIZE_T]

      ! Size of the data for this process.
      block = [INT(nxP, HSIZE_T), INT(nyP, HSIZE_T)]

      ! Dimensions of the complete restart file.
      dimsT = [INT(nxW, HSIZE_T), INT(nyW, HSIZE_T)]

      ! Initialize the Fortran-HDF5 interface.
      CALL H5OPEN_F(ierr)

      ! Set up file access property list with parallel I/O access.
      info = MPI_INFO_NULL
      CALL H5PCREATE_F(H5P_FILE_ACCESS_F, plist_id, ierr)
      CALL H5PSET_FAPL_MPIO_F(plist_id, comm, info, ierr)

      ! Open the file for writing.
      CALL H5FCREATE_F(fName, H5F_ACC_TRUNC_F, file_id, ierr, &
                       access_prp=plist_id)

      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)

      ! Create the 'Header' dataspace.
      CALL H5SCREATE_F(H5S_NULL_F, dataspace_id, ierr)

      ! Create the 'Header' dataset.
      CALL H5DCREATE_F(file_id, 'Header', H5T_NATIVE_INTEGER, dataspace_id, &
                       dataset_id, ierr)

      ! Write header attributes.
      CALL WriteAttribute(dataset_id, H5T_NATIVE_INTEGER, 'nx', nxW)
      CALL WriteAttribute(dataset_id, H5T_NATIVE_INTEGER, 'ny', nyW)

      ! Close the 'Header' dataset
      CALL H5DCLOSE_F(dataset_id, ierr)

      ! Close the 'Header' dataspace
      CALL H5SCLOSE_F(dataspace_id, ierr)

      ! Create the 'FlowData' data group
      CALL H5GCREATE_F(file_id, 'FlowData', group_id, ierr)

      ! Write the data to the file using collective writeout.
      CALL WritePHDF5Dataset(nxP, nyP, group_id, 'Omega', w(1:nxW,1:nyW), &
                             offset, counter, stride, block, dimsT)
      CALL WritePHDF5Dataset(nxP, nyP, group_id, 'Psi', psi(1:nxW,1:nyW), &
                             offset, counter, stride, block, dimsT)
      CALL WritePHDF5Dataset(nxP, nyP, group_id, 'u', u(1:nxW,1:nyW), offset, &
                             counter, stride, block, dimsT)
      CALL WritePHDF5Dataset(nxP, nyP, group_id, 'v', v(1:nxW,1:nyW), offset, &
                             counter, stride, block, dimsT)

      ! Close the 'FlowData' data group.
      CALL H5GCLOSE_F(group_id, ierr)

      ! Close the file.
      CALL H5FCLOSE_F(file_id, ierr)

      ! Close the Fortran-HDF5 interface.
      CALL H5CLOSE_F(ierr)

      ! Write out the xml file needed to view the HDF5 data in Xdmf format.
      !CALL WriteRestartXML(0, nxT, nyT, nzT, fName, xlenT, ylenT, zlenT)
   END SUBROUTINE WriteRestartHDF5

   !> Routine to write a chunk of data with parallel HDF5 to a specified file
   !! space.
   !!
   !> @param[in] nxP Number of cells in the x direction for this process.
   !> @param[in] nyP Number of cells in the y direction for this process.
   !> @param[in] group_id Filespace ID for the data in the HDF5 file.
   !> @param[in] dName Name for the data.
   !> @param[in] data 2D array holding the data.
   !> @param[in] offset Offset for this data in the joined restart file.
   !> @param[in] counter Number of blocks from the hyperslab.
   !> @param[in] stride Strides to take in the HDF5 data in each direction.
   !> @param[in] block Size of the data on this process.
   !> @param[in] dimsT Dimensions of the total data in the restart file.
   SUBROUTINE WritePHDF5Dataset(nxP, nyP, group_id, dName, data, offset, &
                                counter, stride, block, dimsT)
      ! Required modules.
      USE HDF5
      ! Calling arguments.
      INTEGER(KIND=IWPF),INTENT(IN) :: nxP, nyP
      INTEGER(KIND=HID_T),INTENT(IN) :: group_id
      CHARACTER(LEN=*),INTENT(IN) :: dName
      REAL(KIND=RWPC),DIMENSION(1:nxP,1:nyP),INTENT(IN) :: data
      INTEGER(KIND=HSIZE_T),DIMENSION(2),INTENT(IN) :: offset, counter, &
                                                       stride, block, dimsT
      ! Local variables.
      ! HDF5 data space identifier.
      INTEGER(KIND=HID_T) :: filespace
      ! HDF5 data space identifier used for the memory space.
      INTEGER(KIND=HID_T) :: memspace
      ! Property list identifier.
      INTEGER(KIND=HID_T) :: plist_id
      ! HDF5 ID for the header dataset.
      INTEGER(KIND=HID_T) :: dataset_id
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! Create the data space for the data set.
      CALL H5SCREATE_SIMPLE_F(2, dimsT, filespace, ierr)
      CALL H5SCREATE_SIMPLE_F(2, block, memspace, ierr)

      ! Create the chunked data set.
      CALL H5PCREATE_F(H5P_DATASET_CREATE_F, plist_id, ierr)
      CALL H5PSET_CHUNK_F(plist_id, 2, block, ierr)
      CALL H5DCREATE_F(group_id, dName, H5T_NATIVE_DOUBLE, filespace, &
                       dataset_id, ierr, plist_id)
      CALL H5SCLOSE_F(filespace, ierr)

      ! Select the hyperslab in the file.
      CALL H5DGET_SPACE_F(dataset_id, filespace, ierr)
      CALL H5SSELECT_HYPERSLAB_F(filespace, H5S_SELECT_SET_F, offset, counter, &
                                 ierr, stride, block)

      ! Create the property list for collective dataset write.
      CALL H5PCREATE_F(H5P_DATASET_XFER_F, plist_id, ierr)
      CALL H5PSET_DXPL_MPIO_F(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)

      ! Write the data set collectively.
      CALL H5DWRITE_F(dataset_id, H5T_NATIVE_DOUBLE, data, dimsT, ierr, &
                      file_space_id=filespace, mem_space_id=memspace, &
                      xfer_prp=plist_id)

      ! Close the dataset.
      CALL H5DCLOSE_F(dataset_id, ierr)

      ! Close the property list.
      CALL H5PCLOSE_F(plist_id, ierr)

      ! Close the data spaces.
      CALL H5SCLOSE_F(filespace, ierr)
      CALL H5SCLOSE_F(memspace, ierr)
   END SUBROUTINE WritePHDF5Dataset

   !> Procedure to write an integer attribute to an HDF5 dataset.
   !!
   !> @param[in] dset_id HDF5 identifier for the data set.
   !> @param[in] type_id What type of HDF5 variable this is.
   !> @param[in] att_name Name for the attribute.
   !> @param[in] var The data to be written.
   SUBROUTINE WriteAttributeInteger(dset_id, type_id, att_name, var)
      ! Required modules.
      USE HDF5
      IMPLICIT NONE
      ! Calling arguments.
      INTEGER(KIND=HID_T),INTENT(IN) :: dset_id, type_id
      CHARACTER(LEN=*),INTENT(IN) :: att_name
      INTEGER(KIND=IWPF),INTENT(IN) :: var
      ! Local variables.
      ! Dataspace identifier for the attribute.
      INTEGER(KIND=HID_T) :: aspace_id
      ! Dataset identifier for the attribute.
      INTEGER(KIND=HID_T) :: attr_id
      ! Dimensions of the data to be written.
      INTEGER(HSIZE_T),DIMENSION(1),PARAMETER :: dims = [1]
      ! Error handling.
      INTEGER(KIND=IWPF) :: ierr

      ! Create the dataspace for the attribute.
      CALL H5SCREATE_F(H5S_SCALAR_F, aspace_id, ierr)
      !
      ! Create a dataset for the attribute.
      CALL H5ACREATE_F(dset_id, att_name, type_id, aspace_id, attr_id, ierr)
      !
      ! Write to the dataset.
      CALL H5AWRITE_F(attr_id, type_id, var, dims, ierr)
      !
      ! Close the dataspace for the dataset.
      CALL H5SCLOSE_F(aspace_id, ierr)
      !
      ! Close the dataset.
      CALL H5ACLOSE_F(attr_id, ierr)
   END SUBROUTINE WriteAttributeInteger

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
      10 FORMAT (A,A,I6.6,A,A)
   END SUBROUTINE FileNameWithNum

END MODULE IO_m

