! *****************************************************************************
! *                                                                           * 
! *  ECHO-QGP                                                                 * 
! *                                                                           * 
! *  Version: 1.0.00                                                          *
! *                                                                           *
! *  Date: May, 22 - 2015                                                     *
! *                                                                           *
! *  File: timev.f90                                                          *
! *                                                                           *
! *  License: GPL version 2.0 (Please, read the file LICENSE.TXT)             *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  This program is free software; you can redistribute it and/or            *
! *  modify it under the terms of the GNU General Public License              *
! *  as published by the Free Software Foundation; either version 2           *
! *  of the License, or (at your option) any later version.                   *
! *                                                                           *
! *  This program is distributed in the hope that it will be useful,          *
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
! *  GNU General Public License for more details.                             *
! *                                                                           *
! *  You should have received a copy of the GNU General Public License        *
! *  along with this program; if not, write to the Free Software              *
! *  Foundation, Inc., 51 Franklin Street, Fifth Floor,                       *
! *  Boston, MA  02110-1301, USA.                                             *
! *                                                                           *
! *  Authors: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)            *
! *                                                                           *
! *****************************************************************************

      program timev
      implicit none
      integer nx, ny, nz ! number of points along x, y and z (i.e. eta) directions
      integer i !just a counter
      integer start_index, end_index, x_index, y_index, z_index, file_index
      character(4) :: start_index_string, end_index_string, file_index_string
      character(10) :: x_index_string, y_index_string, z_index_string
      character(8),parameter :: grid_data='grid.dat'
      character(8),parameter :: time_data='time.dat'
      character(30) :: prefix
      character(30) :: output_file
      integer filerror, ix, iy, iz, ixcenter, iycenter, izcenter
      real, allocatable, dimension(:,:,:) :: dataset
      real, allocatable, dimension(:) :: xarray, yarray, zarray
      real time
      
009   format(I4.4)

      ! we check that the name of the input file and the chosen coordinates are given as arguments when the program is invoked
      if (command_argument_count() .ne. 7) then
      call print_usage
      call exit(1)
      end if
      
      ! we get the name of the variable prefix
      call get_command_argument(1,prefix)

      ! we choose the initial time index
      call get_command_argument(2,start_index_string)
      read(start_index_string,*) start_index
      
      ! we choose the final time index
      call get_command_argument(3,end_index_string)
      read(end_index_string,*) end_index

      ! we choose the x index of the variable
      call get_command_argument(4,x_index_string)
      read(x_index_string,*) x_index

      ! we choose the y index of the variable
      call get_command_argument(5,y_index_string)
      read(y_index_string,*) y_index

      ! we choose the z index of the variable
      call get_command_argument(6,z_index_string)
      read(z_index_string,*) z_index

      ! we get the name of the output file
      call get_command_argument(7,output_file)
      
      ! now we try to open the file containing the values of the cell coordinates of the grid
      open(2,file=grid_data,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
      ! we try to open the the grid data file assuming that we are into the tools directory
          write(*,*) grid_data, " cannot be opened, so, I'm forced to quit!"
          call exit(2)
      end if
      
      ! now we try to open the file containing the values of the time values corresponding to the output files
      open(11,file=time_data,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
          write(*,*) time_data, " cannot be opened neither in the current or in the ../out or in the ../../out directory."
          write(*,*) "So, I'm forced to quit!"
          call exit(11)
      end if
       
      ! here we read the number of cells along the various directions
      read(2,*) nx, ny, nz
           
      ! we allocate memory for the array containing the values of the cells along x
      allocate(xarray(nx))
      
      ! we allocate memory for the array containing the values of the cells along y
      allocate(yarray(ny))
      
      ! we allocate memory for the array containing the values of the cells along z
      allocate(zarray(nz))

      ! now we fill the two arrays with coordinates values
      read(2,*) xarray
            
      read(2,*) yarray
      
      read(2,*) zarray
      
      ! we've finished with the grid file, so we close it
      close(2)
       
      ! now we try to open the output file
      open(10,file=output_file,status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
       write(*,*) "I cannot open the output file, so I'm forced to quit!"
       call exit(3)
      end if
      
      ! now we allocate memory for the data
      allocate(dataset(nx,ny,nz))

      !now we move into the file containing time values up to the first entry

      if(start_index .gt. 1) then
         do i=1,start_index-1
            read(11,*)
         end do
      end if
      
      do file_index=start_index,end_index
        
       write(file_index_string,009) file_index
       
       ! we try to open the input file
       open(3,file=trim(adjustl(prefix))//file_index_string//'.dat',status='OLD',iostat=filerror)
       if (filerror .ne. 0) then
       write(*,*) trim(adjustl(prefix))//file_index_string//'.dat', " cannot be opened..."
       write(*,*) "Are you sure that it exists into the current directory?"
       write(*,*) "However, I'm forced to quit!"
       call exit(1)
       end if
      
      
       ! now we read the data values
       do iz=1,nz
         do iy=1,ny
           do ix=1,nx
             read(3,*) dataset(ix,iy,iz)
           end do
         end do
       end do
      
       close(3)

       read(11,"(5x,f12.9)") time
      
       write(10,"(e14.9,1x,e14.9)") time, dataset(x_index,y_index,z_index)
         
      end do
      close(10)
      close(11)
      deallocate(xarray)
      deallocate(yarray)
      deallocate(zarray)
      deallocate(dataset)

      contains
      
      subroutine print_usage
       write(*,*) 'Synopsis: ./timev prefix_of_the_variable starting_index final_index x_index y_index z_index output_filename '
       write(*,*) 'Example:'
       write(*,*) './timev T 1 137 51 40 37 pippo'
       write(*,*) 'It prints into the file named "pippo" two columns of values: the first one contains the time, the second one'
       write(*,*) 'the values of the variable (usually T stands for temperature) at cell of indexes x 51, y 40 and z 37'
       write(*,*) 'stored in the files from T0001.dat to T0137.dat'
      end subroutine print_usage
      
      end program
