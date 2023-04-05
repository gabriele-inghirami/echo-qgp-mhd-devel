! *****************************************************************************
! *                                                                           * 
! *  ECHO-QGP                                                                 * 
! *                                                                           * 
! *  Version: 1.0.00                                                          *
! *                                                                           *
! *  Date: May, 22 - 2015                                                     *
! *                                                                           *
! *  File: fromecho.f90                                                       *
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
! *  Acknowledgments: Arturo de Pace, Andrea Beraudo                          *
! *                                                                           *
! *****************************************************************************

      program fromecho
      implicit none
      integer nx, ny, nz ! number of points along x, y and z (i.e. eta) directions
      character(1) :: ax_choice ! flag to choose between x,y or z axis
      character(5) :: x_start_index_string, x_end_index_string, y_start_index_string, y_end_index_string
      character(5) :: z_start_index_string, z_end_index_string 
      character(10),parameter :: grid_data='grid.dat'
      character(30) :: input_file
      character(30) :: output_file
      integer filerror, ix, iy, iz, ixstart,ixend, iystart, iyend, izstart, izend
      real, allocatable, dimension(:,:,:) :: dataset
      real, allocatable, dimension(:) :: xarray, yarray, zarray
      integer :: numarg

      numarg=command_argument_count()
      
      ! we check that the name of the input file and the chosen coordinates are given as arguments when the program is invoked
      if (numarg .lt. 3) then 
        call print_usage
        call exit(1)
      else if (numarg .eq. 3) then
      else if (numarg .ne. 9) then
        call print_usage
        call exit(1)
      end if

      ! we get the name of the input file
      call get_command_argument(1,input_file)
      
      ! we get the name of the output file
      call get_command_argument(2,output_file)
      
      ! we choose which axis to draw ( x ,y, z)
      call get_command_argument(3,ax_choice)
      
      if(numarg .gt. 3) then
        ! we retrieve the starting x index, if provided
        call get_command_argument(4,x_start_index_string)
        read(x_start_index_string,'(i5)') ixstart  

        ! we retrieve the ending x index, if provided
        call get_command_argument(5,x_end_index_string)
        read(x_end_index_string,'(i5)') ixend  

        ! we retrieve the starting y index, if provided
        call get_command_argument(6,y_start_index_string)
        read(y_start_index_string,'(i5)') iystart  

        ! we retrieve the ending y index, if provided
        call get_command_argument(7,y_end_index_string)
        read(y_end_index_string,'(i5)') iyend  

        ! we retrieve the starting z index, if provided
        call get_command_argument(8,z_start_index_string)
        read(z_start_index_string,'(i5)') izstart  

        ! we retrieve the ending z index, if provided
        call get_command_argument(9,z_end_index_string)
        read(z_end_index_string,'(i5)') izend  

      endif

      ! we try to open the input file
      open(1,file=input_file,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
      write(*,*) input_file, " cannot be opened, are you sure that it exists into the current directory?"
      write(*,*) "However, I'm forced to quit!"
      call exit(1)
      end if
      
      ! now we try to open the file containing the values of the cell coordinates of the grid
      open(2,file=grid_data,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
          write(*,*) 'Sorry, it was not possibile to open the file grid.dat, so I quit...'
         call exit(2)
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
      
      ! now we allocate memory for the data
      allocate(dataset(nx,ny,nz))
      
      ! now we read the data values
      do iz=1,nz
        do iy=1,ny
          do ix=1,nx
             read(1,*) dataset(ix,iy,iz)
          end do
        end do
      end do
      
      close(1)
      
      ! now we try to open the output file
      open(3,file=output_file,status='REPLACE', iostat=filerror)
      if (filerror .ne. 0) then
      write(*,*) "I cannot open the output file, so I'm forced to quit!"
      call exit(3)
      end if
       
      !we set the indexes for extracting data
      if(nx .eq. 1) then
          ixstart=1
          ixend=1
      else
          if(numarg .ne. 9) then
            ixstart=nx/2+1
            ixend=nx
          else
            if(( ixstart .gt. ixend) .or. (ixstart .gt. nx) .or. (ixend .gt. nx)) then
              write(*,*) 'Sorry, the x indexes you chose are not allowed, please, check them...'
              call exit(11)
            end if
          end if
      end if
      
      if(ny .eq. 1) then
          iystart=1
          iyend=1
      else
          if(numarg .ne. 9) then
            iystart=ny/2+1
            iyend=ny
          else
            if(( iystart .gt. iyend) .or. (iystart .gt. ny) .or. (iyend .gt. ny)) then
              write(*,*) 'Sorry, the y indexes you chose are not allowed, please, check them...'
              call exit(12)
            end if
          end if
      end if
      
      if(nz .eq. 1) then
          izstart=1
          izend=1
      else
          if(numarg .ne. 9) then
            izstart=nz/2+1
            izend=nz
          else
            if(( izstart .gt. izend) .or. (izstart .gt. nz) .or. (izend .gt. nz)) then
              write(*,*) 'Sorry, the z indexes you chose are not allowed, please, check them...'
              call exit(13)
            end if
          end if
      end if

      ! now we write the data values     
      if(ax_choice .eq. 'x') then     
       if(((iystart .ne. iyend) .or. (izstart .ne. izend)) .and. (numarg .eq. 9)) then
          write(*,*) "Please, choose y index start = y index end and z index start = z index end"
          call exit(14)
       end if
       do ix=ixstart,ixend
         write(3,*) xarray(ix), dataset(ix,iystart,izstart)
       end do
      else if(ax_choice .eq. 'y')  then
       if(((ixstart .ne. ixend) .or. (izstart .ne. izend)) .and. (numarg .eq. 9)) then
          write(*,*) "Please, choose x index start = x index end and z index start = z index end"
          call exit(15)
       end if
       do iy=iystart,iyend
         write(3,*) yarray(iy), dataset(ixstart,iy,izstart)
       end do
      else if(ax_choice .eq. 'z') then
       if(((ixstart .ne. ixend) .or. (iystart .ne. iyend)) .and. (numarg .eq. 9)) then
          write(*,*) "Please, choose x index start = x index end and y index start = y index end"
          call exit(16)
       end if
       do iz=izstart,izend
         write(3,*) zarray(iz), dataset(ixstart,iystart,iz)
       end do
      else
      write(*,*) "Something has gone wrong when choosing along which axis to plot."
      write(*,*)  "Please, check the invocation syntax.."
      call print_usage
      end if
      
           
      close(3)
      
      deallocate(xarray)
      deallocate(yarray)
      deallocate(zarray)
      deallocate(dataset)
      
      contains
      
      subroutine print_usage
       write(*,*) 'Synopsis:'
       write(*,*) 'To extract data from the center of the grid to the right border along x, y or z direction:'
       write(*,*) './fromecho input_file output_file x|y|z'
       write(*,*) 'To extract data using a specific range along one direction and fixing the other ones at a point:'
       write(*,*) './fromecho input_file output_file x|y|z x-indx-start x-indx-end y-indx-start y-indx-end z-indx-start z-indx-end'
       write(*,*) 'Please, note that starting and ending indexes of the fixed directions must be equal'
       write(*,*) 'Usage example: ./fromecho prova_input prova_output y 33 33 2 150 44 44'
       write(*,*) 'This extracts data from file prova_inputs and write the values along y direction from point with index 2 to'
       write(*,*) 'point with index 150, fixing x at index 33 and z at index 44'
      end subroutine print_usage
      
      end program
