! *****************************************************************************
! *                                                                           * 
! *  ECHO-QGP                                                                 * 
! *                                                                           * 
! *  Version: 1.0.00                                                          *
! *                                                                           *
! *  Date: May, 22 - 2015                                                     *
! *                                                                           *
! *  File: fromecho2d.f90                                                     *
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
! *  Contributors: Andrea Beraudo (beraudo@to.infn.it)                        *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************

      program fromecho2d
      implicit none
      integer nx, ny, nz ! number of points along x, y and z (i.e. eta) directions
      character(1) :: ax1_choice ! flag to choose between x,y or z axis
      character(1) :: ax2_choice ! flag to choose between x,y or z axis
      character(16) :: ax3_value ! value at which to cut the slice on the third axis
      character(10),parameter :: grid_data='grid.dat'
      character(30) :: input_file
      character(30) :: output_file
      integer filerror, ix, iy, iz, ixmiddle,iymiddle, izmiddle
      real :: xmiddle, ymiddle, zmiddle
      real, allocatable, dimension(:,:,:) :: dataset
      real, allocatable, dimension(:) :: xarray, yarray, zarray
      integer :: numarg

      numarg=command_argument_count()
      
      ! we check that the name of the input file and the chosen coordinates are given as arguments when the program is invoked
      if ((numarg .lt. 4) .or. (numarg .gt. 5)) then 
        call print_usage
        call exit(1)
      end if

      ! we get the name of the input file
      call get_command_argument(1,input_file)
      
      ! we get the name of the output file
      call get_command_argument(2,output_file)
      
      ! we choose which axis to draw ( x ,y, z) as first
      call get_command_argument(3,ax1_choice)

      ! we choose which axis to draw ( x ,y, z) as second
      call get_command_argument(4,ax2_choice)
      
      if(numarg .eq. 5) then
        ! we optionally choose at which value to cut the slice on the third axis
        call get_command_argument(5,ax3_value)
      end if

      !check that ax1_choice and ax2_choice are different
      if(ax1_choice .eq. ax2_choice) then
        write(*,*) 'Error, the axis defining the plane must be different... Exiting...'
        call exit(3)
      end if

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
      if(nx .gt. 1) then
        if(numarg .eq. 5) then
          read(ax3_value, '(f14.0)' )  xmiddle
        else
          xmiddle=0.
        end if
      else
         if ((ax1_choice .eq. 'x') .or. (ax2_choice .eq. 'x')) then 
              write(*,*) 'Sorry, you chose the x axis, but it has only 1 cell and this is an utility to extract 2D data...'
              write(*,*) 'Quitting...'
              call exit(11)
         end if
      end if
      
      if(ny .gt. 1) then
        if(numarg .eq. 5) then
          read(ax3_value, '(f14.0)' )  ymiddle
        else
          ymiddle=0.
        end if
      else
         if ((ax1_choice .eq. 'y') .or. (ax2_choice .eq. 'y')) then 
              write(*,*) 'Sorry, you chose the y axis, but it has only 1 cell and this is an utility to extract 2D data...'
              write(*,*) 'Quitting...'
              call exit(11)
         end if
      end if
      
      if(nz .gt. 1) then
        if(numarg .eq. 5) then
          read(ax3_value, '(f14.0)' )  zmiddle
        else
          zmiddle=0.
        end if
      else
         if ((ax1_choice .eq. 'z') .or. (ax2_choice .eq. 'z')) then 
              write(*,*) 'Sorry, you chose the z axis, but it has only 1 cell and this is an utility to extract 2D data...'
              write(*,*) 'Quitting...'
              call exit(11)
         end if
      end if
      
      ! now we write the data values     
      if((ax1_choice .eq. 'x') .and. (ax2_choice .eq. 'y')) then     
        call find_index(zarray,nz,zmiddle,izmiddle) 
        do ix=1,nx      
           do iy=1,ny  
             write(3,*) xarray(ix), yarray(iy), dataset(ix,iy,izmiddle)
           end do
           write(3,*)
        end do
      else if((ax1_choice .eq. 'y') .and. (ax2_choice .eq. 'x')) then     
        call find_index(zarray,nz,zmiddle,izmiddle)
        do iy=1,ny      
           do ix=1,nx  
             write(3,*) yarray(iy), xarray(ix), dataset(ix,iy,izmiddle)
           end do
           write(3,*)
        end do
      else if((ax1_choice .eq. 'x') .and. (ax2_choice .eq. 'z'))  then
        call find_index(yarray,ny,ymiddle,iymiddle)
        do ix=1,nx      
           do iz=1,nz  
             write(3,*) xarray(ix), zarray(iz), dataset(ix,iymiddle,iz)
           end do
           write(3,*)
        end do
      else if((ax1_choice .eq. 'z') .and. (ax2_choice .eq. 'x'))  then
        call find_index(yarray,ny,ymiddle,iymiddle)
        do iz=1,nz   
           do ix=1,nx
             write(3,*) zarray(iz), xarray(ix), dataset(ix,iymiddle,iz)
           end do
           write(3,*)
        end do
      else if((ax1_choice .eq. 'y') .and. (ax2_choice .eq. 'z')) then
        call find_index(xarray,nx,xmiddle,ixmiddle)
        do iy=1,ny     
          do iz=1,nz  
             write(3,*) yarray(iy), zarray(iz), dataset(ixmiddle,iy,iz)
          end do
          write(3,*)
        end do
      else if((ax1_choice .eq. 'z') .and. (ax2_choice .eq. 'y')) then
        call find_index(xarray,nx,xmiddle,ixmiddle)
        do iz=1,nz     
          do iy=1,ny  
             write(3,*) zarray(iz), yarray(iy), dataset(ixmiddle,iy,iz)
          end do
          write(3,*)
        end do
      else
        write(*,*) "Something has gone wrong when choosing along the two axis to plot."
        write(*,*)  "Please, check the invocation syntax.."
        call print_usage
      end if
      
           
      close(3)
      
      deallocate(xarray)
      deallocate(yarray)
      deallocate(zarray)
      deallocate(dataset)

      write(*,*) output_file," written." 
      
      contains
      
      subroutine print_usage
       write(*,*) 'Synopsis:'
       write(*,*) 'To produce 2d-plot in the x-y, x-z or y-z planes:'
       write(*,*) './fromecho input_file output_file x|y|z x|y|z (optional: third axis value, default is 0)'
       write(*,*) 'Usage example: ./fromecho prova_input prova_output x y'
       write(*,*) 'This extracts data from file prova_inputs and writes into file output_file the values on the x-y plane for z=0'
       write(*,*) 'Usage example: ./fromecho prova_input prova_output z x 5'
       write(*,*) 'This extracts data from file prova_inputs and writes into file output_file the values on the z-x plane for y=5'
       write(*,*) 'For the optional 5th argument, actually it is selected the grid value which is closest to the chosen value'
       write(*,*) 'i.e. no interpolations are made'
      end subroutine print_usage
     
      subroutine find_index(array_in, array_dim, target_value, index_found)
       implicit none
       real, allocatable, dimension(:) :: array_in
       integer :: array_dim
       real :: target_value
       integer :: index_found
       integer :: i
       
       if(array_dim .eq. 1) then
         index_found=1
         write(*,*) 'Third axis has one value only, selecting index=1 corresponding to value:',array_in(1)
         return
       end if

       do i=1,array_dim-1
          if((array_in(i) .le. target_value) .and. (array_in(i+1) .ge. target_value)) then
            if((abs(array_in(i) - target_value)) .lt. (abs(array_in(i+1) - target_value))) then
              index_found=i
            else
              index_found=i+1
            end if
            write(*,*) 'Index associated to target value:',target_value,' is:', index_found, ' corresponding to grid value ',&
                       &array_in(index_found)
            return
          end if
       end do

       write(*,*) 'Sorry, but your value for the the third axis lies outside the boundaries of the grid... Routine failed...'              
       call exit(5)
      end subroutine find_index
      
      end program
