! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  File: parallel_mpi.f90                                                   *
! *                                                                           *
! *  License: GPL version 2.0 (Please, read the file LICENSE.TXT)             *
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
! *  Authors: Luca Del Zanna (luca.delzanna@unifi.it)                         *
! *                                                                           *
! *  Contributors: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)       *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************

module parallel
!-- Interface for MPI routines

  use mpi 
  implicit none

  !include 'mpif.h'

  logical :: prl,pe0
  integer :: ipe,npe,ipe_prev,ipe_next,ipec
  integer :: mpi_realtype,icomm=mpi_comm_world,istat(mpi_status_size),ierr
  logical :: check_crash
  integer :: nx,ny,nz,mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2

contains

! *****************************************************************************
subroutine parallel_start()
  call mpi_init(ierr)

!  line commented because from May 2017 we switched to 8 bytes precision only
!  if(kind(1.0)==4) mpi_realtype=mpi_real
  mpi_realtype=mpi_double_precision

  call mpi_comm_rank(icomm,ipe,ierr)
  call mpi_comm_size(icomm,npe,ierr)
  pe0=ipe==0
  prl=npe>1

end subroutine parallel_start

! *****************************************************************************

subroutine parallel_grid_setup(nx,ny,nz,mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2)

  integer,intent( in) :: nx,ny,nz
  integer,intent(out) :: mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2


  !check if the number of processor is greater than the number of cells divided by 3
  if(npe .gt. nx/3) then
     write(*,*) 'Sorry, but you cannot run a parallel simulation with more processors than the numb. of cells along x divided by 3'
     write(*,*) 'Number of processors used:', npe
     write(*,*) 'Number of cells along x:', nx
     call exit(1)
  end if
    
  ipe_prev=modulo(ipe-1,npe)
  ipe_next=modulo(ipe+1,npe)

  ipec=mod(nx,npe)  

  if (ipe .lt. ipec) then
     mx=nx/npe+1
  else
     mx=nx/npe
  end if
  my=ny
  mz=nz
 
  if (ipe .lt. ipec) then
     ix0=ipe*mx; ix1=ix0+1; ix2=ix0+mx
  else
     ix0=ipec*(mx+1)+(ipe-ipec)*mx; ix1=ix0+1; ix2=ix0+mx
  end if
  iy0=0     ; iy1=1    ; iy2=ny
  iz0=0     ; iz1=1    ; iz2=nz

end subroutine parallel_grid_setup

! *****************************************************************************

subroutine parallel_end

  call mpi_finalize(ierr)

end subroutine parallel_end

! *****************************************************************************

end module parallel

! *****************************************************************************
