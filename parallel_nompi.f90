! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  File: parallel_nompi.f90                                                 *
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
!-- Module for scalar runs (no calls to MPI routines)

  logical,parameter :: prl=.false.,pe0=.true.

  integer,parameter :: ipe=0,npe=1,ipe_prev=0,ipe_next=0, ipec=0

  integer :: mpi_integer,mpi_realtype,icomm,istat,ierr,mpi_max,mpi_min,mpi_sum,MPI_COMM_WORLD,MPI_LAND,MPI_Logical,MPI_CHARACTER

  logical :: check_crash

  integer :: nx,ny,nz,mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2

contains

! *****************************************************************************

subroutine parallel_start()
  return
end subroutine parallel_start

! *****************************************************************************

subroutine parallel_grid_setup(nx,ny,nz,mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2)

  integer,intent( in) :: nx,ny,nz
  integer,intent(out) :: mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2

  mx=nx
  my=ny
  mz=nz

  ix0=0; ix1=1; ix2=nx
  iy0=0; iy1=1; iy2=ny
  iz0=0; iz1=1; iz2=nz

end subroutine parallel_grid_setup

! *****************************************************************************

subroutine parallel_end

end subroutine parallel_end

! *****************************************************************************

end module parallel

! *****************************************************************************
