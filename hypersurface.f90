! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2016 The ECHO-QGP team                                * 
! *                                                                           *
! *  File: hypersurface.f90                                                   *
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
! *  Authors: Valentina Rolando (rolando@fe.infn.it)                          *         
! *                                                                           *
! *  Contributors: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)       *
! *                                                                           *
! *  Acknowledgments: Giuseppe Pagliara, Raffaele Tripiccione                 *
! *                                                                           *
! *****************************************************************************

module hypersurface
  use eos
  use viscous
  use common
  use parallel
  
  procedure(), pointer ::  find_hypersurface
  
  !  This array must be the field on wich we want to check the hypersurface
  !  e.g. the temperature. The first index is for the time-direction: 0 is the previous step in tau, 
  ! 1 is the next. Every time the check is done
  !  field(1, i, j, k)=  field(2, i, j, k)
  !  and 
  !  field(2,i,j,k) is filled with a new time-step
  real(8), allocatable, dimension (:,:,:,:) :: field
  ! time coordinate for the 2 slices  e.g.:
  ! time(1)=tau
  ! time(2)=tau + d tau
  real(8), dimension(2) :: time_arr
  ! ! Must be a backup of the variables at the previous time step considered
  ! ! (can be void) 
  ! !  primitive must correspond to time(2)
  ! !  v_backup  must correspond to time(1)
  real(8), allocatable, dimension (:,:,:,:) :: v_backup, primitives
  real(8), allocatable, dimension (:,:,:,:) :: der_backup

  integer, parameter :: tau_dir=1, x_dir=2, y_dir=3, z_dir=4

  real(8), dimension(:,:), allocatable :: hy_ide, hy_vis, hy_deriv, hy_mhd

  integer :: nit_mhd !3 or 4 depending if the carge density in the comoving frame is printed or not (print_rho_comov true or false)
  ! always
  ! real(8), intent(out), dimension(2*10*nx*ny*nz, 1:13), allocatable :: hy_ide
  ! only if viscous==1
  ! real(8), intent(out), dimension(2*10*nx*ny*nz, 1:7), allocatable :: hy_vis
  ! in hy_ide the indexes are:
  ! 1-tau, 2-x, 3-y, 4-z, 5-rho, 6-vx, 7-vy, 8-vz, 9-prex, 10-dV_t, 11-dV_x, 12-dV_y, 13-dV_z
  ! = 13 entries 
  ! here hy_vis is just a dummy number as 0.0


! **************************************
contains

! **************************************
subroutine check_processor

  if(.not. pe0) then
    write(*,*) "Sorry, but the subroutine allocate_hypersurface_arrays must be called only from processor 0, while I am proc"
    write(*,*) ipe
    call exit(1)
  end if

end subroutine check_processor

! **************************************

subroutine allocate_hypersurface_arrays
  implicit none
  integer allocate_result

  !we check that only processor 0 is calling this subroutine
  call check_processor()

  allocate(field(1:2,1:nx,1:ny,1:nz),v_backup(1:nx,1:ny,1:nz,1:nv),primitives(1:nx,1:ny,1:nz,1:nv), STAT=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate v_backup, primitives or field arrays into allocate_hypersurface_arrays,"
    write(*,*) "contained into hypersurface.f08"
    write(*,*) "all result is:", allocate_result

    call exit(1)
  end if
  field=0.
  v_backup=0.
  primitives=0.
  
  allocate(hy_ide(4*nx*ny*nz, 1:13), STAT = allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate , primitives or field arrays into allocate_hypersurface_arrays,"
    write(*,*) "contained into hypersurface.f08"
    call exit(1)
  end if
  hy_ide=0.

  if (derivatives_out) then
    ALLOCATE (hy_deriv(4*nx*ny*nz, 1:nderivatives), STAT = allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Error, I can't allocate hy_der into allocate_hypersurface_arrays, contained into hypersurface.f08"
      call exit(1)
    end if
    hy_deriv=0.
    ALLOCATE (der_backup(1:nx,1:ny,1:nz, 1:nderivatives), STAT = allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Error, I can't allocate der_backup into allocate_hypersurface_arrays, contained into hypersurface.f08"
      call exit(1)
    end if
    der_backup=0.
  end if
 
  if (viscosity) then
    ALLOCATE (hy_vis(4*nx*ny*nz, 1:7), STAT = allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Error, I can't allocate hy_vis into allocate_hypersurface_arrays, contained into hypersurface.f08"
      call exit(1)
    end if
    hy_vis=0.
    find_hypersurface => find_hypersurface_visco
  else
    ALLOCATE (hy_vis(1,1), STAT = allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Error, I can't allocate hcy_vis into allocate_hypersurface_arrays, contained into hypersurface.f08"
      call exit(1)
    end if
    hy_vis=0.
    find_hypersurface => find_hypersurface_ideal
  endif

  if (mhd) then
    if(print_rho_comov) then
      nit_mhd=4
    else
      nit_mhd=3
    end if
    ALLOCATE (hy_mhd(4*nx*ny*nz, 1:nit_mhd), STAT = allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Error, I can't allocate hy_mhd into allocate_hypersurface_arrays, contained into hypersurface.f08"
      call exit(1)
    end if
    hy_mhd=0.
    find_hypersurface => find_hypersurface_mhd
  end if

end subroutine allocate_hypersurface_arrays

! **************************************

subroutine find_hypersurface_ideal()
  implicit none
  !  subroutine variables
  real(8) sign_
  integer ix, iy, iz,iv
  real(8), dimension(0:3):: xfo, vfo, dVfo
  
  real(8), dimension (1:13) :: info_ide
  real(8), dimension (1:nderivatives) :: info_der !used for derivatives computation
  real(8)  info_vis

  integer allocate_error

  integer i_hysu
  i_hysu=0
  !we check that only processor 0 is calling this subroutine
  call check_processor()

  do ix=1,nx
    do iy=1,ny
      do iz=1,nz
        ! check on tau (t)---------------------------------------
        if ((field(2,ix,iy,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
          sign_=sign(1.0d0,(field(1,ix,iy,iz)-field(2,ix,iy,iz)))
          call set_fo_ide(1, ix,iy,iz, 2, ix, iy, iz, info_ide, info_vis, info_der, tau_dir)
          call save_hypersurface_ideal(i_hysu, info_ide, info_vis,sign_, info_der)
        endif
        ! check on x -------------------------------------------
        if (ix<nx) then
          if ((field(1,ix+1,iy,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix+1,iy,iz)))
            call set_fo_ide(1, ix,iy,iz, 1, ix+1, iy, iz, info_ide, info_vis, info_der, x_dir)
            call save_hypersurface_ideal(i_hysu, info_ide, info_vis,sign_, info_der)
          endif
        endif
        ! check on y -------------------------------------------
        if (iy<ny) then
          if ((field(1,ix,iy+1,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix,iy+1,iz)))
            call set_fo_ide(1, ix,iy,iz, 1, ix, iy+1, iz, info_ide, info_vis, info_der, y_dir)
            call save_hypersurface_ideal(i_hysu, info_ide, info_vis,sign_, info_der)
          endif
        endif
        ! check on eta (z) -------------------------------------
        if (iz<nz) then
          if ((field(1,ix,iy,iz+1)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix,iy,iz+1)))
            call set_fo_ide(1, ix,iy,iz, 1, ix, iy, iz+1, info_ide, info_vis, info_der, z_dir)
            call save_hypersurface_ideal(i_hysu, info_ide, info_vis, sign_, info_der)
            endif
        endif
      end do
    end do
  end do
  
  call print_hypersurface_ideal(i_hysu)
  
end subroutine find_hypersurface_ideal
! *****************************************************************************
subroutine print_hypersurface_ideal(ind_h) 
!   use eos
  implicit none 
  
  integer, intent(in) :: ind_h

  integer fileErr, i, j, erralert
  real(8) temp, rho, eng, prex, entropy
  
  open (21,file=prefix_dir//'hypersurface.txt',form='formatted',status='old', access='append', iostat=fileErr)
    if (fileErr .ne. 0) then 
      write (*,*) " I cannot open the file ", prefix_dir//'hypersurface.txt' 
      write (*,*) " I have to quit "
      call exit(1)
    endif

  if(derivatives_out) then
   open (23,file=prefix_dir//'hypersurf_deriv.txt',form='formatted',status='old', access='append', iostat=fileErr)
    if (fileErr .ne. 0) then
      write (*,*) " I cannot open the file ", prefix_dir//'hypersurf_deriv.txt'
      write (*,*) " I have to quit "
      call exit(1)
    endif
  endif     

12 format (i6)
21 format (14(ES14.7, 2x))
22 format (3(ES14.7, 2x))
23 format (20(ES14.7, 2x))

    write(21, 12) ind_h 
    do i=1, ind_h
      rho=hy_ide(i,5)
      prex=hy_ide(i,9)

      write(21, 21) (hy_ide(i, j), j=1,13)
      if(derivatives_out) then
        write(23,23) (hy_deriv(i,j), j=1,nderivatives)
      end if
    end do
     close(21)
     if(derivatives_out) then
       close(23)
     end if
  return
end subroutine print_hypersurface_ideal

! *****************************************************************************

subroutine save_hypersurface_ideal(ind_hys, ide_fo, vis_fo, si, der_fo)
  use eos
  use viscous
  implicit none
  integer, intent(inout):: ind_hys
  real(8), intent(in), dimension (1:13) :: ide_fo
  real(8), intent(in), dimension (1:nderivatives) :: der_fo !derivatives for derivatives computation
  real(8), intent(in) :: vis_fo
  real(8), intent(in) :: si

  real(8) temperature, zz, dens, press, energy, entropy, vx, vy, vz, time, glf,gcov3
  real(8) dutdt, dutdx, dutdy, dutdz, duxdt, duxdx, duxdy, duxdz, duydt, duydx, duydy, duydz, duzdt, duzdx, duzdy, duzdz, dtedt,&
      &dtedx, dtedy, dtedz
  integer errors

  integer i

  errors=0

  ind_hys=ind_hys+1
  do i=1, 9
    hy_ide(ind_hys, i)=ide_fo(i)
  end do
  do i=10, 13
    hy_ide(ind_hys, i)=si*ide_fo(i)
  end do
! !  just a dummy variable
  hy_vis=0.0

  if(derivatives_out) then
    dutdt=der_fo(dtt)
    dutdx=der_fo(dtx)
    dutdy=der_fo(dty)
    dutdz=der_fo(dtz)
    duxdt=der_fo(dxt)
    duxdx=der_fo(dxx)
    duxdy=der_fo(dxy)
    duxdz=der_fo(dxz)
    duydt=der_fo(dyt)
    duydx=der_fo(dyx)
    duydy=der_fo(dyy)
    duydz=der_fo(dyz)
    duzdt=der_fo(dzt)
    duzdx=der_fo(dzx)
    duzdy=der_fo(dzy)
    duzdz=der_fo(dzz)
    dtedt=der_fo(dtet)
    dtedx=der_fo(dtex)
    dtedy=der_fo(dtey)
    dtedz=der_fo(dtez)
    dens=ide_fo(5)
    press=ide_fo(9)
    temperature=freeze_value
    time=ide_fo(1)
    zz=ide_fo(4)
    vx=ide_fo(6)
    vy=ide_fo(7)
    vz=ide_fo(8)
    if(coordinates .eq. MINKOWSKI) then
      gcov3=1.
    else
      gcov3=time**2.
    end if
    glf=1./sqrt(1.-vx*vx-vy*vy-vz*vz*gcov3)
    do i=1,nderivatives
      hy_deriv(ind_hys,i)=der_fo(i)
    end do
  end if

  return
end subroutine save_hypersurface_ideal

! *****************************************************************************

subroutine set_fo_ide(it, ix,iy,iz, iit, iix, iiy, iiz, fo_ide, fo_vis, fo_der, dir_flag)
use common, only : x,y,z, g_cov, gp,  gm, &
  & krh, kvx, kvy, kvz, kpr, d_val
implicit none
  integer, intent(in) :: it, ix, iy, iz, iit, iix, iiy, iiz

  integer :: idx !just a counter
  
  real(8), intent(out), dimension (1:13) :: fo_ide
  real(8), intent(out), dimension (1:nderivatives) :: fo_der
! here everything is ideal, fo_vis is just a dummy variable
  real(8), intent(out) :: fo_vis 

  ! subroutine variables
  real(8) weight, betasq
  real(8), dimension(0:3) :: vfo, xfo
  real(8) dx, dy, dz, dt
  real(8) ddx, ddy, ddz, ddt

  integer, intent(in) :: dir_flag


  fo_vis=0.0
  
  weight= (field(it,ix,iy,iz)-freeze_value)/(field(it,ix,iy,iz)-field(iit,iix,iiy,iiz))

  xfo(0)=time_arr(it)+weight*(time_arr(iit)-time_arr(it))
  xfo(1)=x(ix)+weight*(x(iix)-x(ix))
  xfo(2)=y(iy)+weight*(y(iiy)-y(iy))
  xfo(3)=z(iz)+weight*(z(iiz)-z(iz))

  if(dir_flag .eq. tau_dir) then
    vfo(1)=v_backup(ix,iy,iz,kvx)+weight*(primitives(iix,iiy,iiz,kvx)-v_backup(ix,iy,iz,kvx))
    vfo(2)=v_backup(ix,iy,iz,kvy)+weight*(primitives(iix,iiy,iiz,kvy)-v_backup(ix,iy,iz,kvy))
    vfo(3)=v_backup(ix,iy,iz,kvz)+weight*(primitives(iix,iiy,iiz,kvz)-v_backup(ix,iy,iz,kvz))
  else
    vfo(1)=v_backup(ix,iy,iz,kvx)+weight*(v_backup(iix,iiy,iiz,kvx)-v_backup(ix,iy,iz,kvx))
    vfo(2)=v_backup(ix,iy,iz,kvy)+weight*(v_backup(iix,iiy,iiz,kvy)-v_backup(ix,iy,iz,kvy))
    vfo(3)=v_backup(ix,iy,iz,kvz)+weight*(v_backup(iix,iiy,iiz,kvz)-v_backup(ix,iy,iz,kvz))
  end if

  betasq=vfo(1)*vfo(1)*(g_cov(1))+vfo(2)*vfo(2)*(g_cov(2))+vfo(3)*vfo(3)*(g_cov(3))
  if (betasq>1.0) then 
    call renormalize(vfo, xfo)
    print *, "you have to fix this"
    call exit(1)
  endif

  ! tau  or t 
  fo_ide(1)=xfo(0)
  ! x
  fo_ide(2)=xfo(1)
  ! y
  fo_ide(3)=xfo(2)
  ! eta or z
  fo_ide(4)=xfo(3)

  ! rho
  if(dir_flag .eq. tau_dir) then
    fo_ide(5)=v_backup(ix,iy,iz,krh)+weight*(primitives(iix,iiy,iiz,krh)-v_backup(ix,iy,iz,krh))
  else
    fo_ide(5)=v_backup(ix,iy,iz,krh)+weight*(v_backup(iix,iiy,iiz,krh)-v_backup(ix,iy,iz,krh))
  end if

  ! v^x
  fo_ide(6)=vfo(1)
  ! v^y
  fo_ide(7)=vfo(2)
  ! v^eta or v^z
  fo_ide(8)=vfo(3)

  ! prex
  if(dir_flag .eq. tau_dir) then
    fo_ide(9)=v_backup(ix,iy,iz,kpr)+weight*(primitives(iix,iiy,iiz,kpr)-v_backup(ix,iy,iz,kpr))
  else
    fo_ide(9)=v_backup(ix,iy,iz,kpr)+weight*(v_backup(iix,iiy,iiz,kpr)-v_backup(ix,iy,iz,kpr))
  end if

  dt=abs(time_arr(2)-time_arr(1))
  dx= d_val(1)
  dy= d_val(2)
  dz= d_val(3)

! dV perp t
  fo_ide(10)=real((iit-it),8)*dx*dy*dz*xfo(0)
! dV perp x
  fo_ide(11)=real((iix-ix),8)*dt*dy*dz*xfo(0)
! dV perp y
  fo_ide(12)=real((iiy-iy),8)*dx*dt*dz*xfo(0)
! dV perp z
  fo_ide(13)=real((iiz-iz),8)*dx*dy*dt*xfo(0)

 ! HERE STARTS THE SECTION ADDED FOR STORING DERIVATIVES NEEDED FOR COMPUTING TERMAL VORTICITY
  if(derivatives_out) then
    if(dir_flag .eq. tau_dir) then
      do idx=1,nderivatives
         fo_der(idx)=der_backup(ix,iy,iz,idx)+weight*(derivatives_all(iix,iiy,iiz,idx)-der_backup(ix,iy,iz,idx))
      end do
    else
      do idx=1,nderivatives
         fo_der(idx)=der_backup(ix,iy,iz,idx)+weight*(der_backup(iix,iiy,iiz,idx)-der_backup(ix,iy,iz,idx))
      end do
    end if
  end if


end subroutine set_fo_ide
! *****************************************************************************
! *****************************************************************************


! *****************************************************************************
subroutine find_hypersurface_mhd()
implicit none
!  subroutine variables
real(8) sign_
integer ix, iy, iz, iv
real(8), dimension(0:3):: xfo, vfo, dVfo

real(8), dimension (1:13) :: info_ide
real(8), dimension (1:4) ::  info_mhd !!Here we allocate it differently
real(8), dimension (1:nderivatives) :: info_der !used for derivatives computation
  
integer i_hysu

  i_hysu=0
  !we check that only processor 0 is calling this subroutine
  call check_processor()
  do ix=1, nx
    do iy=1, ny
      do iz=1, nz
	! check on tau (t)---------------------------------------
	if ((field(2,ix,iy,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) > 0.0) then ! FO 
	  sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(2,ix,iy,iz)))
	  call set_fo_mhd(1,ix,iy,iz, 2,ix,iy,iz, info_ide, info_mhd, info_der, tau_dir)
	  call save_hypersurface_mhd(i_hysu, info_ide, info_mhd, sign_, info_der)
	endif
	! check on x -------------------------------------------
	if (ix<nx) then
	  if ((field(1,ix+1,iy,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix+1,iy,iz)))
	    call set_fo_mhd(1,ix,iy,iz, 1,ix+1,iy,iz, info_ide, info_mhd, info_der, x_dir)
	    call save_hypersurface_mhd(i_hysu, info_ide, info_mhd, sign_, info_der)
	  endif
	endif
	! check on y -------------------------------------------
	if (iy<ny) then
	  if ((field(1,ix,iy+1,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix,iy+1,iz)))
	    call set_fo_mhd(1,ix,iy,iz, 1,ix,iy+1,iz, info_ide, info_mhd, info_der, y_dir)
	    call save_hypersurface_mhd(i_hysu, info_ide, info_mhd, sign_, info_der)
	  endif
	endif
	! check on eta (z) -------------------------------------
	if (iz<nz) then
	  if ((field(1,ix,iy,iz+1)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix,iy,iz+1)))
	    call set_fo_mhd(1,ix,iy,iz, 1,ix,iy,iz+1, info_ide, info_mhd, info_der, z_dir)
	    call save_hypersurface_mhd(i_hysu, info_ide, info_mhd, sign_, info_der)
	  endif
	endif
      end do
    end do
  end do
  
  call print_hypersurface_mhd(i_hysu)
  return
end subroutine find_hypersurface_mhd
! *****************************************************************************
subroutine print_hypersurface_mhd(ind_h)
  use common, only: kbx, kby, kbz, nv
  use eos
!   use viscous
  implicit none 

  integer, intent(in) :: ind_h

  integer fileErr, i, j

  open (21,file=prefix_dir//'hypersurface.txt',form='formatted',status='old', access='append', iostat=fileErr)
    if (fileErr .ne. 0) then 
      write (*,*) " I cannot open the file ", prefix_dir//'hypersurface.txt' 
      write (*,*) " I have to quit "
      call exit(1)
    endif
  if(derivatives_out) then
   open (23,file=prefix_dir//'hypersurf_deriv.txt',form='formatted',status='old', access='append', iostat=fileErr)
    if (fileErr .ne. 0) then
      write (*,*) " I cannot open the file ", prefix_dir//'hypersurf_deriv.txt'
      write (*,*) " I have to quit "
      call exit(1)
    endif
  endif 

12 format (i6)
20 format (16(ES14.7, 2x))
21 format (17(ES14.7, 2x))
22 format (3(ES14.7, 2x))
23 format (20(ES14.7, 2x))

    write(21, 12) ind_h 
    do i=1, ind_h
      write(21, 21) (hy_ide(i, j), j=1,13), (hy_mhd(i, j), j=1,nit_mhd)
      if(derivatives_out) then
        write(23,23) (hy_deriv(i,j), j=1,nderivatives)
      end if
    end do
      close(21)
      if(derivatives_out) then
        close(23)
      end if
  return
end subroutine print_hypersurface_mhd
! *****************************************************************************
subroutine save_hypersurface_mhd(ind_hys, ide_fo, mhd_fo, si, der_fo)
use viscous
  implicit none
  integer, intent(inout):: ind_hys 
  real(8), intent(in), dimension (1:13) :: ide_fo
  real(8), intent(in), dimension (1:4) ::  mhd_fo !only 3 are used when the charge density in the comoving frame is not printed
  real(8), intent(in), dimension (1:nderivatives) :: der_fo !derivatives for derivatives computation  
  real(8), intent(in) :: si
  
  real(8) temperature, dens, press, energy, entropy,vx,vy,vz
  real(8) time,glf,gcov3,zz
  integer errors
  integer i

  real(8) dutdt, dutdx, dutdy, dutdz, duxdt, duxdx, duxdy, duxdz, duydt, duydx, duydy, duydz, duzdt, duzdx, duzdy, duzdz, dtedt,&
      &dtedx, dtedy, dtedz
 
  errors=0
  ind_hys=ind_hys+1
  do i=1, nit_mhd
    hy_mhd(ind_hys, i)=mhd_fo(i)
  end do
  do i=1,9
    hy_ide(ind_hys, i)=ide_fo(i)
  end do
  do i=10,13
    hy_ide(ind_hys, i)=si*ide_fo(i)    
  end do
 
  if(derivatives_out) then
    dutdt=der_fo(dtt)
    dutdx=der_fo(dtx)
    dutdy=der_fo(dty)
    dutdz=der_fo(dtz)
    duxdt=der_fo(dxt)
    duxdx=der_fo(dxx)
    duxdy=der_fo(dxy)
    duxdz=der_fo(dxz)
    duydt=der_fo(dyt)
    duydx=der_fo(dyx)
    duydy=der_fo(dyy)
    duydz=der_fo(dyz)
    duzdt=der_fo(dzt)
    duzdx=der_fo(dzx)
    duzdy=der_fo(dzy)
    duzdz=der_fo(dzz)
    dtedt=der_fo(dtet)
    dtedx=der_fo(dtex)
    dtedy=der_fo(dtey)
    dtedz=der_fo(dtez)

    dens=ide_fo(5)
    press=ide_fo(9)
    
    temperature=freeze_value
    time=ide_fo(1)
    zz=ide_fo(4)
    vx=ide_fo(6)
    vy=ide_fo(7)
    vz=ide_fo(8)
    
    if(coordinates .eq. MINKOWSKI) then
      gcov3=1.
    else
      gcov3=time*time 
    end if
    
    glf=1./sqrt(1.-vx*vx-vy*vy-vz*vz*gcov3)
    do i=1,nderivatives
      hy_deriv(ind_hys,i)=der_fo(i)
    end do
  end if
  return
end subroutine save_hypersurface_mhd

! *****************************************************************************

subroutine set_fo_mhd(it, ix,iy,iz, iit, iix, iiy, iiz, fo_ide, fo_mhd, fo_der, dir_flag)
use common, only : x,y,z, g_cov, gp,  gm, &
  & krh, kvx, kvy, kvz, kpr, d_val &
  ,kbx,kby,kbz

implicit none
! passed variables
!  WARNING  remember to pass always v_t2 as current and v_backup as backup
!  when we are checking on directins different  than time 
!  v_backup and v_t2 are passed as the same  field (backup)
!  and the indexes change (for instance ix=1  iix=2)
!   real(8), intent(in), dimension (nx, ny, nz, nv) :: v_backup, v_t2
  integer, intent(in) :: it, ix, iy, iz, iit, iix, iiy, iiz
  
  real(8), intent(out), dimension (1:13) :: fo_ide
  real(8), intent(out), dimension (1:4) :: fo_mhd 
  real(8), intent(out), dimension (1:nderivatives) :: fo_der
  integer, intent(in) :: dir_flag

  integer :: idx !just a counter
  
  ! subroutine variables
  real(8) weight, betasq
  
  real(8) dx, dy, dz, dt
  real(8) ddx, ddy, ddz, ddt
 
  weight= (field(it,ix,iy,iz)-freeze_value)/(field(it,ix,iy,iz)-field(iit,iix,iiy,iiz))
  ! tau  or t 
  fo_ide(1)=time_arr(it)+weight*(time_arr(iit)-time_arr(it))
  ! x
  fo_ide(2)=x(ix)+weight*(x(iix)-x(ix))
  ! y
  fo_ide(3)=y(iy)+weight*(y(iiy)-y(iy))
  ! eta or z
  fo_ide(4)=z(iz)+weight*(z(iiz)-z(iz))


  if(dir_flag .eq. tau_dir) then
    !rho
    fo_ide(5)=v_backup(ix,iy,iz,krh)+weight*(primitives(iix,iiy,iiz,krh)-v_backup(ix,iy,iz,krh))
    !vx
    fo_ide(6)=v_backup(ix,iy,iz,kvx)+weight*(primitives(iix,iiy,iiz,kvx)-v_backup(ix,iy,iz,kvx))
    !vy
    fo_ide(7)=v_backup(ix,iy,iz,kvy)+weight*(primitives(iix,iiy,iiz,kvy)-v_backup(ix,iy,iz,kvy))
    !vz
    fo_ide(8)=v_backup(ix,iy,iz,kvz)+weight*(primitives(iix,iiy,iiz,kvz)-v_backup(ix,iy,iz,kvz))
    !prex
    fo_ide(9)=v_backup(ix,iy,iz,kpr)+weight*(primitives(iix,iiy,iiz,kpr)-v_backup(ix,iy,iz,kpr))
    
    ! components of the B field vector
    ! bx
    fo_mhd(1)=v_backup(ix,iy,iz,kbx)+weight*(primitives(iix,iiy,iiz,kbx)-v_backup(ix,iy,iz,kbx))    
    ! by
    fo_mhd(2)=v_backup(ix,iy,iz,kby)+weight*(primitives(iix,iiy,iiz,kby)-v_backup(ix,iy,iz,kby))
    ! bz
    fo_mhd(3)=v_backup(ix,iy,iz,kbz)+weight*(primitives(iix,iiy,iiz,kbz)-v_backup(ix,iy,iz,kbz))
    if(print_rho_comov) then
      ! rho in comoving frame
      fo_mhd(4)=v_backup(ix,iy,iz,krc)+weight*(primitives(iix,iiy,iiz,krc)-v_backup(ix,iy,iz,krc))
    end if  
  else
    !rho
    fo_ide(5)=v_backup(ix,iy,iz,krh)+weight*(v_backup(iix,iiy,iiz,krh)-v_backup(ix,iy,iz,krh))  
    !vx
    fo_ide(6)=v_backup(ix,iy,iz,kvx)+weight*(v_backup(iix,iiy,iiz,kvx)-v_backup(ix,iy,iz,kvx))
    !vy
    fo_ide(7)=v_backup(ix,iy,iz,kvy)+weight*(v_backup(iix,iiy,iiz,kvy)-v_backup(ix,iy,iz,kvy))
    !vz
    fo_ide(8)=v_backup(ix,iy,iz,kvz)+weight*(v_backup(iix,iiy,iiz,kvz)-v_backup(ix,iy,iz,kvz))
    !prex
    fo_ide(9)=v_backup(ix,iy,iz,kpr)+weight*(v_backup(iix,iiy,iiz,kpr)-v_backup(ix,iy,iz,kpr))

    ! components of the B field vector
    ! bx
    fo_mhd(1)=v_backup(ix,iy,iz,kbx)+weight*(v_backup(iix,iiy,iiz,kbx)-v_backup(ix,iy,iz,kbx))    
    ! by
    fo_mhd(2)=v_backup(ix,iy,iz,kby)+weight*(v_backup(iix,iiy,iiz,kby)-v_backup(ix,iy,iz,kby))
    ! bz
    fo_mhd(3)=v_backup(ix,iy,iz,kbz)+weight*(v_backup(iix,iiy,iiz,kbz)-v_backup(ix,iy,iz,kbz))
    if(print_rho_comov) then
      ! rho in comoving frame
      fo_mhd(4)=v_backup(ix,iy,iz,krc)+weight*(v_backup(iix,iiy,iiz,krc)-v_backup(ix,iy,iz,krc))
    end if  
  end if

  betasq=fo_ide(6)*fo_ide(6)*(g_cov(1))&
	+fo_ide(7)*fo_ide(7)*(g_cov(2))&
	+fo_ide(8)*fo_ide(8)*(g_cov(3))
    if (betasq>1.0) then 
      call renormalize(fo_ide(6:8), fo_ide(1:4))
      print *, "you have to fix this"
      call exit(1)
    endif

  
  dt=abs(time_arr(2)-time_arr(1))
  dx= d_val(1)
  dy= d_val(2)
  dz= d_val(3)

! dV perp t
   fo_ide(10)=real((iit-it),8)*dx*dy*dz*fo_ide(1)
! dV perp x
   fo_ide(11)=real((iix-ix),8)*dt*dy*dz*fo_ide(1)
! dV perp y
   fo_ide(12)=real((iiy-iy),8)*dx*dt*dz*fo_ide(1)
! dV perp z
   fo_ide(13)=real((iiz-iz),8)*dx*dy*dt*fo_ide(1)
   
 ! HERE STARTS THE SECTION ADDED FOR STORING DERIVATIVES NEEDED FOR COMPUTING TERMAL VORTICITY
  if(derivatives_out) then
    if(dir_flag .eq. tau_dir) then
      do idx=1,nderivatives
         fo_der(idx)=der_backup(ix,iy,iz,idx)+weight*(derivatives_all(iix,iiy,iiz,idx)-der_backup(ix,iy,iz,idx))
      end do
    else
      do idx=1,nderivatives
         fo_der(idx)=der_backup(ix,iy,iz,idx)+weight*(der_backup(iix,iiy,iiz,idx)-der_backup(ix,iy,iz,idx))
      end do
    end if
  end if

  return
end subroutine set_fo_mhd

! !******************************************************************


! *****************************************************************************
subroutine find_hypersurface_visco()
implicit none
!  subroutine variables
real(8) sign_
integer ix, iy, iz, iv
real(8), dimension(0:3):: xfo, vfo, dVfo

real(8), dimension (1:13) :: info_ide
real(8), dimension (1:7) ::  info_vis !!Here we allocate it differently
real(8), dimension (1:nderivatives) :: info_der !used for derivatives computation
  
integer i_hysu

  i_hysu=0
  !we check that only processor 0 is calling this subroutine
  call check_processor()
  do ix=1, nx
    do iy=1, ny
      do iz=1, nz
	! check on tau (t)---------------------------------------
	if ((field(2,ix,iy,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) > 0.0) then ! FO 
	  sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(2,ix,iy,iz)))
	  call set_fo_vis(1,ix,iy,iz, 2,ix,iy,iz, info_ide, info_vis, info_der, tau_dir)
	  call save_hypersurface_visco(i_hysu, info_ide, info_vis, sign_, info_der)
	endif
	! check on x -------------------------------------------
	if (ix<nx) then
	  if ((field(1,ix+1,iy,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix+1,iy,iz)))
	    call set_fo_vis(1,ix,iy,iz, 1,ix+1,iy,iz, info_ide, info_vis, info_der, x_dir)
	    call save_hypersurface_visco(i_hysu, info_ide, info_vis, sign_, info_der)
	  endif
	endif
	! check on y -------------------------------------------
	if (iy<ny) then
	  if ((field(1,ix,iy+1,iz)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix,iy+1,iz)))
	    call set_fo_vis(1,ix,iy,iz, 1,ix,iy+1,iz, info_ide, info_vis, info_der, y_dir)
	    call save_hypersurface_visco(i_hysu, info_ide, info_vis, sign_, info_der)
	  endif
	endif
	! check on eta (z) -------------------------------------
	if (iz<nz) then
	  if ((field(1,ix,iy,iz+1)-freeze_value)*(freeze_value-field(1,ix,iy,iz)) >0.0) then ! FO 
            sign_=1.0d0*sign(1.0d0,(field(1,ix,iy,iz)-field(1,ix,iy,iz+1)))
	    call set_fo_vis(1,ix,iy,iz, 1,ix,iy,iz+1, info_ide, info_vis, info_der, z_dir)
	    call save_hypersurface_visco(i_hysu, info_ide, info_vis, sign_, info_der)
	  endif
	endif
      end do
    end do
  end do
  
  call print_hypersurface_visco(i_hysu)
  return
end subroutine find_hypersurface_visco
! *****************************************************************************
subroutine print_hypersurface_visco(ind_h)
  use common, only: kpixx, kpiyy, kpizz, kpixy, kpixz, kpiyz, nv, kvx,kvy,kvz
  use eos
!   use viscous
  implicit none 

  integer, intent(in) :: ind_h

  integer fileErr, i, j, erralert
  real(8) temp, rho, eng, prex, entropy, fracden

  real(8) titt, titx, tity, titz
  real(8) tizz
  real(8) time_fo,gcv3

  real(8), dimension (nv) :: fake_v

  open (21,file=prefix_dir//'hypersurface.txt',form='formatted',status='old', access='append', iostat=fileErr)
    if (fileErr .ne. 0) then 
      write (*,*) " I cannot open the file ", prefix_dir//'hypersurface.txt' 
      write (*,*) " I have to quit "
      call exit(1)
    endif
  if(derivatives_out) then
   open (23,file=prefix_dir//'hypersurf_deriv.txt',form='formatted',status='old', access='append', iostat=fileErr)
    if (fileErr .ne. 0) then
      write (*,*) " I cannot open the file ", prefix_dir//'hypersurf_deriv.txt'
      write (*,*) " I have to quit "
      call exit(1)
    endif
  endif 

12 format (i6)
21 format (25(ES14.7, 2x))
22 format (3(ES14.7, 2x))
23 format (20(ES14.7, 2x))

    write(21, 12) ind_h 
    do i=1, ind_h
      rho=hy_ide(i,5)
      prex=hy_ide(i,9)
      call get_derived_data(rho, prex, eng, temp, entropy, erralert)

      fake_v=0.0
      time_fo=hy_ide(i, 1)
      gcv3=time_fo*time_fo
      fake_v(kvx)=hy_ide(i, 6)
      fake_v(kvy)=hy_ide(i, 7)
      fake_v(kvz)=hy_ide(i, 8)
      fake_v(kpixx)=hy_vis(i, 5)
      fake_v(kpiyy)=hy_vis(i, 6)
      fake_v(kpizz)=hy_vis(i, 7)
      fake_v(kpixy)=hy_vis(i, 2)
      fake_v(kpixz)=hy_vis(i, 3)
      fake_v(kpiyz)=hy_vis(i, 4)
      if(obtained .eq. 'no') then
        call get_derived_pi_gcv3(fake_v,titt,titx,tity,titz,gcv3)
      else
        call get_derived_pi_zz_gcv3(fake_v,titt,titx,tity,titz,tizz,gcv3)
        fake_v(kpizz)=tizz
      end if 

      fracden=1.0/(temp*temp*(eng+prex))
      write(21, 21) (hy_ide(i, j), j=1,13), (hy_vis(i, j), j=1,7), titt, titx, tity, titz, fracden
      if(derivatives_out) then
        write(23,23) (hy_deriv(i,j), j=1,nderivatives)
      end if
    end do
      close(21)
      if(derivatives_out) then
        close(23)
      end if
  return
end subroutine print_hypersurface_visco
! *****************************************************************************
subroutine save_hypersurface_visco(ind_hys, ide_fo, vis_fo, si, der_fo)
use viscous
  implicit none
  integer, intent(inout):: ind_hys 
  real(8), intent(in), dimension (1:13) :: ide_fo
  real(8), intent(in), dimension (1:7) ::  vis_fo
  real(8), intent(in), dimension (1:nderivatives) :: der_fo !derivatives for derivatives computation  
  real(8), intent(in) :: si
  
  real(8) temperature, dens, press, energy, entropy, vx, vy, vz 
  real(8) time, zz, glf,gcov3       
  integer errors
  integer i

  real(8) dutdt, dutdx, dutdy, dutdz, duxdt, duxdx, duxdy, duxdz, duydt, duydx, duydy, duydz, duzdt, duzdx, duzdy, duzdz, dtedt,&
      &dtedx, dtedy, dtedz
  
  errors=0
  ind_hys=ind_hys+1
  do i=1, 7
    hy_ide(ind_hys, i)=ide_fo(i)
    hy_vis(ind_hys, i)=vis_fo(i)
  end do
  hy_ide(ind_hys, 8)=ide_fo(8)
  hy_ide(ind_hys, 9)=ide_fo(9)
  do i=10,13
    hy_ide(ind_hys, i)=si*ide_fo(i)    
  end do
 
  if(derivatives_out) then
    dutdt=der_fo(dtt)
    dutdx=der_fo(dtx)
    dutdy=der_fo(dty)
    dutdz=der_fo(dtz)
    duxdt=der_fo(dxt)
    duxdx=der_fo(dxx)
    duxdy=der_fo(dxy)
    duxdz=der_fo(dxz)
    duydt=der_fo(dyt)
    duydx=der_fo(dyx)
    duydy=der_fo(dyy)
    duydz=der_fo(dyz)
    duzdt=der_fo(dzt)
    duzdx=der_fo(dzx)
    duzdy=der_fo(dzy)
    duzdz=der_fo(dzz)
    dtedt=der_fo(dtet)
    dtedx=der_fo(dtex)
    dtedy=der_fo(dtey)
    dtedz=der_fo(dtez)

    dens=ide_fo(5)
    press=ide_fo(9)
    
    temperature=freeze_value
    time=ide_fo(1)
    zz=ide_fo(4)
    vx=ide_fo(6)
    vy=ide_fo(7)
    vz=ide_fo(8)
    
    if(coordinates .eq. MINKOWSKI) then
      gcov3=1.
    else
      gcov3=time*time 
    end if
    
    glf=1./sqrt(1.-vx*vx-vy*vy-vz*vz*gcov3)
    do i=1,nderivatives
      hy_deriv(ind_hys,i)=der_fo(i)
    end do
  end if
  return
end subroutine save_hypersurface_visco

! *****************************************************************************

subroutine set_fo_vis(it, ix,iy,iz, iit, iix, iiy, iiz, fo_ide, fo_vis, fo_der, dir_flag)
use common, only : x,y,z, g_cov, gp,  gm, &
  & krh, kvx, kvy, kvz, kpr, d_val &
! we need other terms
  ,kpibu, kpixy, kpixz,kpiyz,kpixx,kpiyy,kpizz

implicit none
! passed variables
!  WARNING  remember to pass always v_t2 as current and v_backup as backup
!  when we are checking on directins different  than time 
!  v_backup and v_t2 are passed as the same  field (backup)
!  and the indexes change (for instance ix=1  iix=2)
!   real(8), intent(in), dimension (nx, ny, nz, nv) :: v_backup, v_t2
  integer, intent(in) :: it, ix, iy, iz, iit, iix, iiy, iiz
  
  real(8), intent(out), dimension (1:13) :: fo_ide
  real(8), intent(out), dimension (1:7) :: fo_vis 
  real(8), intent(out), dimension (1:nderivatives) :: fo_der
  integer, intent(in) :: dir_flag

  integer :: idx !just a counter
  
  ! subroutine variables
  real(8) weight, betasq
  
  real(8) dx, dy, dz, dt
  real(8) ddx, ddy, ddz, ddt
 
  weight= (field(it,ix,iy,iz)-freeze_value)/(field(it,ix,iy,iz)-field(iit,iix,iiy,iiz))
  ! tau  or t 
  fo_ide(1)=time_arr(it)+weight*(time_arr(iit)-time_arr(it))
  ! x
  fo_ide(2)=x(ix)+weight*(x(iix)-x(ix))
  ! y
  fo_ide(3)=y(iy)+weight*(y(iiy)-y(iy))
  ! eta or z
  fo_ide(4)=z(iz)+weight*(z(iiz)-z(iz))


  if(dir_flag .eq. tau_dir) then
    !rho
    fo_ide(5)=v_backup(ix,iy,iz,krh)+weight*(primitives(iix,iiy,iiz,krh)-v_backup(ix,iy,iz,krh))
    !vx
    fo_ide(6)=v_backup(ix,iy,iz,kvx)+weight*(primitives(iix,iiy,iiz,kvx)-v_backup(ix,iy,iz,kvx))
    !vy
    fo_ide(7)=v_backup(ix,iy,iz,kvy)+weight*(primitives(iix,iiy,iiz,kvy)-v_backup(ix,iy,iz,kvy))
    !vz
    fo_ide(8)=v_backup(ix,iy,iz,kvz)+weight*(primitives(iix,iiy,iiz,kvz)-v_backup(ix,iy,iz,kvz))
    !prex
    fo_ide(9)=v_backup(ix,iy,iz,kpr)+weight*(primitives(iix,iiy,iiz,kpr)-v_backup(ix,iy,iz,kpr))
    
    ! bulk viscosity
    fo_vis(1)=v_backup(ix,iy,iz,kpibu)+weight*(primitives(iix,iiy,iiz,kpibu)-v_backup(ix,iy,iz,kpibu))    
    ! components of the shear stress tensor
    ! xy
    fo_vis(2)=v_backup(ix,iy,iz,kpixy)+weight*(primitives(iix,iiy,iiz,kpixy)-v_backup(ix,iy,iz,kpixy))
    ! xz
    fo_vis(3)=v_backup(ix,iy,iz,kpixz)+weight*(primitives(iix,iiy,iiz,kpixz)-v_backup(ix,iy,iz,kpixz))
    ! yz
    fo_vis(4)=v_backup(ix,iy,iz,kpiyz)+weight*(primitives(iix,iiy,iiz,kpiyz)-v_backup(ix,iy,iz,kpiyz))
    ! xx
    fo_vis(5)=v_backup(ix,iy,iz,kpixx)+weight*(primitives(iix,iiy,iiz,kpixx)-v_backup(ix,iy,iz,kpixx))
    ! yy
    fo_vis(6)=v_backup(ix,iy,iz,kpiyy)+weight*(primitives(iix,iiy,iiz,kpiyy)-v_backup(ix,iy,iz,kpiyy))
    ! zz 
    fo_vis(7)=v_backup(ix,iy,iz,kpizz)+weight*(primitives(iix,iiy,iiz,kpizz)-v_backup(ix,iy,iz,kpizz))    
  else
    !rho
    fo_ide(5)=v_backup(ix,iy,iz,krh)+weight*(v_backup(iix,iiy,iiz,krh)-v_backup(ix,iy,iz,krh))  
    !vx
    fo_ide(6)=v_backup(ix,iy,iz,kvx)+weight*(v_backup(iix,iiy,iiz,kvx)-v_backup(ix,iy,iz,kvx))
    !vy
    fo_ide(7)=v_backup(ix,iy,iz,kvy)+weight*(v_backup(iix,iiy,iiz,kvy)-v_backup(ix,iy,iz,kvy))
    !vz
    fo_ide(8)=v_backup(ix,iy,iz,kvz)+weight*(v_backup(iix,iiy,iiz,kvz)-v_backup(ix,iy,iz,kvz))
    !prex
    fo_ide(9)=v_backup(ix,iy,iz,kpr)+weight*(v_backup(iix,iiy,iiz,kpr)-v_backup(ix,iy,iz,kpr))
    
    ! bulk viscosity
    fo_vis(1)=v_backup(ix,iy,iz,kpibu)+weight*(v_backup(iix,iiy,iiz,kpibu)-v_backup(ix,iy,iz,kpibu)) 
    ! components of the shear stress tensor
    ! xy
    fo_vis(2)=v_backup(ix,iy,iz,kpixy)+weight*(v_backup(iix,iiy,iiz,kpixy)-v_backup(ix,iy,iz,kpixy))
    ! xz
    fo_vis(3)=v_backup(ix,iy,iz,kpixz)+weight*(v_backup(iix,iiy,iiz,kpixz)-v_backup(ix,iy,iz,kpixz))
    ! yz
    fo_vis(4)=v_backup(ix,iy,iz,kpiyz)+weight*(v_backup(iix,iiy,iiz,kpiyz)-v_backup(ix,iy,iz,kpiyz))
    ! xx
    fo_vis(5)=v_backup(ix,iy,iz,kpixx)+weight*(v_backup(iix,iiy,iiz,kpixx)-v_backup(ix,iy,iz,kpixx))
    ! yy
    fo_vis(6)=v_backup(ix,iy,iz,kpiyy)+weight*(v_backup(iix,iiy,iiz,kpiyy)-v_backup(ix,iy,iz,kpiyy))
    ! zz 
    fo_vis(7)=v_backup(ix,iy,iz,kpizz)+weight*(v_backup(iix,iiy,iiz,kpizz)-v_backup(ix,iy,iz,kpizz))      
  end if

  betasq=fo_ide(6)*fo_ide(6)*(g_cov(1))&
	+fo_ide(7)*fo_ide(7)*(g_cov(2))&
	+fo_ide(8)*fo_ide(8)*(g_cov(3))
    if (betasq>1.0) then 
      call renormalize(fo_ide(6:8), fo_ide(1:4))
      print *, "you have to fix this"
      call exit(1)
    endif

  
  dt=abs(time_arr(2)-time_arr(1))
  dx= d_val(1)
  dy= d_val(2)
  dz= d_val(3)

! dV perp t
   fo_ide(10)=real((iit-it),8)*dx*dy*dz*fo_ide(1)
! dV perp x
   fo_ide(11)=real((iix-ix),8)*dt*dy*dz*fo_ide(1)
! dV perp y
   fo_ide(12)=real((iiy-iy),8)*dx*dt*dz*fo_ide(1)
! dV perp z
   fo_ide(13)=real((iiz-iz),8)*dx*dy*dt*fo_ide(1)
   
 ! HERE STARTS THE SECTION ADDED FOR STORING DERIVATIVES NEEDED FOR COMPUTING TERMAL VORTICITY
  if(derivatives_out) then
    if(dir_flag .eq. tau_dir) then
      do idx=1,nderivatives
         fo_der(idx)=der_backup(ix,iy,iz,idx)+weight*(derivatives_all(iix,iiy,iiz,idx)-der_backup(ix,iy,iz,idx))
      end do
    else
      do idx=1,nderivatives
         fo_der(idx)=der_backup(ix,iy,iz,idx)+weight*(der_backup(iix,iiy,iiz,idx)-der_backup(ix,iy,iz,idx))
      end do
    end if
  end if

  return
end subroutine set_fo_vis

! !******************************************************************
! !******************************************************************
! !******************************************************************
      subroutine renormalize(vel, coord)
      use common, only :g_cov
      implicit none
	real(8), dimension(1:3), intent(inout) :: vel
	real(8), dimension(0:3), intent(in) :: coord
	real(8) modulus
	integer i 
	real(8), parameter :: epsil=1e-16
	modulus=vel(1)*vel(1)*(g_cov(1)) + vel(2)*vel(2)*(g_cov(2)) &
		+ vel(3)*vel(3)*(g_cov(3)) - epsil
	do i=1,3
	  vel(i)=vel(i)/modulus
	end do
      end subroutine renormalize
! !*****************************************************************
subroutine get_hysufile_ready()
  use eos
  use viscous
  use common, only:restarted_run
  implicit none 
  integer fileErr
  integer erralert
  
  real(8) rho_f, prex_f, eng_f, temp_f, entropy_f,dprdrh_f,dprden_f
 
  if(restarted_run) then
    open (21,file=prefix_dir//'hypersurface.txt',form='formatted',status="old", position="append", action="write",iostat=fileErr)
  else
    open (21,file=prefix_dir//'hypersurface.txt',form='formatted',status='replace', iostat=fileErr)
  end if
    if (fileErr .ne. 0) then 
       write (*,*) " I cannot open the file ", prefix_dir//'hypersurface.txt' 
       write (*,*) " I have to quit "
       call exit(1)
    endif
      
! !     eos_energy		in :: pr,rh // out :: en
! !     eos_pressure		in :: rh,en // out :: pr,dprdrh,dprden
! !     eos_pressure0		in :: pr_input // out :: en_output
! !     get_derived_data	in :: rh, pressure //  out :: energy_density, temperature, entropy_density 
! !     find_pressure_from_temperature 		in :: temp_input // out :: press_output
       
   select case (freeze_type) 
    case(0)
!     Have temperature, want energy density
      temp_f=freeze_value
      call find_pressure_from_temperature(temp_f, prex_f)
      call eos_energy(rho_f,eng_f, prex_f,erralert)
    case(1)
!     Have  energy density, want temperature
      eng_f=freeze_value
      call eos_pressure(rho_f,eng_f,prex_f,dprdrh_f,dprden_f,erralert)
      call get_derived_data(rho_f, prex_f, eng_f, temp_f, entropy_f,erralert )
    case default
      print *, "The switch value for the freeze-out hypersurface can only be"
      print *, "0 - constant temperature, 1- constant energy "
      print *, "Your choice was", freeze_type
      print *, "I am going to quit"
      call exit 
   end select
   
! save in the hypersurface file the criterion, the energy and the temperature  
   if(.not. restarted_run) write (21,*) freeze_type, temp_f, eng_f
   close(21)

   if(derivatives_out) then
    if(restarted_run) then
     open(23,file=prefix_dir//'hypersurf_deriv.txt',form='formatted',status="old", position="append", action="write",iostat=fileErr)
    else
     open(23,file=prefix_dir//'hypersurf_deriv.txt',form='formatted',status='replace', iostat=fileErr)
    end if
     if (fileErr .ne. 0) then 
	write (*,*) " I cannot open the file ", prefix_dir//'hypersurf_deriv.txt' 
	write (*,*) " I have to quit "
	call exit(1)
     endif
    endif

end subroutine get_hysufile_ready

! !*****************************************************************

end module hypersurface
