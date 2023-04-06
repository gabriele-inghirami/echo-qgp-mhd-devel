! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  File: work.f90                                                           *
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

module work
!-- Riemann solver and interface for high-order procedures

  use common, only: derivatives_out,aflux

  implicit none

  real(8),allocatable,dimension(:,:) :: w0,w1,wl,wr,wd,wdr,wdl

contains

! *****************************************************************************

subroutine work_alloc(i1,i2,nv)
!-- Allocate main 1-D arrays

  integer,intent(in) :: i1,i2,nv

  allocate(w0(i1:i2,nv),wl(i1:i2,nv),wr(i1:i2,nv),w1(i1:i2,nv+2))
  if(derivatives_out) then 
    allocate(wd(i1:i2,5),wdl(i1:i2,5),wdr(i1:i2,5))
  else 
    allocate(wd(i1:i2,4),wdl(i1:i2,4),wdr(i1:i2,4))
  end if
  w0=0.;wl=0.;wr=0.;w1=0.;wd=0.;wdl=0.;wdr=0.;
  

end subroutine work_alloc

! *****************************************************************************

subroutine work_glmx(n)
!-- Calculate fluid upwind fluxes

  use common,only: solver
  use system_eqgp

  implicit none

  integer               ,intent(in) :: n
  integer :: i


    do i=0,n
   
      wl(i,kglm)=0.5*(wl(i,kglm)+wr(i,kglm))-(wr(i,kbx)-wl(i,kbx))*glm_ch/2
      wr(i,kglm)=wl(i,kglm)
      wl(i,kbx)=0.5*(wl(i,kbx)+wr(i,kbx))-(wr(i,kglm)-wl(i,kglm))/(2*glm_ch)
      wr(i,kbx)=wl(i,kbx)
    end do


end subroutine work_glmx

! *****************************************************************************

subroutine work_glmy(n)
!-- Calculate fluid upwind fluxes

  use common,only: solver
  use system_eqgp

  implicit none

  integer               ,intent(in) :: n
  integer :: i

    do i=0,n
   
      wl(i,kglm)=0.5*(wl(i,kglm)+wr(i,kglm))-(wr(i,kby)-wl(i,kby))*glm_ch/2
      wr(i,kglm)=wl(i,kglm)
      wl(i,kby)=0.5*(wl(i,kby)+wr(i,kby))-(wr(i,kglm)-wl(i,kglm))/(2*glm_ch)
      wr(i,kby)=wl(i,kby)

    end do


end subroutine work_glmy

! *****************************************************************************

subroutine work_glmz(n)
!-- Calculate fluid upwind fluxes

  use common,only: solver
  use system_eqgp

  implicit none

  integer               ,intent(in) :: n
  integer :: i

    do i=0,n
   
      wl(i,kglm)=0.5*(wl(i,kglm)+wr(i,kglm))-(wr(i,kbz)-wl(i,kbz))*glm_ch/2
      wr(i,kglm)=wl(i,kglm)
      wl(i,kbz)=0.5*(wl(i,kbz)+wr(i,kbz))-(wr(i,kglm)-wl(i,kglm))/(2*glm_ch)
      wr(i,kbz)=wl(i,kbz)

    end do



end subroutine work_glmz

! *****************************************************************************

subroutine work_fupx(n,x1,x2,x3,errcode)
!-- Calculate fluid upwind fluxes

  use common,only: solver
  use system_eqgp

  integer               ,intent(in) :: n
  real(8)   ,dimension(0:n),intent(in) :: x1
  real(8)                  ,intent(in) :: x2,x3

  real(8),dimension( 3) :: x
  real(8),dimension(nv) :: ul,ur,fl,fr
  real(8),dimension( 2) :: vfl,vfr

  integer :: i
  real(8)    :: ap,am,a1,a

  real(8),dimension(nv) :: u_hll,f_hll,us
  real(8) :: vzl,vzr,sr,sl,aa,bb,cc,ss,ps
  integer :: errcode

  errcode=0

  select case(solver)

!-- Local Lax Friedrichs (Rusanov) symmetric two-wave solver
  case('LLF')

    do i=0,n
      call system_flux(wl(i,1:nv),ul,fl,vfl,1,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupx subroutine when computing left fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1(i),"x2=",x2,"x3=",x3
        return
      end if
      call system_flux(wr(i,1:nv),ur,fr,vfr,1,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupx subroutine when computing right fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1(i),"x2=",x2,"x3=",x3
        return
      end if
      ap=max(0., vfl(1), vfr(1))
      am=max(0.,-vfl(2),-vfr(2))
      a =max(ap,am)
      w1(i,1:nv)=0.5*(fl+fr-a*(ur-ul))
      w1(i,nv+1)=a
      w1(i,nv+2)=a
    end do

!-- Harten - Lax - van Leer upwind two-wave solver
  case('HLL')
  
    do i=0,n
   
      call system_flux(wl(i,1:nv),ul,fl,vfl,1,errcode)
      !DDD write(*,*) "Call left", fl(9)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupx subroutine when computing left fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1(i),"x2=",x2,"x3=",x3
        return
      end if
      call system_flux(wr(i,1:nv),ur,fr,vfr,1,errcode)
      !DDD write(*,*) "Call right", fr(9)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupx subroutine when computing right fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1(i),"x2=",x2,"x3=",x3
        return
      end if
      ap=max(0., vfl(1), vfr(1))
      am=max(0.,-vfl(2),-vfr(2))
      aflux(1)=max(ap,am)

      if ((ap .eq. 0.0) .and. (am .eq. 0.0)) then
      a1=1.
      else
      a1=1./(ap+am)
      end if
!        w1(i,kbx)=0.5*(ul(kbx)+ur(kbx))-(ur(kglm)+ul(kglm))/(2*glm_ch)
!        w1(i,kglm)=0.5*(ul(kglm)+ur(kglm))-(ur(kbx)+ul(kbx))*glm_ch/2
      w1(i,1:nv)=a1*(ap*fl+am*fr-ap*am*(ur-ul))
      w1(i,nv+1)=ap
      w1(i,nv+2)=am
    end do

  end select

end subroutine work_fupx

! *****************************************************************************

subroutine work_fupy(n,x1,x2,x3,errcode)
!-- Calculate fluid upwind fluxes

  use common,only: solver
  use system_eqgp

  integer               ,intent(in) :: n
  real(8)   ,dimension(0:n),intent(in) :: x2
  real(8)                  ,intent(in) :: x1,x3

  real(8),dimension( 3) :: x
  real(8),dimension(nv) :: ul,ur,fl,fr
  real(8),dimension( 2) :: vfl,vfr

  integer :: i
  real(8)    :: ap,am,a1,a

  real(8),dimension(nv) :: u_hll,f_hll,us
  real(8) :: vzl,vzr,sr,sl,aa,bb,cc,ss,ps
  integer :: errcode

  errcode=0

  select case(solver)

!-- Local Lax Friedrichs (Rusanov) symmetric two-wave solver
  case('LLF')

    do i=0,n
      call system_flux(wl(i,1:nv),ul,fl,vfl,2,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupy subroutine when computing left fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2(i),"x3=",x3
        return
      end if
      
      call system_flux(wr(i,1:nv),ur,fr,vfr,2,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupy subroutine when computing right fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2(i),"x3=",x3
        return
      end if
      ap=max(0., vfl(1), vfr(1))
      am=max(0.,-vfl(2),-vfr(2))
      a =max(ap,am)
      w1(i,1:nv)=0.5*(fl+fr-a*(ur-ul))
      w1(i,nv+1)=a
      w1(i,nv+2)=a
    end do

!-- Harten - Lax - van Leer upwind two-wave solver
  case('HLL')

    do i=0,n
      call system_flux(wl(i,1:nv),ul,fl,vfl,2,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupy subroutine when computing left fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2(i),"x3=",x3
        return
      end if
      call system_flux(wr(i,1:nv),ur,fr,vfr,2,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupy subroutine when computing right fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2(i),"x3=",x3
        return
      end if
      ap=max(0., vfl(1), vfr(1))
      am=max(0.,-vfl(2),-vfr(2))

      aflux(2)=max(ap,am)

      if ((ap .eq. 0.0) .and. (am .eq. 0.0)) then
      a1=1.
      else
      a1=1./(ap+am)
      end if

!        w1(i,kby)=0.5*(ul(kby)+ur(kby))-(ur(kglm)+ul(kglm))/(2*glm_ch)
!        w1(i,kglm)=0.5*(ul(kglm)+ur(kglm))-(ur(kby)+ul(kby))*glm_ch/2
      w1(i,1:nv)=a1*(ap*fl+am*fr-ap*am*(ur-ul))
      w1(i,nv+1)=ap
      w1(i,nv+2)=am
    end do

  end select

end subroutine work_fupy

! *****************************************************************************

subroutine work_fupz(n,x1,x2,x3,errcode)
!-- Calculate fluid upwind fluxes

  use common,only: solver
  use system_eqgp

  integer               ,intent(in) :: n
  real(8)   ,dimension(0:n),intent(in) :: x3
  real(8)                  ,intent(in) :: x1,x2

  real(8),dimension( 3) :: x
  real(8),dimension(nv) :: ul,ur,fl,fr
  real(8),dimension( 2) :: vfl,vfr

  integer :: i
  real(8)    :: ap,am,a1,a

  real(8),dimension(nv) :: u_hll,f_hll,us
  real(8) :: vzl,vzr,sr,sl,aa,bb,cc,ss,ps
  integer :: errcode

  errcode=0
  select case(solver)

!-- Local Lax Friedrichs (Rusanov) symmetric two-wave solver
  case('LLF')

    do i=0,n
      call system_flux(wl(i,1:nv),ul,fl,vfl,3,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupz subroutine when computing left fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2,"x3=",x3(i)
        return
      end if
      call system_flux(wr(i,1:nv),ur,fr,vfr,3,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupz subroutine when computing right fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2,"x3=",x3(i)
        return
      end if
      ap=max(0., vfl(1), vfr(1))
      am=max(0.,-vfl(2),-vfr(2))
      a =max(ap,am)
      w1(i,1:nv)=0.5*(fl+fr-a*(ur-ul))
      w1(i,nv+1)=a
      w1(i,nv+2)=a
    end do

!-- Harten - Lax - van Leer upwind two-wave solver
  case('HLL')

    do i=0,n
      call system_flux(wl(i,1:nv),ul,fl,vfl,3,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupz subroutine when computing left fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2,"x3=",x3(i)
        return
      end if
      call system_flux(wr(i,1:nv),ur,fr,vfr,3,errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into work_fupz subroutine when computing right fluxes."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: x1=",x1,"x2=",x2,"x3=",x3(i)
        return
      end if
      ap=max(0., vfl(1), vfr(1))
      am=max(0.,-vfl(2),-vfr(2))
      aflux(3)=max(ap,am)
      if ((ap .eq. 0.0) .and. (am .eq. 0.0)) then
      a1=1.
      else
      a1=1./(ap+am)
      end if
!        w1(i,kbz)=0.5*(ul(kbz)+ur(kbz))-(ur(kglm)+ul(kglm))/(2*glm_ch)
!        w1(i,kglm)=0.5*(ul(kglm)+ur(kglm))-(ur(kbz)+ul(kbz))*glm_ch/2
      w1(i,1:nv)=a1*(ap*fl+am*fr-ap*am*(ur-ul))
      w1(i,nv+1)=ap
      w1(i,nv+2)=am
    end do

  end select

end subroutine work_fupz

! *****************************************************************************

subroutine work_rec(n,k1,k2,flag)
!-- Interface for REC procedures

  use common,only: recal
  use holib

  integer,intent(in) :: n,k1,k2
  integer            :: ngc,k, flag

  select case(recal)
  case ('TVD2 ')
    ngc=2
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_rectvd2(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_rectvd2(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
   end if  
  case ('CENO3')
    ngc=3
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recceno3(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recceno3(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
   end if 
  case ('WENO3')
    ngc=2
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recweno3(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recweno3(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
    end if
  case ('WENO5')
    ngc=3
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recweno5(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recweno5(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
    end if 
  case ('PPM4 ')
    ngc=2
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recppm4(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recppm4(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
    end if
  case ('MPE3 ')
    ngc=3
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recmpe3(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recmpe3(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
    end if
 
  case ('MPE5 ')
    ngc=3
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recmpe5(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recmpe5(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
    end if
  case ('MPE7 ')
    ngc=4
    if (flag .eq. 0) then
      do k=k1,k2
        call holib_recmpe7(w0(1-ngc:n+ngc,k),wl(0:n,k),wr(0:n,k),n)
      end do
    else
      do k=k1,k2
        call holib_recmpe7(wd(1-ngc:n+ngc,k),wdl(0:n,k),wdr(0:n,k),n)
      end do
    end if
 
  case default
    write(*,*) 'REC undefined'
    call exit(1)
  end select

end subroutine work_rec

! *****************************************************************************

subroutine work_der(n,k1,k2,flag)
!-- Interface for DER procedures

  use common,only: der
  use holib

  integer,intent(in) :: n,k1,k2
  integer            :: ngc,k, flag

  select case(der)
  case ('    DER-E2')
    do k=k1,k2
      call holib_dere2(w0(0:n,k),w1(1:n,k),n)
    end do
  case ('    DER-E4')
    ngc=1
    do k=k1,k2
      call holib_dere4(w0(-ngc:n+ngc,k),w1(1:n,k),n)
    end do
  case ('    DER-E6')
    ngc=2
    do k=k1,k2
      call holib_dere6(w0(-ngc:n+ngc,k),w1(1:n,k),n)
    end do
  case ('    DER-E8')
    ngc=3
    do k=k1,k2
      call holib_dere8(w0(-ngc:n+ngc,k),w1(1:n,k),n)
    end do
  case default
    stop 'DER undefined'
  end select

end subroutine work_der

! ****************************************************************************

end module work

! *****************************************************************************
