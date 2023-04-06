! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2016,2018,2019 The ECHO-QGP team                      * 
! *                                                                           *
! *  File: evolve.f90                                                         *
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
! *           Gabriele Inghirami (gabriele.g.inghirami@jyu.fi)                *
! *                                                                           *
! *  Contributors:                                                            *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************

module evolve
!-- Main file
  use parallel
  use common, only: w,u0,du,ax,ay,az,wsend,wrecv,errcode,run_crashed,du_stiff

  implicit none


contains

! *****************************************************************************

subroutine evolve_alloc
!-- Allocate arrays

  use common
  use work,only: work_alloc
  
  call work_alloc(1-ngc,max(mx,my,mz)+ngc,nv)
    

  allocate( w(ix1-ngc:ix2+ngc,iy1-ngc:iy2+ngc,iz1-ngc:iz2+ngc,nv), &
           u0(ix1:ix2,iy1:iy2,iz1:iz2,nu),    &
           ax(ix0:ix2,iy0:iy2+1,iz0:iz2+1,2), &
           ay(ix0:ix2+1,iy0:iy2,iz0:iz2+1,2), &
           az(ix0:ix2+1,iy0:iy2+1,iz0:iz2,2))

  if(nrk .eq. SSP) then !if we are using imex schems we need more variables
    allocate(du(ix1:ix2,iy1:iy2,iz1:iz2,nu,4),    &
    &      du_stiff(ix1:ix2,iy1:iy2,iz1:iz2,nu,4))
  else
    allocate(du(ix1:ix2,iy1:iy2,iz1:iz2,nu,1),    &
    &      du_stiff(ix1:ix2,iy1:iy2,iz1:iz2,nu,1))
  end if
  u0=0
  w=0.
  du=0. 
  du_stiff=0. 
  ax=0.
  ay=0.
  az=0.

end subroutine evolve_alloc

! *****************************************************************************

subroutine evolve_glm_source(dt_glm)
use common
use system_eqgp

real(8) :: dt_glm, expfact
integer qq,hh

expfact=exp(-glm_ch**2*glm_alpha*dt_glm)

u(ix1:ix2,iy1:iy2,iz1:iz2,kglm)=u(ix1:ix2,iy1:iy2,iz1:iz2,kglm)*expfact

! DDD
!write(*,*)
!write(*,*) "******************************************"
!write(*,*) "glm_ch, glm_ch^2,glm_alpha,dt_glm,expfact"
!write(*,*) glm_ch,glm_ch**2,glm_alpha,dt_glm,expfact

!write(*,*)

!do hh=1,81
!   do qq=1,81
!      write(*,*) hh, qq, u(hh,qq,1,kglm), du(hh,qq,1,kglm)
!   end do
!end do
! DDD


end subroutine evolve_glm_source

! *****************************************************************************

subroutine evolve_main(irk,dt)
!-- Time evolution (RK sub-step)


  use common
  use system_eqgp
  integer :: num_out

  real(8),dimension(2,3,3) ::  &
    crk=reshape((/ 0.,1.,   0.,   0.,   0.,   0.,   &
                   0.,1.,1./2.,1./2.,   0.,   0.,   &
                   0.,1.,3./4.,1./4.,1./3.,2./3. /),shape=(/2,3,3/))

  real(8), parameter :: alpha_imex=0.24169426078821d0
  real(8), parameter :: beta_imex=0.06042356519705d0
  real(8), parameter :: eta_imex=0.12915286960590d0
  real(8), parameter :: delta_imex=0.5d0-beta_imex-eta_imex-alpha_imex

  integer,intent( in) :: irk
  real(8)   ,intent(out) :: dt
  
  real(8) :: c1,c2

  integer :: ix,iy,iz,ik
  
  real(8) :: dt_local

  integer :: errmpi
  logical :: not_crash_flag

  if(prl) check_crash=.false.

!DDD  write(*,*) "IRK is ", irk
!DDD  write(*,*) "u(kpr) at the beginning of the step is", u(2:3,2:3,:,kpr)

  if(nrk .ne. SSP) then !explicit RUNGE-KUTTA integration

    c1=crk(1,irk,nrk)
    c2=crk(2,irk,nrk)

!-- Calculate divergence of fluid fluxes and source terms ( -> du )

    du=0.
 
    glm_ch=glm_ch_old
    call evolve_dfx(1)     !-- dU=dU+dFx/dx
    call evolve_dfy(1)     !-- dU=dU+dFy/dy
    call evolve_dfz(1)     !-- dU=dU+dFz/dz
    call evolve_src(1)
    glm_ch_old=maxval(aflux)

    if(prl) then
       not_crash_flag=(.not. run_crashed)
       call MPI_Allreduce(not_crash_flag,check_crash,1,MPI_Logical,MPI_LAND,MPI_COMM_WORLD,errmpi)
       call MPI_Barrier(icomm,ierr)
       if(.not. check_crash) then
          run_crashed=.true.
          return
       end if
    else !serial run case
       if(run_crashed) return
    end if


!-- Calculate timestep ( -> dt )
    if (irk==1) then

       call evolve_dt(dt_local)
       dt_local=min(dt_local, maxdt)
       dt_local=min(dt_local, tout-t)
       if (abs(t+dt-tout) .lt. 1.e-9) then
         dt_local=tout-t
       else
         if((t+2.*dt_local - tout) .gt. 1.e-9) dt_local=(tout-t)/2.
       end if

       if(prl) then
         call MPI_Barrier(icomm,ierr)
         call MPI_Allreduce(dt_local, dt, 1, mpi_realtype, MPI_MIN, icomm,ierr)
       else
         dt=dt_local
       end if
       if(divclean) call  evolve_glm_source(dt/2.)
       u0=u
    end if  
    u=c1*u0+c2*(u-dt*du(:,:,:,:,1))
!-- Calculate primitive variables ( -> v )
    if(divclean .and. (irk .eq. nrk)) call  evolve_glm_source(dt/2.)
!DDD    write(*,*) "u(kpr) at the end of the step is", u(2:3,2:3,:,kpr)
    call evolve_u2v
  else
    du(:,:,:,:,irk)=0.
    du_stiff(:,:,:,:,irk)=0
    imex_alpha_coeff=alpha_imex

    if (irk==1) then

       call evolve_dt(dt_local)
       dt_local=min(dt_local, maxdt)
       dt_local=min(dt_local, tout-t)
       if (abs(t+dt-tout) .lt. 1.e-9) then
         dt_local=tout-t
       else
         if((t+2.*dt_local - tout) .gt. 1.e-9) dt_local=(tout-t)/2.
       end if

       if(prl) then
         call MPI_Barrier(icomm,ierr)
         call MPI_Allreduce(dt_local, dt, 1, mpi_realtype, MPI_MIN, icomm,ierr)
       else
         dt=dt_local
       end if
       u0=u 
       dt_int=dt
       glm_ch=glm_ch_old
       call evolve_dfx(irk)     !-- dU=dU+dFx/dx
       call evolve_dfy(irk)     !-- dU=dU+dFy/dy
       call evolve_dfz(irk)     !-- dU=dU+dFz/dz
       call evolve_src(irk)
       glm_ch_old=maxval(aflux)
       if(divclean) call  evolve_glm_source(dt/2.)
    else if(irk==2) then
       u=u0+dt*alpha_imex*du_stiff(:,:,:,:,1)
       call evolve_u2v
       glm_ch=glm_ch_old
       call evolve_dfx(irk)     !-- dU=dU+dFx/dx
       call evolve_dfy(irk)     !-- dU=dU+dFy/dy
       call evolve_dfz(irk)     !-- dU=dU+dFz/dz
       call evolve_src(irk)
       glm_ch_old=maxval(aflux)
    else if(irk==3) then
       u=u0-dt*(du(:,:,:,:,2)+(1.d0-alpha_imex)*du_stiff(:,:,:,:,2))
       call evolve_u2v
       glm_ch=glm_ch_old
       call evolve_dfx(irk)     !-- dU=dU+dFx/dx
       call evolve_dfy(irk)     !-- dU=dU+dFy/dy
       call evolve_dfz(irk)     !-- dU=dU+dFz/dz
       call evolve_src(irk)
       glm_ch_old=maxval(aflux)
    else if(irk==4) then 
       u=u0-dt*(du(:,:,:,:,2)/4.d0+du(:,:,:,:,3)/4.d0+beta_imex*du_stiff(:,:,:,:,1)+eta_imex*du_stiff(:,:,:,:,2)+&
                     &delta_imex*du_stiff(:,:,:,:,3))
       call evolve_u2v
       glm_ch=glm_ch_old
       call evolve_dfx(irk)     !-- dU=dU+dFx/dx
       call evolve_dfy(irk)     !-- dU=dU+dFy/dy
       call evolve_dfz(irk)     !-- dU=dU+dFz/dz
       call evolve_src(irk)
       glm_ch_old=maxval(aflux)
       u=u0-dt*(du(:,:,:,:,2)/6.d0+du(:,:,:,:,3)/6.d0+du(:,:,:,:,4)*2.d0/3.d0+du_stiff(:,:,:,:,2)/6.d0+&
        &du_stiff(:,:,:,:,3)/6.d0+du_stiff(:,:,:,:,4)*2.d0/3.d0)
       if(divclean) call  evolve_glm_source(dt/2.)
       call evolve_u2v
    end if  
!DDD    write(*,*) "u(kpr) at the end of the step is", u(2:3,2:3,:,kpr)
  end if

  if(prl) then
     not_crash_flag=(.not. run_crashed)
     call MPI_Allreduce(not_crash_flag,check_crash,1,MPI_Logical,MPI_LAND,MPI_COMM_WORLD,errmpi)
     if(.not. check_crash) then
        run_crashed=.true.
        return
     end if
  else !serial run case
     if(run_crashed) return
  end if

  if(viscosity .and. enable_smooth_viscosity) then
     call evolve_smooth_viscosity
     !if(viscosity .eq. 1) call evolve_smooth_viscosity_2013
     if(prl) then
        not_crash_flag=(.not. run_crashed)
        call MPI_Allreduce(not_crash_flag,check_crash,1,MPI_Logical,MPI_LAND,MPI_COMM_WORLD,errmpi)
        if(.not. check_crash) then
           run_crashed=.true.
           return
        end if
     else !serial run case
        if(run_crashed) return
     end if
  end if
 

end subroutine evolve_main

! *****************************************************************************

subroutine evolve_dfx(iteration)

  use common
  use work

  integer :: ix,iy,iz,k
  real(8) :: v2

  integer :: error_flag
  integer :: errmpi
  logical :: not_crash_flag
  integer :: iteration

  error_flag=0

  if (nx>1) then

    w(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)
    call evolve_bcx(ix1,iy1,iz1,1,nv,kv(1:nv))
    do iz=iz1,iz2
    do iy=iy1,iy2
      w0(1-ngc:mx+ngc,1:nv)=w(ix1-ngc:ix2+ngc,iy,iz,1:nv)
      wd(1-ngc:mx+ngc,1)=1./sqrt(1.-w(ix1-ngc:ix2+ngc,iy,iz,kvx)**2.-w(ix1-ngc:ix2+ngc,iy,iz,kvy)**2.-g_cov(3)*w(ix1-ngc:ix2+&
                         &ngc,iy,iz,kvz)**2.)
      wd(1-ngc:mx+ngc,2)=wd(1-ngc:mx+ngc,1)*w(ix1-ngc:ix2+ngc,iy,iz,kvx)
      wd(1-ngc:mx+ngc,3)=wd(1-ngc:mx+ngc,1)*w(ix1-ngc:ix2+ngc,iy,iz,kvy)
      wd(1-ngc:mx+ngc,4)=wd(1-ngc:mx+ngc,1)*w(ix1-ngc:ix2+ngc,iy,iz,kvz)

      call work_rec(mx,1,nv,0)
      call work_rec(mx,1,4,1)

      
      do ix=0,mx
        if(wr(ix,kpr) .le. 0) then
          write(*,*) "Evolve dfx: pr<=0 in wr:", wr(ix,kpr),"at ix,iy,iz:",ix,iy,iz
          wr(ix,kpr)=w0(ix,kpr)
          write(*,*) "After fix now pr is:", wr(ix,kpr)
        end if
        v2=sum(wr(ix,kvx:kvz)**2*g_cov(1:3))
        if (v2>=1.) then
           write(*,*) "Evolve dfx: v2>0 in wr:",v2, "at ix,iy,iz:",ix,iy,iz
           wr(ix,kvx:kvz)=w0(ix,kvx:kvz)
           v2=sum(wr(ix,kvx:kvz)**2*g_cov(1:3))
           write(*,*) "After fix now v2 is:", v2
           if(v2>=1) then
             write(*,*) "New reduction of v2 to:",maxspeed
             w0(ix,kvx:kvz)=w0(ix,kvx:kvz)*maxspeed/(sqrt(v2))
             wr(ix,kvx:kvz)=w0(ix,kvx:kvz)
           end if
        end if
        if(wl(ix,kpr) .le. 0) then
          write(*,*) "Evolve dfx: pr<=0 in wr:", wl(ix,kpr),"at ix,iy,iz:",ix,iy,iz
          wl(ix,kpr)=w0(ix,kpr)
          write(*,*) "After fix now pr is:", wl(ix,kpr)
        end if
        v2=sum(wl(ix,kvx:kvz)**2*g_cov(1:3))
        if (v2>=1.) then
           write(*,*) "Evolve dfx: v2>0 in wl:",v2, "at ix,iy,iz:",ix,iy,iz
           wl(ix,kvx:kvz)=w0(ix,kvx:kvz)
           v2=sum(wl(ix,kvx:kvz)**2*g_cov(1:3))
           write(*,*) "After fix now v2 is:", v2
           if(v2>=1) then
             write(*,*) "New reduction of v2 to:",maxspeed
             w0(ix,kvx:kvz)=w0(ix,kvx:kvz)*maxspeed/(sqrt(v2))
             wl(ix,kvx:kvz)=w0(ix,kvx:kvz)
           end if
        end if
      end do

      deriv(ix1:ix2,iy,iz,dtx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,1) + wdr(1:mx,1) - wdl(0:mx-1,1) - wdr(0:mx-1,1))
      deriv(ix1:ix2,iy,iz,dxx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,2) + wdr(1:mx,2) - wdl(0:mx-1,2) - wdr(0:mx-1,2))
      deriv(ix1:ix2,iy,iz,dyx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,3) + wdr(1:mx,3) - wdl(0:mx-1,3) - wdr(0:mx-1,3))
      deriv(ix1:ix2,iy,iz,dzx)=ddx(ix1:ix2)*0.5*( wdl(1:mx,4) + wdr(1:mx,4) - wdl(0:mx-1,4) - wdr(0:mx-1,4))

      if(divclean) call work_glmx(mx)
      call work_fupx(mx,x0,y(iy),z(iz),error_flag)

      if (ho) then
        w(ix0:ix2,iy,iz,1:nu)=w1(0:mx,1:nu)
      else
        do k=1,nu
           du(ix1:ix2,iy,iz,k,iteration)=du(ix1:ix2,iy,iz,k,iteration) + ddx(ix1:ix2)*(w1(1:mx,k)-w1(0:mx-1,k))
        end do
        !DDD  write(*,*) du(ix1:ix2,iy,iz,kby)
      end if
      ax(ix0:ix2,iy,iz,1:2)=w1(0:mx,nv+1:nv+2)
    end do
    end do

    if(error_flag .gt. 0) then
       write(*,*) "Error when computing fluxes along x"
       errcode=error_flag
       write(*,*) "Error code: ", errcode
       run_crashed=.true.
    end if

    if(prl) then
       not_crash_flag=(.not. run_crashed)
       call MPI_Allreduce(not_crash_flag,check_crash,1,MPI_Logical,MPI_LAND,MPI_COMM_WORLD,errmpi)
       if(.not. check_crash) then
          run_crashed=.true.
          return
       end if
    else !serial run case
       if(run_crashed) return
    end if
    
    if (ho) then
      call evolve_bcx(ix0,iy1,iz1,1,nu,ku(1:nu))
      do iz=iz1,iz2
      do iy=iy1,iy2
        w0(1-ngc:mx+ngc,1:nu)=w(ix1-ngc:ix2+ngc,iy,iz,1:nu)
        call work_der(mx,1,nu,0)
        do k=1,nu
          du(ix1:ix2,iy,iz,k,iteration)=du(ix1:ix2,iy,iz,k,iteration)+ddx(ix1:ix2)*w1(1:mx,k)
        end do
      end do
      end do
    end if 
   !!!!DDD for debugging only, it can be safely removed
   !write(*,*) "Time: ",t
   !do ix=ix1,ix2
   !   write(*,*) x(ix),v(ix,1,1,kvx), v(ix,1,1,kby), du(ix,1,1,kby)
   !end do
  end if

end subroutine evolve_dfx

! *****************************************************************************

subroutine evolve_dfy(iteration)

  use common
  use work

  integer :: ix,iy,iz,k
  real(8) :: v2

  integer :: error_flag
  integer :: errmpi

  logical :: not_crash_flag
  integer :: iteration

  error_flag=0

  if (ny>1) then

    w(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)
    call evolve_bcy(ix1,iy1,iz1,1,nv,kv(1:nv))
    do iz=iz1,iz2
    do ix=ix1,ix2
      w0(1-ngc:my+ngc,1:nv)=w(ix,iy1-ngc:iy2+ngc,iz,1:nv)
      wd(1-ngc:my+ngc,1)=1./sqrt(1.-w(ix,iy1-ngc:iy2+ngc,iz,kvx)**2.-w(ix,iy1-ngc:iy2+ngc,iz,kvy)**2.-g_cov(3)*w(ix,iy1-ngc:iy2+&
                         &ngc,iz,kvz)**2.)
      wd(1-ngc:my+ngc,2)=wd(1-ngc:my+ngc,1)*w(ix,iy1-ngc:iy2+ngc,iz,kvx)
      wd(1-ngc:my+ngc,3)=wd(1-ngc:my+ngc,1)*w(ix,iy1-ngc:iy2+ngc,iz,kvy)
      wd(1-ngc:my+ngc,4)=wd(1-ngc:my+ngc,1)*w(ix,iy1-ngc:iy2+ngc,iz,kvz)
      call work_rec(my,1,nv,0)
      call work_rec(my,1,4,1)

      do iy=0,my
        if(wr(iy,kpr) .le. 0) then
          write(*,*) "Evolve dfy: pr<=0 in wr:", wr(iy,kpr),"at ix,iy,iz:",ix,iy,iz
          wr(iy,kpr)=w0(iy,kpr)
          write(*,*) "After fix now pr is:", wr(iy,kpr)
        end if
        v2=sum(wr(iy,kvx:kvz)**2*g_cov(1:3))
        if (v2>=1.) then
           write(*,*) "Evolve dfy: v2>0 in wr:",v2, "at ix,iy,iz:",ix,iy,iz
           wr(iy,kvx:kvz)=w0(iy,kvx:kvz)
           v2=sum(wr(iy,kvx:kvz)**2*g_cov(1:3))
           write(*,*) "After fix now v2 is:", v2
           if(v2>=1) then
             write(*,*) "New reduction of v2 to:",maxspeed
             w0(iy,kvx:kvz)=w0(iy,kvx:kvz)*maxspeed/(sqrt(v2))
             wr(iy,kvx:kvz)=w0(iy,kvx:kvz)
           end if
        end if
        if(wl(iy,kpr) .le. 0) then
          write(*,*) "Evolve dfy: pr<=0 in wl:", wl(iy,kpr),"at ix,iy,iz:",ix,iy,iz
          wl(iy,kpr)=w0(iy,kpr)
          write(*,*) "After fix now pr is:", wl(iy,kpr)
        end if
        v2=sum(wl(iy,kvx:kvz)**2*g_cov(1:3))
        if (v2>=1.) then
           write(*,*) "Evolve dfy: v2>0 in wl:",v2, "at ix,iy,iz:",ix,iy,iz
           wl(iy,kvx:kvz)=w0(iy,kvx:kvz)
           v2=sum(wl(iy,kvx:kvz)**2*g_cov(1:3))
           write(*,*) "After fix now v2 is:", v2
           if(v2>=1) then
             write(*,*) "New reduction of v2 to:",maxspeed
             w0(iy,kvx:kvz)=w0(iy,kvx:kvz)*maxspeed/(sqrt(v2))
             wl(iy,kvx:kvz)=w0(iy,kvx:kvz)
           end if
        end if
      end do

      deriv(ix,iy1:iy2,iz,dty)=ddy(iy1:iy2)*0.5*( wdl(1:my,1) + wdr(1:my,1) - wdl(0:my-1,1) - wdr(0:my-1,1))
      deriv(ix,iy1:iy2,iz,dxy)=ddy(iy1:iy2)*0.5*( wdl(1:my,2) + wdr(1:my,2) - wdl(0:my-1,2) - wdr(0:my-1,2))
      deriv(ix,iy1:iy2,iz,dyy)=ddy(iy1:iy2)*0.5*( wdl(1:my,3) + wdr(1:my,3) - wdl(0:my-1,3) - wdr(0:my-1,3))
      deriv(ix,iy1:iy2,iz,dzy)=ddy(iy1:iy2)*0.5*( wdl(1:my,4) + wdr(1:my,4) - wdl(0:my-1,4) - wdr(0:my-1,4))

      if(divclean) call work_glmy(my)
      call work_fupy(my,x(ix),y0,z(iz),error_flag)

      if (ho) then
        w(ix,iy0:iy2,iz,1:nu)=w1(0:my,1:nu)
      else
        do k=1,nu
          du(ix,iy1:iy2,iz,k,iteration)=du(ix,iy1:iy2,iz,k,iteration) + &
            ddy(iy1:iy2)*(w1(1:my,k)-w1(0:my-1,k))
        end do
      end if
      ay(ix,iy0:iy2,iz,1:2)=w1(0:my,nv+1:nv+2)
      if(error_flag .gt. 0) exit
    end do
    if(error_flag .gt. 0) exit
    end do

    if(error_flag .gt. 0) then
      write(*,*) "Error when computing fluxes along y"
      errcode=error_flag
      write(*,*) "Error code: ", errcode
      run_crashed=.true.
    end if
    
    if(prl) then
       not_crash_flag=(.not. run_crashed)
       call MPI_Allreduce(not_crash_flag,check_crash,1,MPI_Logical,MPI_LAND,MPI_COMM_WORLD,errmpi)
       if(.not. check_crash) then
          run_crashed=.true.
          return
       end if
    else !serial run case
       if(run_crashed) return
    end if

    if (ho) then
      call evolve_bcy(ix1,iy0,iz1,1,nu,ku(1:nu))
      do iz=iz1,iz2
      do ix=ix1,ix2
        w0(1-ngc:my+ngc,1:nu)=w(ix,iy1-ngc:iy2+ngc,iz,1:nu)
        call work_der(my,1,nu,0)
        do k=1,nu
          du(ix,iy1:iy2,iz,k,iteration)=du(ix,iy1:iy2,iz,k,iteration)+ddy(iy1:iy2)*w1(1:my,k)
        end do
      end do
      end do
    end if

  end if

end subroutine evolve_dfy

! *****************************************************************************

subroutine evolve_dfz(iteration)

  use common
  use work

  integer :: ix,iy,iz,k
  real(8) :: v2

  integer error_flag
  integer :: errmpi
  logical :: not_crash_flag
  integer :: iteration

  error_flag=0

  if (nz>1) then

    w(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)=v(ix1:ix2,iy1:iy2,iz1:iz2,1:nv)
    call evolve_bcz(ix1,iy1,iz1,1,nv,kv(1:nv))
    do iy=iy1,iy2
    do ix=ix1,ix2
      w0(1-ngc:mz+ngc,1:nv)=w(ix,iy,iz1-ngc:iz2+ngc,1:nv)
      wd(1-ngc:mz+ngc,1)=1./sqrt(1.-w(ix,iy,iz1-ngc:iz2+ngc,kvx)**2.-w(ix,iy,iz1-ngc:iz2+ngc,kvy)**2.-g_cov(3)*w(ix,iy,iz1-ngc:&
                         &iz2+ngc,kvz)**2.)
      wd(1-ngc:mz+ngc,2)=wd(1-ngc:mz+ngc,1)*w(ix,iy,iz1-ngc:iz2+ngc,kvx)
      wd(1-ngc:mz+ngc,3)=wd(1-ngc:mz+ngc,1)*w(ix,iy,iz1-ngc:iz2+ngc,kvy)
      wd(1-ngc:mz+ngc,4)=wd(1-ngc:mz+ngc,1)*w(ix,iy,iz1-ngc:iz2+ngc,kvz)
      call work_rec(mz,1,nv,0)
      call work_rec(mz,1,4,1)

      do iz=0,mz
        if(wr(iz,kpr) .le. 0) then
          write(*,*) "Evolve dfz: pr<=0 in wr:", wr(iz,kpr),"at ix,iy,iz:",ix,iy,iz
          wr(iz,kpr)=w0(iz,kpr)
          write(*,*) "After fix now pr is:", wr(iz,kpr)
        end if
        v2=sum(wr(iz,kvx:kvz)**2*g_cov(1:3))
        if (v2>=1.) then
           write(*,*) "Evolve dfz: v2>0 in wr:",v2, "at ix,iy,iz:",ix,iy,iz
           wr(iz,kvx:kvz)=w0(iz,kvx:kvz)
           v2=sum(wr(iz,kvx:kvz)**2*g_cov(1:3))
           write(*,*) "After fix now v2 is:", v2
           if(v2>=1) then
             write(*,*) "New reduction of v2 to:",maxspeed
             w0(iz,kvx:kvz)=w0(iz,kvx:kvz)*maxspeed/(sqrt(v2))
             wr(iz,kvx:kvz)=w0(iz,kvx:kvz)
           end if
        end if
        if(wl(iz,kpr) .le. 0) then
          write(*,*) "Evolve dfz: pr<=0 in wl:", wl(iz,kpr),"at ix,iy,iz:",ix,iy,iz
          wl(iz,kpr)=w0(iz,kpr)
          write(*,*) "After fix now pr is:", wl(iz,kpr)
        end if
        v2=sum(wl(iz,kvx:kvz)**2*g_cov(1:3))
        if (v2>=1.) then
           write(*,*) "Evolve dfz: v2>0 in wl:",v2, "at ix,iy,iz:",ix,iy,iz
           wl(iz,kvx:kvz)=w0(iz,kvx:kvz)
           v2=sum(wl(iz,kvx:kvz)**2*g_cov(1:3))
           write(*,*) "After fix now v2 is:", v2
           if(v2>=1) then
             write(*,*) "New reduction of v2 to:",maxspeed
             w0(iz,kvx:kvz)=w0(iz,kvx:kvz)*maxspeed/(sqrt(v2))
             wl(iz,kvx:kvz)=w0(iz,kvx:kvz)
           end if
        end if
      end do

      deriv(ix,iy,iz1:iz2,dtz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,1) + wdr(1:mz,1) - wdl(0:mz-1,1) - wdr(0:mz-1,1))
      deriv(ix,iy,iz1:iz2,dxz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,2) + wdr(1:mz,2) - wdl(0:mz-1,2) - wdr(0:mz-1,2))
      deriv(ix,iy,iz1:iz2,dyz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,3) + wdr(1:mz,3) - wdl(0:mz-1,3) - wdr(0:mz-1,3))
      deriv(ix,iy,iz1:iz2,dzz)=ddz(iz1:iz2)*0.5*( wdl(1:mz,4) + wdr(1:mz,4) - wdl(0:mz-1,4) - wdr(0:mz-1,4))

      if(divclean) call work_glmz(mz)
      call work_fupz(mz,x(ix),y(iy),z0,error_flag)

      if (ho) then
        w(ix,iy,iz0:iz2,1:nu)=w1(0:mz,1:nu)
      else
        do k=1,nu
          du(ix,iy,iz1:iz2,k,iteration)=du(ix,iy,iz1:iz2,k,iteration) + &
            ddz(iz1:iz2)*(w1(1:mz,k)-w1(0:mz-1,k))
        end do
      end if
      az(ix,iy,iz0:iz2,1:2)=w1(0:mz,nv+1:nv+2)
    end do
    end do

    if(error_flag .gt. 0) then
       write(*,*) "Error when computing fluxes along z"
       errcode=error_flag
       write(*,*) "Error code: ", errcode
       run_crashed=.true.
    end if
    
    if(prl) then
       not_crash_flag=(.not. run_crashed)
       call MPI_Allreduce(not_crash_flag,check_crash,1,MPI_Logical,MPI_LAND,MPI_COMM_WORLD,errmpi)
       if(.not. check_crash) then
          run_crashed=.true.
          return
       end if
    else !serial run case
       if(run_crashed) return
    end if

    if (ho) then
      call evolve_bcz(ix1,iy1,iz0,1,nu,ku(1:nu))
      do iy=iy1,iy2
      do ix=ix1,ix2
        w0(1-ngc:mz+ngc,1:nu)=w(ix,iy,iz1-ngc:iz2+ngc,1:nu)
        call work_der(mz,1,nu,0)
        do k=1,nu
          du(ix,iy,iz1:iz2,k,iteration)=du(ix,iy,iz1:iz2,k,iteration)+ddz(iz1:iz2)*w1(1:mz,k)
        end do
      end do
      end do
    end if
    if(error_flag .gt. 0) then
      write(*,*) "Error when computing fluxes along z"
      errcode=error_flag
      write(*,*) "Error code: ", errcode
      run_crashed=.true.
    end if

  end if

end subroutine evolve_dfz

! *****************************************************************************

subroutine evolve_src(iteration)

  use common
  use system_eqgp

  real(8),dimension(nv) :: vv
  real(8),dimension(nv) :: src,src_stiff
  integer :: error_flag
  integer :: ix,iy,iz
  integer :: iteration

  error_flag=0

  call dts(vnewstep,vold)
  do iz=iz1,iz2
  do iy=iy1,iy2
  do ix=ix1,ix2
    vv=v(ix,iy,iz,1:nv)
    call system_source(ix,iy,iz,v,vv,src,src_stiff,error_flag)
    if(error_flag .gt. 0) then
      errcode=error_flag
      return
    end if
    du(ix,iy,iz,1:nu,iteration)=du(ix,iy,iz,1:nu,iteration)-src(ku(1:nu))
    du_stiff(ix,iy,iz,1:nu,iteration)=du_stiff(ix,iy,iz,1:nu,iteration)-src_stiff(ku(1:nu))
  end do
  end do
  end do

end subroutine evolve_src


! *****************************************************************************

subroutine evolve_dt(dt)

  use parallel
  use common
  use viscous

  real(8),intent(inout) :: dt

  integer :: ix,iy,iz
  real(8)    :: amax

    if(divclean) then
     amax=1
    else
     amax=0.

     do iz=iz1,iz2
      do iy=iy1,iy2
       do ix=ix1,ix2
        if (nx>1) amax=max(amax,ddx(ix)*maxval(ax(ix-1:ix,iy,iz,1:2)))
        if (ny>1) amax=max(amax,ddy(iy)*maxval(ay(ix,iy-1:iy,iz,1:2)))
        if (nz>1) amax=max(amax,ddz(iz)*maxval(az(ix,iy,iz-1:iz,1:2)))
       end do
      end do
     end do

     if(amax .eq. 0) then
      do iz=iz1,iz2
       do iy=iy1,iy2
        do ix=ix1,ix2
         if (nx>1) amax=max(amax,ddx(ix))
         if (ny>1) amax=max(amax,ddy(iy))
         if (nz>1) amax=max(amax,ddz(iz))
        end do
       end do
      end do
     end if
    end if !end if divclean
    
    dt=cfl/amax

end subroutine evolve_dt


! *****************************************************************************

subroutine evolve_u2v

  use common
  use system_eqgp
  use viscous
  implicit none
  real(8),dimension(nv) :: uu,vv

  integer :: ix,iy,iz
  
  real(8) :: dutdt, dutdx, dutdy, dutdz, duxdx, duydy, duzdz, duzdt, duzdx, duzdy, duxdt,&
                      &duxdy, duxdz, duydt, duydx, duydz

  integer :: error_flag
 
  real(8) pxl, pxr, pyl, pyr, pzl, pzr, v2
  real(8),parameter :: prlimit=0.01 !it issues a warning when fixing problems and pressure is above this threshold

  error_flag=0
  
  do iz=iz1,iz2
  do iy=iy1,iy2
  do ix=ix1,ix2
    vv(1:nv)=v(ix,iy,iz,1:nv)
    uu(1:nu)=u(ix,iy,iz,1:nu)
    call system_cons2prim(uu,vv,ierr, error_flag) 
    if(error_flag .ne. 0) then
      if(nx>0 .and. ny>0 .and. nz>0) then
        error_flag=0
        run_crashed=.false.
        write(*,*) "Error detected in, trying to fix"
        write(*,*) "Position (index and coord) is:", ix, iy, iz, x(ix), y(iy), z(iz)
        u(ix,iy,iz,1:nu)=(u(ix-1,iy,iz,1:nu)+u(ix+1,iy,iz,1:nu)+u(ix,iy-1,iz,1:nu)+u(ix,iy+1,iz,1:nu)+u(ix,iy,iz-1,1:nu)+&
                        &u(ix,iy,iz+1,1:nu))/6.
        uu(1:nu)=u(ix,iy,iz,1:nu)
        call system_cons2prim(uu,vv,ierr, error_flag)
        v2=vv(kvx)**2+vv(kvy)**2+vv(kvz)**2*g_cov(3)
        if(v2 .gt. 1) then
          vv(kvx:kvz)=maxspeed/sqrt(v2)*vv(kvx:kvz)
          write(*,*) "Superluminal speed, I reduced it to",maxspeed," from",sqrt(v2)
        end if
        error_flag=0
        run_crashed=.false.

        if(vv(kpr) .lt. 0) then
          write(*,*) "Pressure is negative, trying to fix"
          pxl=v(ix-1,iy,iz,kpr)
          pxr=v(ix+1,iy,iz,kpr)
          pyl=v(ix,iy-1,iz,kpr)
          pyr=v(ix,iy+1,iz,kpr)
          pzl=v(ix,iy,iz-1,kpr)
          pzr=v(ix,iy,iz+1,kpr)
          if((pxl<prlimit) .and. (pxr<prlimit) .and. (pyl<prlimit) .and. (pyr<prlimit) .and. (pzl<prlimit) .and. (pzr<prlimit)) then
            vv(kpr)=(pxl+pxr+pyl+pyr+pzl+pzr)/6.
          else
            vv(kpr)=(pxl+pxr+pyl+pyr+pzl+pzr)/6.
            write(*,*) "Warning, I used the average value of pressure in the neighbouring cells, but in some of them is > prlimit!"
          end if
        end if
        write(*,*) "Executing prim2cons"
        call system_prim2cons(vv,uu,error_flag)
        call system_cons2prim(uu,vv,ierr, error_flag)
        write(*,*) "Executing cons2prim and checking the results"
        v2=vv(kvx)**2+vv(kvy)**2+vv(kvz)**2*g_cov(3)
        if((vv(kpr)>0) .and. (v2<1)) then
             error_flag=0
             run_crashed=.false.
        else
             error_flag=1
             run_crashed=.true.
             write(*,*) "Sorry, I tried to fix the error, but I failed..."
        end if 
      end if
    end if

    if(viscosity .and. (obtained .eq. 'zz')) then
      v(ix,iy,iz,kpizz)=vv(kpizz)
    end if

    u(ix,iy,iz,1:nu)=uu(ku(1:nu)) !because in the cons2prim sometimes the conserved values are changed
    if (error_flag>0) then
      print*,'Problem with cons2prim pointed out by evolve_u2v'
      print*,'time:',t
      errcode=error_flag
      print*,'iteration or error code = ',errcode 
      print*,'ix,iy,iz = ',ix,iy,iz
      print*,'x,y,z = ',x(ix),y(iy),z(iz)
      print*,'primitives = ',vv
      run_crashed=.true.
      return
    end if

     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)
      
  end do
  end do
  end do

end subroutine evolve_u2v

! *****************************************************************************

subroutine evolve_bcx(ixl,iyl,izl,kl,kr,kbc)
!  0 = periodic
!  1 = reflecting (positive)
! -1 = reflecting (negative)
!  2 = zeroth order extrapolation
!  3 = third order extrapolation

  use parallel
  use common

  integer,intent(in) :: ixl,iyl,izl,kl,kr,kbc(kl:kr)
  integer            :: i0,l0,l,k,lenw

  if (ixl==ix1) then
    i0=1
    l0=ngc
  else
    i0=0
    l0=ngc-1
  end if

  if (prl) then ! Circular data transfer if parallel

    allocate(wsend(l0,iyl:iy2,izl:iz2,kl:kr))
    allocate(wrecv(l0,iyl:iy2,izl:iz2,kl:kr))
    lenw=size(wsend)

    wsend=w(ix2-ngc+1:ix2-ngc+l0,iyl:iy2,izl:iz2,kl:kr)
    call mpi_sendrecv(wsend,lenw,mpi_realtype,ipe_next,1, &
                      wrecv,lenw,mpi_realtype,ipe_prev,1, &
                      icomm,istat,ierr)
    w(ix0-ngc+1:ix0-ngc+l0,iyl:iy2,izl:iz2,kl:kr)=wrecv

    wsend=w(ix0+1:ix0+l0,iyl:iy2,izl:iz2,kl:kr)
    call mpi_sendrecv(wsend,lenw,mpi_realtype,ipe_prev,2, &
                      wrecv,lenw,mpi_realtype,ipe_next,2, &
                      icomm,istat,ierr)
    w(ix2+1:ix2+l0,iyl:iy2,izl:iz2,kl:kr)=wrecv

    deallocate(wsend,wrecv)

  end if
 
  do k=kl,kr

    if (ipe==0) then ! Only the first processor

    select case(ibc(kbc(k),1,1))
    case(0)
      do l=1,l0
        w(-ngc+l,iyl:iy2,izl:iz2,k)=w(nx-ngc+l,iyl:iy2,izl:iz2,k)
      end do
    case(1)
      do l=1,l0
        w(-ngc+l,iyl:iy2,izl:iz2,k)=w(ngc+i0-l,iyl:iy2,izl:iz2,k)
      end do
    case(-1)
      do l=1,l0
        w(-ngc+l,iyl:iy2,izl:iz2,k)=-w(ngc+i0-l,iyl:iy2,izl:iz2,k)
      end do
    case(2)
      do l=1,l0
        w(-ngc+l,iyl:iy2,izl:iz2,k)=w(i0,iyl:iy2,izl:iz2,k)
      end do
    case(3)
      do l=1,l0
        w(i0-l,iyl:iy2,izl:iz2,k)=w(i0-l+3,iyl:iy2,izl:iz2,k) &
          -3.*(w(i0-l+2,iyl:iy2,izl:iz2,k)-w(i0-l+1,iyl:iy2,izl:iz2,k))
      end do
    end select

    end if

    if (ipe==npe-1) then ! Only the last processor

    select case(ibc(kbc(k),2,1))
    case(0)
      do l=1,l0
        w(nx+l,iyl:iy2,izl:iz2,k)=w(l,iyl:iy2,izl:iz2,k)
      end do
    case(1)
      do l=1,l0
        w(nx+l,iyl:iy2,izl:iz2,k)=w(nx+i0-l,iyl:iy2,izl:iz2,k)
      end do
    case(-1)
      do l=1,l0
        w(nx+l,iyl:iy2,izl:iz2,k)=-w(nx+i0-l,iyl:iy2,izl:iz2,k)
      end do
    case(2)
      do l=1,l0
        w(nx+l,iyl:iy2,izl:iz2,k)=w(nx,iyl:iy2,izl:iz2,k)
      end do
    case(3)
      do l=1,l0
        w(nx+l,iyl:iy2,izl:iz2,k)=w(nx+l-3,iyl:iy2,izl:iz2,k) &
          -3.*(w(nx+l-2,iyl:iy2,izl:iz2,k)-w(nx+l-1,iyl:iy2,izl:iz2,k))
      end do
    end select

    end if

  end do

end subroutine evolve_bcx

! *****************************************************************************

subroutine evolve_bcy(ixl,iyl,izl,kl,kr,kbc)
!  0 = periodic
!  1 = reflecting (positive)
! -1 = reflecting (negative)
!  2 = zeroth order extrapolation
!  3 = third order extrapolation

  use common

  integer,intent(in) :: ixl,iyl,izl,kl,kr,kbc(kl:kr)
  integer            :: i0,l0,l,k,lenw

  if (iyl==iy1) then
    i0=1
    l0=ngc
  else
    i0=0
    l0=ngc-1
  end if

  do k=kl,kr

    select case(ibc(kbc(k),1,2))
    case(0)
      do l=1,l0
        w(ixl:ix2,-ngc+l,izl:iz2,k)=w(ixl:ix2,ny-ngc+l,izl:iz2,k)
      end do
    case(1)
      do l=1,l0
        w(ixl:ix2,-ngc+l,izl:iz2,k)=w(ixl:ix2,ngc+i0-l,izl:iz2,k)
      end do
    case(-1)
      do l=1,l0
        w(ixl:ix2,-ngc+l,izl:iz2,k)=-w(ixl:ix2,ngc+i0-l,izl:iz2,k)
      end do
    case(2)
      do l=1,l0
        w(ixl:ix2,-ngc+l,izl:iz2,k)=w(ixl:ix2,i0,izl:iz2,k)
      end do
    case(3)
      do l=1,l0
        w(ixl:ix2,i0-l,izl:iz2,k)=w(ixl:ix2,i0-l+3,izl:iz2,k) &
          -3.*(w(ixl:ix2,i0-l+2,izl:iz2,k)-w(ixl:ix2,i0-l+1,izl:iz2,k))
      end do
    end select

    select case(ibc(kbc(k),2,2))
    case(0)
      do l=1,l0
        w(ixl:ix2,ny+l,izl:iz2,k)=w(ixl:ix2,l,izl:iz2,k)
      end do
    case(1)
      do l=1,l0
        w(ixl:ix2,ny+l,izl:iz2,k)=w(ixl:ix2,ny+i0-l,izl:iz2,k)
      end do
    case(-1)
      do l=1,l0
        w(ixl:ix2,ny+l,izl:iz2,k)=-w(ixl:ix2,ny+i0-l,izl:iz2,k)
      end do
    case(2)
      do l=1,l0
        w(ixl:ix2,ny+l,izl:iz2,k)=w(ixl:ix2,ny,izl:iz2,k)
      end do
    case(3)
      do l=1,l0
        w(ixl:ix2,ny+l,izl:iz2,k)=w(ixl:ix2,ny+l-3,izl:iz2,k) &
          -3.*(w(ixl:ix2,ny+l-2,izl:iz2,k)-w(ixl:ix2,ny+l-1,izl:iz2,k))
      end do
    end select

  end do

end subroutine evolve_bcy

! *****************************************************************************

subroutine evolve_bcz(ixl,iyl,izl,kl,kr,kbc)
!  0 = periodic
!  1 = reflecting (positive)
! -1 = reflecting (negative)
!  2 = zeroth order extrapolation
!  3 = third order extrapolation

  use common

  integer,intent(in) :: ixl,iyl,izl,kl,kr,kbc(kl:kr)
  integer            :: i0,l0,l,k,lenw

  if (izl==iz1) then
    i0=1
    l0=ngc
  else
    i0=0
    l0=ngc-1
  end if

  do k=kl,kr

    select case(ibc(kbc(k),1,3))
    case(0)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,-ngc+l,k)=w(ixl:ix2,iyl:iy2,nz-ngc+l,k)
      end do
    case(1)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,-ngc+l,k)=w(ixl:ix2,iyl:iy2,ngc+i0-l,k)
      end do
    case(-1)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,-ngc+l,k)=-w(ixl:ix2,iyl:iy2,ngc+i0-l,k)
      end do
    case(2)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,-ngc+l,k)=w(ixl:ix2,iyl:iy2,i0,k)
      end do
    case(3)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,i0-l,k)=w(ixl:ix2,iyl:iy2,i0-l+3,k) &
          -3.*(w(ixl:ix2,iyl:iy2,i0-l+2,k)-w(ixl:ix2,iyl:iy2,i0-l+1,k))
      end do
    end select

    select case(ibc(kbc(k),2,3))
    case(0)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,nz+l,k)=w(ixl:ix2,iyl:iy2,l,k)
      end do
    case(1)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,nz+l,k)=w(ixl:ix2,iyl:iy2,nz+i0-l,k)
      end do
    case(-1)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,nz+l,k)=-w(ixl:ix2,iyl:iy2,nz+i0-l,k)
      end do
    case(2)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,nz+l,k)=w(ixl:ix2,iyl:iy2,nz,k)
      end do
    case(3)
      do l=1,l0
        w(ixl:ix2,iyl:iy2,nz+l,k)=w(ixl:ix2,iyl:iy2,nz+l-3,k) &
          -3.*(w(ixl:ix2,iyl:iy2,nz+l-2,k)-w(ixl:ix2,iyl:iy2,nz+l-1,k))
      end do
    end select

  end do

end subroutine evolve_bcz

! *****************************************************************************
subroutine evolve_check_xy_corners(zc)
use common
implicit none
real(8), dimension(1:nv) :: average
integer q !just a counter
integer, intent(in) :: zc
average(1:nv)=(u(ix1+1,iy1,zc,1:nv)+u(ix1,iy1+1,zc,1:nv)+u(ix1+1,iy1+1,zc,1:nv))/3.
do q=1,nv
   u(ix1,iy1,zc,q)=average(q)
end do

average(1:nv)=(u(ix2-1,iy1,zc,1:nv)+u(ix2,iy1+1,zc,1:nv)+u(ix2-1,iy1+1,zc,1:nv))/3.
do q=1,nv
   u(ix2,iy1,zc,q)=average(q)
end do

average(1:nv)=(u(ix1,iy2-1,zc,1:nv)+u(ix1+1,iy2,zc,1:nv)+u(ix1+1,iy2-1,zc,1:nv))/3.
do q=1,nv
   u(ix1,iy2,zc,q)=average(q)
end do

average(1:nv)=(u(ix2-1,iy2,zc,1:nv)+u(ix2,iy2-1,zc,1:nv)+u(ix2-1,iy2-1,zc,1:nv))/3.
do q=1,nv
   u(ix2,iy2,zc,q)=average(q)
end do

end subroutine evolve_check_xy_corners

! *****************************************************************************
subroutine evolve_smooth_viscosity
use common
use system_eqgp
implicit none
integer :: l,m,n
logical,save :: first_time=.true.
real(8) :: smooth_energy
real(8),save :: smooth_pressure=0.
real(8), parameter :: rh_fake=1.
integer :: error_flag

error_flag=0

if(first_time) then
   call find_pressure_from_temperature(smooth_temp, smooth_pressure)
   call eos_energy(rh_fake,smooth_energy,smooth_pressure,error_flag)
   if(error_flag .gt. 0) then
     write(*,*) "An error occurred when calling eos_energy inside evolve_smooth_viscosity (file evolve.f08)"
     errcode=error_flag
     write(*,*) "Error code:", errcode
     return
   end if
   if(pe0) write(*,*) "Pressure treshold (GeV/fm^3):", smooth_pressure,' - Energy density treshold (Gev/fm^3):',smooth_energy
   first_time=.false.
end if

do n=iz1,iz2
   do m=iy1,iy2
      do l=ix1,ix2
         if(v(l,m,n,kpr) .lt. smooth_pressure) then
            v(l,m,n,kpibu:kpizz)=v(l,m,n,kpibu:kpizz)*(v(l,m,n,kpr)/smooth_pressure)**2.
            if(obtained .eq. 'zz') then
              call get_derived_pi_zz(v(l,m,n,:),pitt,pitx,pity,pitz,pizz)
              v(l,m,n,kpizz)=pizz
            end if
            call system_prim2cons(v(l,m,n,:),u(l,m,n,:),error_flag)
            if(error_flag .gt. 0) then
              write(*,*) "An error occurred when trying to dim viscous tensor components at low energies."
              errcode=error_flag
              write(*,*) "Error code:", errcode
              write(*,*) "Position: ix=",l," iy=",m," iz=",n
              run_crashed=.true.
              return
            end if
         end if
      end do
   end do
end do
end subroutine evolve_smooth_viscosity
! *****************************************************************************
subroutine evolve_smooth_viscosity_2013
use common
use system_eqgp
use eos
implicit none
integer :: l,m,n
real(8) :: limitator, dens,press,energy_dens,temp,entropy_dens
integer :: error_flag

error_flag=0

do n=iz1,iz2
   do m=iy1,iy2
      do l=ix1,ix2
         press=v(l,m,n,kpr)
         if(press .lt. 1.e-2) then
            if(obtained .eq. 'no') then
              call get_derived_pi(v(l,m,n,1:nv),pitt,pitx,pity,pitz)
            else
              call get_derived_pi_zz(v(l,m,n,1:nv),pitt,pitx,pity,pitz,pizz)
              v(l,m,n,kpizz)=pizz
            end if
            call get_derived_data(dens, press, energy_dens, temp, entropy_dens,error_flag)
            if(error_flag .gt. 0) then
               write(*,*) "An error occurred when calling get_derived_data inside evolve_smooth_viscosity_2013 (file evolve.f08)"
               errcode=error_flag
               write(*,*) "Error code:", errcode
               write(*,*) "Position: ix=",l," iy=",m," iz=",n
               run_crashed=.true.
               return
            end if
            limitator=0.1*sqrt((v(l,m,n,kpixy)**2.+v(l,m,n,kpixx)**2.+v(l,m,n,kpiyy)**2.+pitt**2.+(v(l,m,n,kpizz)*g_cov(3))**2.+&
            &g_cov0*g_cov(3)*pitz**2.+g_cov0*pitx**2.+g_cov0*pity**2.+g_cov(3)*v(l,m,n,kpixz)**2.+&
            &g_cov(3)*v(l,m,n,kpiyz)**2.)/(3.*press**2.+energy_dens**2.))
            v(l,m,n,kpibu:kpizz)=v(l,m,n,kpibu:kpizz)*tanh(limitator)/limitator
            if(obtained .eq. 'zz') then
              call get_derived_pi_zz(v(l,m,n,1:nv),pitt,pitx,pity,pitz,pizz)
              v(l,m,n,kpizz)=pizz
            end if
             
            call system_prim2cons(v(l,m,n,:),u(l,m,n,:),error_flag)
            if(error_flag .gt. 0) then
              write(*,*) "An error occurred when trying to dim viscous tensor components at low energies."
              errcode=error_flag
              write(*,*) "Error code:", errcode
              write(*,*) "Position: ix=",l," iy=",m," iz=",n
              run_crashed=.true.
              return
            end if
         end if
      end do
   end do
end do
end subroutine evolve_smooth_viscosity_2013
! *****************************************************************************

end module evolve

! *****************************************************************************
