! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  File: viscous.f90                                                        *
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
! *  Authors: Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)            *
! *                                                                           *
! *  Contributors: Vinod Chandra (vchandra@iitgn.ac.in)                       *
! *                                                                           *
! *  Acknowledgments: Luca Del Zanna (delzanna@unifi.it)                      *
! *                                                                           *
! *****************************************************************************

module viscous
use parallel, only: ipe
use common, only: nv,krh,kvx,kvy,kvz,kpr,kpibu,kpizz,kpixz,kpixy,kpiyz,kpixx,kpiyy,g_cov,pitt,timeold,pitx,pity,pitz
use common, only: t, coordinates,nx,ny,nz, timeinterval, hbar,viscosity,g_cov0,init_type,g_cov_3_old,bulkvis
use common, only:  dtx, dty, dtz, dxx, dxy, dxz, dyx, dyy, dyz, dzx, dzy, dzz, deriv, dtt, dxt, dyt, dzt, dtet
use common, only: ix1,ix2,iy1,iy2,iz1,iz2,run_crashed,derivatives_out,eta_over_s,tau_pi_coeff
use common, only: GUBSER_VISCOUS, SHEAR_VISCOUS_1D, GLAUBER_GEOMETRIC, TABULATED_INIT, MINKOWSKI, BJORKEN, GLAUBER_MONTECARLO


implicit none
  integer, parameter :: nsigma=6 !number of sigma components needed to compute the viscosity tensor
  integer, parameter :: sxz=5, sxy=4, syz=6, sxx=1, syy=2, szz=3 !stt=7, stx=8, sty=9, stz=10
  real(8) :: tau_pi, eta_vis, zeta_vis, tau_pi_big
  real(8), dimension(1:nsigma) :: sigma

contains

! ***********************************************************************************

  ! surbroutine dts. derive vs time simple
  subroutine dts(quantity_array, quantity_array_old)
    use eos
    implicit none
    integer :: ix, iy, iz
    real(8), allocatable, dimension(:,:,:,:) :: quantity_array, quantity_array_old
    
    real(8) :: gamma2,gamma1 

    real(8) :: temperature_old, temperature_cur

    real(8) :: energy_bb, entropy_bb

    integer :: errcode !AAA DDD this error code should be passed as a return value to the caller subroutine

    if(timeinterval .eq. 0) then
      deriv(:,:,:,dtt:dzt)=0.
      return
    end if

    do iz=iz1,iz2
     do iy=iy1,iy2
      do ix=ix1,ix2
       gamma1=1./sqrt(1.-g_cov(1)*quantity_array_old(ix,iy,iz,kvx)**2.-g_cov(2)*quantity_array_old(ix,iy,iz,kvy)**2.-&
       &g_cov_3_old*quantity_array_old(ix,iy,iz,kvz)**2.)
       gamma2=1./sqrt(1.-g_cov(1)*quantity_array(ix,iy,iz,kvx)**2.-g_cov(2)*quantity_array(ix,iy,iz,kvy)**2.-&
       &g_cov(3)*quantity_array(ix,iy,iz,kvz)**2.)
      
       deriv(ix,iy,iz,dtt)=(gamma2-gamma1)/(timeinterval)

       deriv(ix,iy,iz,dxt:dzt)=(gamma2*quantity_array(ix,iy,iz,kvx:kvz)-gamma1*quantity_array_old(ix,iy,iz,kvx:kvz))/timeinterval
 
       call get_derived_data(quantity_array_old(ix,iy,iz,krh), quantity_array_old(ix,iy,iz,kpr), energy_bb, temperature_old, &
            &entropy_bb,errcode)
       call get_derived_data(quantity_array(ix,iy,iz,krh), quantity_array(ix,iy,iz,kpr), energy_bb, temperature_cur, entropy_bb,&
            &errcode)

       if(derivatives_out) deriv(ix,iy,iz,dtet)=(temperature_cur-temperature_old)/timeinterval
      end do       
     end do       
    end do       

  end subroutine dts

! ***********************************************************************************

  subroutine viscous_initio(ix,iy,iz,nv,v,tauzero,errcode)
    use eos
    implicit none
    integer,intent(in) :: ix,iy,iz, nv
    real(8) :: dens, press, vx, vy, vz, tauzero
    real(8) :: entropy_dens, energy_dens, temp
    real(8), dimension(:,:,:,:), allocatable, intent(inout) :: v
    integer :: errcode
    
    errcode=0
    dens=v(ix,iy,iz,krh)
    vx=v(ix,iy,iz,kvx)
    vy=v(ix,iy,iz,kvy)
    vz=v(ix,iy,iz,kvz)
    press=v(ix,iy,iz,kpr)

    call get_derived_data(dens, press, energy_dens, temp, entropy_dens,errcode)
    call get_viscous_parameters_0(eta_vis,zeta_vis,tau_pi, tau_pi_big, dens, press,errcode)
    if(errcode .gt. 0) then
      write(*,*) "An error occurred into the viscous_initio subroutine"
      write(*,*) "Error code:", errcode
      write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
      run_crashed=.true.
      return
    end if
    
    v(ix,iy,iz,kpibu:kpizz)=0.
    if((init_type .eq. GLAUBER_GEOMETRIC) .or. (init_type .eq. TABULATED_INIT) .or. (init_type .eq. GLAUBER_MONTECARLO)) then
      v(ix,iy,iz,kpizz)=-g_cov(3)*4.*eta_vis/(3.*tauzero)
      v(ix,iy,iz,kpixx)=2.*eta_vis/(3.*tauzero)
      v(ix,iy,iz,kpiyy)=2.*eta_vis/(3.*tauzero)
      v(ix,iy,iz,kpibu)=-zeta_vis/tauzero
    end if


  end subroutine viscous_initio

! ***********************************************************************************

  subroutine get_viscous_parameters_0(eta_vis,zeta_vis,tau_pi, tau_pi_big, dens, press,errcode)
    use eos
    implicit none
    real(8) :: entropy_dens, energy_dens, temp,cs2
    real(8), intent(in) :: dens, press
    real(8), intent(out) :: eta_vis, zeta_vis, tau_pi, tau_pi_big
    integer, intent(out) :: errcode

    errcode=0

    call get_derived_data(dens, press, energy_dens, temp, entropy_dens,errcode)
    call eos_sound(dens,energy_dens,press,cs2,errcode)

    eta_vis=hbar*eta_over_s*entropy_dens

    if(init_type .eq. SHEAR_VISCOUS_1D) eta_vis=eta_over_s !for 1D viscous shear flow init_type

    if(bulkvis) then
      zeta_vis=2.*eta_vis*(1./3.-cs2)
    else
      zeta_vis=0.
    end if
    tau_pi=tau_pi_coeff*eta_vis/(temp*entropy_dens)

    if(init_type .eq. GUBSER_VISCOUS) tau_pi=5.*eta_vis/(temp*entropy_dens)

    tau_pi_big=tau_pi

  end subroutine get_viscous_parameters_0

! ***********************************************************************************

  subroutine get_viscous_parameters(eta_vis,zeta_vis,tau_pi, tau_pi_big, dens, press,errcode)
    use eos
    implicit none
    real(8) :: entropy_dens, energy_dens, temp, cs2
    real(8), intent(in) :: dens, press
    real(8), intent(out) :: eta_vis, zeta_vis, tau_pi, tau_pi_big
    integer, intent(out) :: errcode

    call get_derived_data(dens, press, energy_dens, temp, entropy_dens,errcode)
    call eos_sound(dens,energy_dens,press,cs2,errcode)
    
     eta_vis=hbar*eta_over_s*entropy_dens

    if(init_type .eq. SHEAR_VISCOUS_1D) eta_vis=eta_over_s ! for 1D viscous shear flow init_type
    if(bulkvis) then
      zeta_vis=2.*eta_vis*(1./3.-cs2)
    else
      zeta_vis=0.
    end if
    if((temp  .eq. 0) .or. (entropy_dens .eq. 0) ) then
      tau_pi=0.
    else
    tau_pi=tau_pi_coeff*eta_vis/(temp*entropy_dens) 
    end if
    if(init_type .eq. GUBSER_VISCOUS) tau_pi=5.*eta_vis/(temp*entropy_dens)
    tau_pi_big=tau_pi

  end subroutine get_viscous_parameters
  
! ***********************************************************************************
  
  subroutine get_derived_pi_gcv3(p,pitt,pitx,pity,pitz,gcv3)
    implicit none
    real(8), intent(in), dimension(:) :: p
    real(8), intent(inout) :: pitt,pitx,pity,pitz
    real(8), intent(in) :: gcv3

    !pitt=(g_cov(1)*p(kpixx)+g_cov(2)*p(kpiyy)+gcv3*p(kpizz))/(-g_cov0) !traceless condition
    pitx=p(kpixx)*p(kvx)*g_cov(1)+p(kpixy)*p(kvy)*g_cov(2)+p(kpixz)*p(kvz)*gcv3
    pity=p(kpixy)*p(kvx)*g_cov(1)+p(kpiyy)*p(kvy)*g_cov(2)+p(kpiyz)*p(kvz)*gcv3
    pitz=p(kpixz)*p(kvx)*g_cov(1)+p(kpiyz)*p(kvy)*g_cov(2)+p(kpizz)*p(kvz)*gcv3
    pitt=pitx*p(kvx)*g_cov(1)+pity*p(kvy)*g_cov(2)+pitz*p(kvz)*gcv3
    
  end subroutine get_derived_pi_gcv3

! ***********************************************************************************
  
  subroutine get_derived_pi_zz_gcv3(p,pitt,pitx,pity,pitz,pizz,gcv3)
    implicit none
    real(8), intent(in), dimension(:) :: p
    real(8), intent(inout) :: pitt,pitx,pity,pitz,pizz
    real(8), intent(in) :: gcv3

    pitt=1./(1.-p(kvz)*p(kvz)*gcv3)*(p(kvx)*g_cov(1)*p(kpixy)*p(kvy)*g_cov(2)+p(kvx)*p(kvx)*g_cov(1)*g_cov(1)*p(kpixx)+&
        &p(kvy)*g_cov(1)*p(kpiyz)*p(kvz)*gcv3+2.*p(kvx)*g_cov(1)*p(kpixz)*p(kvz)*gcv3+p(kvy)*g_cov(1)*g_cov(1)*p(kpixy)*&
        &p(kvx)+p(kvy)*p(kvy)*g_cov(1)*p(kpiyy)*g_cov(2)+p(kvz)*gcv3*p(kpiyz)*p(kvy)*g_cov(2)-p(kpixx)*g_cov(1)*p(kvz)*p(kvz)*&
        &gcv3-p(kpiyy)*g_cov(2)*p(kvz)*p(kvz)*gcv3)
    pitx=p(kpixx)*p(kvx)*g_cov(1)+p(kpixy)*p(kvy)*g_cov(2)+p(kpixz)*p(kvz)*gcv3
    pity=p(kpixy)*p(kvx)*g_cov(1)+p(kpiyy)*p(kvy)*g_cov(2)+p(kpiyz)*p(kvz)*gcv3
    pitz=1./(1.-p(kvz)*p(kvz)*gcv3)*(p(kpixz)*p(kvx)*g_cov(1)+p(kpixz)*p(kvx)*g_cov(1)*p(kvz)*p(kvz)*gcv3+p(kpiyz)*p(kvy)*&
        &g_cov(2)+p(kvz)*p(kvx)*g_cov(1)*p(kpixy)*p(kvy)*g_cov(2)-p(kvz)*p(kpixx)*g_cov(1)-p(kvz)*p(kpiyy)*g_cov(2)+p(kvz)*p(kvx)*&
        &p(kvx)*g_cov(1)*g_cov(1)*p(kpixx)+p(kvy)*g_cov(1)*p(kpiyz)*p(kvz)*p(kvz)*gcv3+p(kvz)*p(kvy)*g_cov(1)*g_cov(1)*&
        &p(kpixy)*p(kvx)+p(kvz)*p(kvy)*p(kvy)*g_cov(1)*p(kpiyy)*g_cov(2))
    pizz=1./(gcv3*(1.-p(kvz)*p(kvz)*gcv3))*(p(kvx)*g_cov(1)*p(kpixy)*p(kvy)*g_cov(2)-p(kpixx)*g_cov(1)-p(kpiyy)*g_cov(2)+&
        &p(kvx)*p(kvx)*g_cov(1)*g_cov(1)*p(kpixx)+p(kvy)*g_cov(1)*p(kpiyz)*p(kvz)*gcv3+2.*p(kvx)*g_cov(1)*p(kpixz)*p(kvz)*&
        &gcv3+p(kvy)*g_cov(1)*g_cov(1)*p(kpixy)*p(kvx)+p(kvy)*p(kvy)*g_cov(1)*p(kpiyy)*g_cov(2)+p(kvz)*gcv3*p(kpiyz)*&
        &p(kvy)*g_cov(2))
    
  end subroutine get_derived_pi_zz_gcv3

! ***********************************************************************************
  
  subroutine get_derived_pi(p,pitt,pitx,pity,pitz)
    implicit none
    real(8), intent(in), dimension(:) :: p
    real(8), intent(inout) :: pitt,pitx,pity,pitz

    !pitt=(g_cov(1)*p(kpixx)+g_cov(2)*p(kpiyy)+g_cov(3)*p(kpizz))/(-g_cov0) !traceless condition
    pitx=p(kpixx)*p(kvx)*g_cov(1)+p(kpixy)*p(kvy)*g_cov(2)+p(kpixz)*p(kvz)*g_cov(3)
    pity=p(kpixy)*p(kvx)*g_cov(1)+p(kpiyy)*p(kvy)*g_cov(2)+p(kpiyz)*p(kvz)*g_cov(3)
    pitz=p(kpixz)*p(kvx)*g_cov(1)+p(kpiyz)*p(kvy)*g_cov(2)+p(kpizz)*p(kvz)*g_cov(3)
    pitt=pitx*p(kvx)*g_cov(1)+pity*p(kvy)*g_cov(2)+pitz*p(kvz)*g_cov(3)
    
  end subroutine get_derived_pi

! ***********************************************************************************
  
  subroutine get_derived_pi_zz(p,pitt,pitx,pity,pitz,pizz)
    implicit none
    real(8), intent(in), dimension(:) :: p
    real(8), intent(inout) :: pitt,pitx,pity,pitz,pizz

    pitt=1./(1.-p(kvz)*p(kvz)*g_cov(3))*(p(kvx)*g_cov(1)*p(kpixy)*p(kvy)*g_cov(2)+p(kvx)*p(kvx)*g_cov(1)*g_cov(1)*p(kpixx)+&
        &p(kvy)*g_cov(1)*p(kpiyz)*p(kvz)*g_cov(3)+2.*p(kvx)*g_cov(1)*p(kpixz)*p(kvz)*g_cov(3)+p(kvy)*g_cov(1)*g_cov(1)*p(kpixy)*&
        &p(kvx)+p(kvy)*p(kvy)*g_cov(1)*p(kpiyy)*g_cov(2)+p(kvz)*g_cov(3)*p(kpiyz)*p(kvy)*g_cov(2)-p(kpixx)*g_cov(1)*p(kvz)*p(kvz)*&
        &g_cov(3)-p(kpiyy)*g_cov(2)*p(kvz)*p(kvz)*g_cov(3))
    pitx=p(kpixx)*p(kvx)*g_cov(1)+p(kpixy)*p(kvy)*g_cov(2)+p(kpixz)*p(kvz)*g_cov(3)
    pity=p(kpixy)*p(kvx)*g_cov(1)+p(kpiyy)*p(kvy)*g_cov(2)+p(kpiyz)*p(kvz)*g_cov(3)
    pitz=1./(1.-p(kvz)*p(kvz)*g_cov(3))*(p(kpixz)*p(kvx)*g_cov(1)+p(kpixz)*p(kvx)*g_cov(1)*p(kvz)*p(kvz)*g_cov(3)+p(kpiyz)*p(kvy)*&
        &g_cov(2)+p(kvz)*p(kvx)*g_cov(1)*p(kpixy)*p(kvy)*g_cov(2)-p(kvz)*p(kpixx)*g_cov(1)-p(kvz)*p(kpiyy)*g_cov(2)+p(kvz)*p(kvx)*&
        &p(kvx)*g_cov(1)*g_cov(1)*p(kpixx)+p(kvy)*g_cov(1)*p(kpiyz)*p(kvz)*p(kvz)*g_cov(3)+p(kvz)*p(kvy)*g_cov(1)*g_cov(1)*&
        &p(kpixy)*p(kvx)+p(kvz)*p(kvy)*p(kvy)*g_cov(1)*p(kpiyy)*g_cov(2))
    pizz=1./(g_cov(3)*(1.-p(kvz)*p(kvz)*g_cov(3)))*(p(kvx)*g_cov(1)*p(kpixy)*p(kvy)*g_cov(2)-p(kpixx)*g_cov(1)-p(kpiyy)*g_cov(2)+&
        &p(kvx)*p(kvx)*g_cov(1)*g_cov(1)*p(kpixx)+p(kvy)*g_cov(1)*p(kpiyz)*p(kvz)*g_cov(3)+2.*p(kvx)*g_cov(1)*p(kpixz)*p(kvz)*&
        &g_cov(3)+p(kvy)*g_cov(1)*g_cov(1)*p(kpixy)*p(kvx)+p(kvy)*p(kvy)*g_cov(1)*p(kpiyy)*g_cov(2)+p(kvz)*g_cov(3)*p(kpiyz)*&
        &p(kvy)*g_cov(2))
    
  end subroutine get_derived_pi_zz

! ***********************************************************************************

  subroutine get_sigma(ix,iy,iz,u,sigma, dutdt, dutdx, dutdy, dutdz, duxdx, duydy, duzdz, duzdt, duzdx, duzdy, duxdt,&
                      &duxdy, duxdz, duydt, duydx, duydz)
  implicit none
  integer, intent(in) :: ix,iy,iz
  real(8), dimension(:) :: sigma
  real(8), dimension(:,:,:,:), allocatable :: u
  real(8) :: ut, ux, uy, uz  
  real(8), intent(out) :: dutdt, dutdx, dutdy, dutdz, duxdx, duydy, duzdz, duzdt, duzdx, duzdy, duxdt, duxdy, duxdz, duydt, duydx,&
                       & duydz
  
  real(8) :: ortho_check, trace_check
  
  !these u are the contravariant primitive v
  ux=u(ix,iy,iz,kvx)
  uy=u(ix,iy,iz,kvy)
  uz=u(ix,iy,iz,kvz)
  ut=1./sqrt(1.-g_cov(1)*ux*ux-g_cov(2)*uy*uy-g_cov(3)*uz*uz)
  !not these us are the contravariant u
  ux=ux*ut
  uy=uy*ut
  uz=uz*ut
  dutdt=deriv(ix,iy,iz,dtt)
  dutdx=deriv(ix,iy,iz,dtx)
  dutdy=deriv(ix,iy,iz,dty)
  dutdz=deriv(ix,iy,iz,dtz)
  duxdx=deriv(ix,iy,iz,dxx)
  duydy=deriv(ix,iy,iz,dyy)
  duzdz=deriv(ix,iy,iz,dzz)
  duzdt=deriv(ix,iy,iz,dzt)
  duzdx=deriv(ix,iy,iz,dzx)
  duzdy=deriv(ix,iy,iz,dzy)
  duxdt=deriv(ix,iy,iz,dxt)
  duxdy=deriv(ix,iy,iz,dxy)
  duxdz=deriv(ix,iy,iz,dxz)
  duydt=deriv(ix,iy,iz,dyt)
  duydx=deriv(ix,iy,iz,dyx)
  duydz=deriv(ix,iy,iz,dyz)


  !AAA please, compute only sigmas which are really used (for example, only the spatial part)
  if(coordinates .eq. MINKOWSKI) then !Minkoski coordinates
   
   !THIS SIGMAS ARE NOT USED, BUT THEY ARE LEFT HERE FOR A POSSIBLE FUTURE USE 
   !sigma(stt)=-(2./3.)*dutdt +(2./3.)*dutdt*(ut**2.)+ut*ux*dutdx+ut*uy*dutdy+ut*uz*dutdz+duxdx/3.+duydy/3.+duzdz/3.-(ut**2.)*duxdx&
   !           &/3.-(ut**2.)*duydy/3.-(ut**2.)*duzdz/3.
  
   !sigma(stx)=-0.5*duxdt+0.5*duxdt*(ut**2.)+ut*ux*duxdx/6.+ut*uy*duxdy/2.+ut*uz*duxdz/2.+ut*ux*dutdt/6.+0.5*dutdx+0.5*dutdx*ux**2.&
   !          &+0.5*ux*uy*dutdy+0.5*ux*uz*dutdz-ut*ux*duydy/3.-ut*ux*duzdz/3.
        
   !sigma(sty)=-0.5*duydt+0.5*duydt*(ut**2.)+ut*ux*duydx/2.+ut*uy*duydy/6.+ut*uz*duydz/2.+ut*uy*dutdt/6.+0.5*ux*uy*dutdx+0.5*&
   !          &dutdy+0.5*dutdy*(uy**2.)+0.5*uy*uz*dutdz-ut*uy*duxdx/3.-ut*uy*duzdz/3.

   !sigma(stz)=-0.5*duzdt+0.5*duzdt*(ut**2.)+ut*ux*duzdx/2.+ut*uy*duzdy/2.+ut*uz*duzdz/6.+ut*uz*dutdt/6.+0.5*ux*uz*dutdx+0.5*&
   !          &dutdz+0.5*dutdz*(uz**2.)+0.5*uy*uz*dutdy-ut*uz*duxdx/3.-ut*uz*duydy/3.
 
   sigma(szz)=ut*uz*duzdt+uz*ux*duzdx+uz*uy*duzdy+(2./3.)*duzdz+(2./3.)*duzdz*(uz**2.)-dutdt/3.-duxdx/3.-duydy/3.-(uz**2.)*dutdt/3.&
             &-(uz**2.)*duxdx/3.-(uz**2.)*duydy/3.

   sigma(sxz)=0.5*ut*ux*duzdt+0.5*duzdx+0.5*(ux**2.)*duzdx+0.5*ux*uy*duzdy+ux*uz*duzdz/6.+0.5*ut*uz*duxdt+ux*uz*duxdx/6.+0.5*uy*uz*&
             &duxdy+0.5*duxdz+0.5*duxdz*uz**2.-ux*uz*dutdt/3.-ux*uz*duydy/3.

   sigma(sxy)=0.5*ut*ux*duydt+0.5*duydx+0.5*(ux**2.)*duydx+ux*uy*duydy/6.+ux*uz*duydz/2.+0.5*ut*uy*duxdt+ux*uy*duxdx/6.+0.5*duxdy+&
             &0.5*duxdy*(uy**2.)+uy*uz*duxdz/2.-ux*uy*dutdt/3.-ux*uy*duzdz/3.

   sigma(syz)=0.5*ut*uy*duzdt+ux*uy*duzdx/2.+duzdy/2.+duzdy*(uy**2.)/2.+uy*uz*duzdz/6.+ut*uz*duydt/2.+ux*uz*duydx/2.+uy*uz*duydy/6.&
             &+duydz/2.+duydz/2.*(uz**2.)-uy*uz*dutdt/3.-uy*uz*duxdx/3.

   sigma(sxx)=ut*ux*duxdt+(2./3.)*duxdx+(2./3.)*duxdx*(ux**2.)+ux*uy*duxdy+ux*uz*duxdz-dutdt/3.-duydy/3.-duzdz/3.-(ux**2.)*&
             &dutdt/3.-(ux**2.)*duydy/3.-(ux**2.)*duzdz/3.
            
   sigma(syy)=ut*uy*duydt+ux*uy*duydx+(2./3.)*duydy+(2./3.)*duydy*(uy**2.)+uy*uz*duydz-dutdt/3.-duxdx/3.-duzdz/3.-(uy**2.)*&
             &dutdt/3.-(uy**2.)*duxdx/3.-(uy**2.)*duzdz/3.
  
  else !Bjorken coordinates

   !THIS SIGMAS ARE NOT USED, BUT THEY ARE LEFT HERE FOR A POSSIBLE FUTURE USE 
   ! sigma(stt)=-(2./3.)*dutdt +(2./3.)*dutdt*(ut**2.)+ut*ux*dutdx+ut*uy*dutdy+ut*uz*dutdz+duxdx/3.+duydy/3.+duzdz/3.+ut/(3.*t)-&
   !           &(ut**2.)*duxdx/3.-(ut**2.)*duydy/3.-(ut**2.)*duzdz/3.-ut**3./(3.*t)+t*ut*(uz)**2.
  
   ! sigma(stx)=-0.5*duxdt+0.5*duxdt*(ut**2.)+ut*ux*duxdx/6.+ut*uy*duxdy/2.+ut*uz*duxdz/2.+ut*ux*dutdt/6.+0.5*dutdx+0.5*dutdx*&
   !          &(ux**2.)+0.5*ux*uy*dutdy+0.5*ux*uz*dutdz-ut*ux*duydy/3.-ut*ux*duzdz/3.-ux*(ut**2.)/(3.*t)+t/2.*(uz**2.)*ux
    
   ! sigma(sty)=-0.5*duydt+0.5*duydt*(ut**2.)+ut*ux*duydx/2.+ut*uy*duydy/6.+ut*uz*duydz/2.+ut*uy*dutdt/6.+0.5*ux*uy*dutdx+0.5*&
   !          &dutdy+0.5*dutdy*(uy**2.)+0.5*uy*uz*dutdz-ut*uy*duxdx/3.-ut*uy*duzdz/3.-uy*(ut**2.)/(3.*t)+t/2.*(uz**2.)*uy

   ! sigma(stz)=-0.5*duzdt+0.5*duzdt*(ut**2.)+ut*ux*duzdx/2.+ut*uy*duzdy/2.+ut*uz*duzdz/6.+ut*uz*dutdt/6.+0.5*ux*uz*dutdx+0.5*&
   !          &dutdz/(t**2.)+0.5*dutdz*(uz**2.)+0.5*uy*uz*dutdy-ut*uz*duxdx/3.-ut*uz*duydy/3.-(ut**2.)*uz/(3.*t)+0.5*(uz**3.)*t&
   !          &+(ut**2.)/t*uz
   
   sigma(szz)=ut*uz*duzdt+uz*ux*duzdx+uz*uy*duzdy+(2.*duzdz)/(3.*t**2.)+(2./3.)*duzdz*(uz**2.)-dutdt/(3.*(t**2.))-duxdx/(3.*&
             &(t**2.))-duydy/(3.*(t**2.))-(uz**2.)*dutdt/3.-(uz**2.)*duxdx/3.-(uz**2.)*duydy/3.&
             &+2./3.*ut/(t**3.)+5./3.*(ut/t)*(uz**2.)
   
   sigma(sxz)=0.5*ut*ux*duzdt+0.5*duzdx+0.5*(ux**2.)*duzdx+0.5*ux*uy*duzdy+ux*uz*duzdz/6.+0.5*ut*uz*duxdt+ux*uz*duxdx/6.+0.5*uy*uz*&
             &duxdy+duxdz/(2.*t**2.)+0.5*duxdz*uz**2.-ux*uz*dutdt/3.-ux*uz*duydy/3.+(2./3.)*ut*uz*ux/t
  
   sigma(sxy)=0.5*ut*ux*duydt+0.5*duydx+0.5*(ux**2.)*duydx+ux*uy*duydy/6.+ux*uz*duydz/2.+0.5*ut*uy*duxdt+ux*uy*duxdx/6.+0.5*duxdy+&
             &0.5*duxdy*(uy**2.)+uy*uz*duxdz/2.-ux*uy*dutdt/3.-ux*uy*duzdz/3.-ut*ux*uy/(3.*t)

   sigma(syz)=0.5*ut*uy*duzdt+ux*uy*duzdx/2.+duzdy/2.+duzdy*(uy**2.)/2.+uy*uz*duzdz/6.+ut*uz*duydt/2.+ux*uz*duydx/2.+uy*uz*duydy/6.&
             &+duydz/(2.*(t**2.))+duydz/2.*(uz**2.)-uy*uz*dutdt/3.-uy*uz*duxdx/3.+(2./3.)*ut*uz*uy/t

   sigma(sxx)=ut*ux*duxdt+(2./3.)*duxdx+(2./3.)*duxdx*(ux**2.)+ux*uy*duxdy+ux*uz*duxdz-dutdt/3.-duydy/3.-duzdz/3.-ut/(3.*t)-&
             &(ux**2.)*dutdt/3.-(ux**2.)*duydy/3.-(ux**2.)*duzdz/3.-(ux**2.)*ut/(3.*t)

   sigma(syy)=ut*uy*duydt+ux*uy*duydx+(2./3.)*duydy+(2./3.)*duydy*(uy**2.)+uy*uz*duydz-dutdt/3.-duxdx/3.-duzdz/3.-ut/(3.*t)-&
             &(uy**2.)*dutdt/3.-(uy**2.)*duxdx/3.-(uy**2.)*duzdz/3.-(uy**2.)*ut/(3.*t)
  end if
  

  end subroutine get_sigma
! ***********************************************************************************
  
  subroutine viscous_renormalize(comp)
  implicit none
  real(8), dimension(:) :: comp
  real(8) :: normalization_term, pitt, pitx, pity, pitz
  
  return

 
  if(comp(kpr) .ne. 0) then
    call get_derived_pi(comp,pitt,pitx,pity,pitz)
    normalization_term=(1.+4*(pitt**2.+comp(kpixx)**2.+comp(kpiyy)**2.+comp(kpizz)**2.+2.*pitx**2.+2.*pity**2.+2.*&
                       &pitz**2.+2.*comp(kpixy)**2.+2.*comp(kpixz)**2.+2.*comp(kpiyz)**2.)**2./(9.*comp(kpr)**4.))**0.25
  else
    normalization_term=1.
  end if
  comp(kpixy:kpizz)=comp(kpixy:kpizz)/normalization_term

  end subroutine viscous_renormalize
! ***********************************************************************************
end module viscous
