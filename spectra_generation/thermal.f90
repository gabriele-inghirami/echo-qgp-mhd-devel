! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2019,2021,2023 The ECHO-QGP team                      * 
! *                                                                           *
! *  File: analysis/thermal.f90                                               *
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
! *  Contributors: Giuseppe Pagliara (pagliara@fe.infn.it)                    *
! *                                                                           *
! *  Other modifications: Gabriele Inghirami (g.inghirami@gsi.de)             *
! *                                                                           *
! *  Acknowledgments: Alessandro Drago                                        *
! *                                                                           *
! *****************************************************************************

module work 
  implicit none
  real, dimension(0:3) :: xfo, ufo, dV
  real, dimension(1:3) :: vfo
  real rho, prex
  
  real pre_ex3, ex0, ex1, ex2, ex3, esp, umupmu, den, num(0:3), cnum(0:3)
  real ch, sh, expnophi
  real numtot
  integer ipt, ipart, iphi, irap

  contains
    
! !------------------------------------------------------------------  
! !------------------------------------------------------------------ 
! !------------------------------------------------------------------  
 
   subroutine work_thermal_viscous_3D(hyper_i, hyper_v, bulk, fracden) !! 3D
    use common, only : Tfo,Efo, npart, m, mu, g, bf, pdg_number, name, Energy_int
    use common_thermal, only : nrap, nphi, npt, &
    & phi, cosphi, senphi, &
    ! viscous here !
    & cosphi2, sencosphi, senphi2, &
    !-- 
    & rapidity,pt, mt, &
    & maxpar, &
    & width, charge, isospin, baryon, strange, charme, bottom, &
    & n_decays, &
    & spectra, spectra_tau, spectra_x, spectra_y, spectra_eta, &
    ! debugging here
    & s_0_x, s_0_y, s_0_tau, s_0_eta, spectra_0
    
    implicit none 
    real, intent(in), dimension(1:13) :: hyper_i
    real, intent(in), dimension(1:10) :: hyper_v
    real, intent(in) :: bulk, fracden
    real one_ov_epT2
    real pixx, piyy, pizz, pitt, pitx, pity, pitz, pixy, pixz, piyz
    real f, f0, delta_f, pimunupmupnu
    real sh, sh2, ch, ch2, chsh, pt2, mt2, ptmt 
    real betasq, tausq, w, umudV, pi0mudV
    integer i
    pixy=hyper_v(1)
    pixz=hyper_v(2)
    piyz=hyper_v(3)
    pixx=hyper_v(4)
    piyy=hyper_v(5)
    pizz=hyper_v(6)    
    pitt=hyper_v(7)
    pitx=hyper_v(8)
    pity=hyper_v(9)
    pitz=hyper_v(10)
    ! tau  or t 
    xfo(0)=hyper_i(1)
    ! x
    xfo(1)=hyper_i(2)
    ! y
    xfo(2)=hyper_i(3)
    ! eta or z
    xfo(3)=hyper_i(4)
    
    tausq=xfo(0)*xfo(0)
    
!     baryon density
    rho=hyper_i(5)
    
!     velocities
    vfo(1)=hyper_i(6)
    vfo(2)=hyper_i(7)
    vfo(3)=hyper_i(8)    
    
!     pressure
    prex=hyper_i(9)
    
!     Volume element
    dV(0)=hyper_i(10)	!fm^3
    dV(1)=hyper_i(11)	!fm^3
    dV(2)=hyper_i(12)	!fm^3
!     YOU MUST CHEK in ECHO IF THIS LAST COMPONENT IS 
! 	tau dtau dx dy      (= fm ^ 4)
! 	OR
! 	1/tau dtau dx dy       (= fm ^ 2)
! 	currently it is fm^4 so:
!     dV(3)=hyper_i(13)/xfo(0)
    dV(3)=hyper_i(13)/xfo(0)
!       otherwise     
!     dV(3)=hyper_i(13)	!fm^2
    
    betasq=vfo(1)*vfo(1)+vfo(2)*vfo(2)+vfo(3)*vfo(3)*tausq
    ufo(0)=1.0/sqrt(1-betasq)
    do i=1, 3
      ufo(i)=ufo(0)*vfo(i)
    end do
! ! ---    
    one_ov_epT2=0.5/((Tfo*Tfo)*(efo+prex))
    ! fm^3 GeV^-3
    ! while p^mu p^nu pi^{mu nu} = GeV^2 GeV fm^{-3}
    ! hence delta_f is adimensional:		correct
! ! ---
 
    
!   u^3 p^3 g_33 = u^3 mt/tau sinh(y-eta) tau*tau = u^3 mt sinh(y-eta) tau
!   = (u^3 * tau) mt sinh(y-eta) = pre_ex3 * mt sinh(y-eta)
!   = [pre_ex3 * sinh(y-eta)] mt
!   = ex3 * mt
    pre_ex3=ufo(3)*xfo(0)
    do ipart=1,npart
      !$OMP PARALLEL
      !$OMP DO PRIVATE(irap,ch,sh,ch2,sh2,chsh,ex0,ex3,cnum,ipt,ex1,ex2,expnophi,num,mt2,pt2,ptmt,iphi,umupmu,esp,den,numtot)&
      !$OMP& PRIVATE(pimunupmupnu,f0,delta_f,f)
      do irap=1, nrap
      ch=cosh(rapidity(irap)-xfo(3))
      sh=sinh(rapidity(irap)-xfo(3))
      !
      ch2=ch*ch
      sh2=sh*sh
      chsh=ch*sh
      !
      ex0=ufo(0)*ch
      ex3=pre_ex3*sh
      cnum(0)=ch*dV(0)
      cnum(3)=sh*dV(3)
      do ipt=1,npt
	ex1=ufo(1)*pt(ipt)
	ex2=ufo(2)*pt(ipt)
	expnophi=(ex0-ex3)*mt(ipart, ipt)
	cnum(1)=pt(ipt)*dV(1)
	cnum(2)=pt(ipt)*dV(2)
	
	num(0)=cnum(0)*mt(ipart,ipt)
	num(3)=cnum(3)*mt(ipart,ipt)
	!
	mt2=mt(ipart,ipt)*mt(ipart,ipt)
	pt2=pt(ipt)*pt(ipt)
	ptmt=mt(ipart,ipt)*pt(ipt)
	!
	do iphi=1, nphi
	  umupmu= expnophi - ex1*cosphi(iphi) - ex2*senphi(iphi)			  
	  esp=(umupmu-mu(ipart))/Tfo
!!!!          den=exp(esp)			!! Maxwell
          den=exp(esp)-(1.0*bf(ipart))  !! Bose
	  
	  num(1)=cnum(1)*cosphi(iphi)
	  num(2)=cnum(2)*senphi(iphi)	  
	  
	  numtot=num(0)+num(1)+num(2)+num(3)

    ! 	  ------- debugging purpose
	  spectra_0(ipart, irap, ipt, iphi)= spectra_0(ipart, irap, ipt, iphi) + (numtot/den)
	  s_0_tau(ipart, irap, ipt, iphi)= s_0_tau(ipart, irap, ipt, iphi) + (num(0)/den)
	  s_0_x(ipart, irap, ipt, iphi)= s_0_x(ipart, irap, ipt, iphi) + (num(1)/den)
	  s_0_y(ipart, irap, ipt, iphi)= s_0_y(ipart, irap, ipt, iphi) + (num(2)/den)
	  s_0_eta(ipart, irap, ipt, iphi)= s_0_eta(ipart, irap, ipt, iphi) + (num(3)/den)
! !     	  ------- 

	  ! p=particle momentum pi=shear stress tensor
	  ! this is 
	  !p^mu p^nu pi^alpha pi^beta g_{alpha mu}g_{beta nu}
	  pimunupmupnu=   pitt*ch2*mt2  &
		      & + pixx*pt2*cosphi2(iphi)  &
		      & + piyy*pt2*senphi2(iphi)  &
		      & + pizz*tausq*sh2*mt2      &
		      & -2.0*pitx*ch*ptmt*cosphi(iphi) &
		      & -2.0*pity*ch*ptmt*senphi(iphi) &
		      & -2.0*xfo(0)*pitz*ch*sh*mt2 &
		      & +2.0*pixy*pt2*sencosphi(iphi) &
		      & +2.0*xfo(0)*pixz*sh*ptmt*cosphi(iphi) &
		      & +2.0*xfo(0)*piyz*sh*ptmt*senphi(iphi)       
	  
	  f0=1.0/den
! 	  DOI: 10.1103/PhysRevC.78.034915 ---- eq.(19)
! 	       10.1103/PhysRevC.79.039903, 
! 	       10.1103/PhysRevC.78.034915, 
! 	       10.1103/PhysRevC.79.039903
	  delta_f=pimunupmupnu*f0*(1.0+(f0*bf(ipart)))*one_ov_epT2    ! Quantum
!!!       delta_f=pimunupmupnu*f0*one_ov_epT2        ! Maxwell
	  f=f0+delta_f
	  spectra(ipart, irap, ipt, iphi)= spectra(ipart, irap, ipt, iphi) + (numtot*f)
	  !spectra_tau(ipart, irap, ipt, iphi)= spectra_tau(ipart, irap, ipt, iphi) + (num(0)*f)
	  !spectra_x(ipart, irap, ipt, iphi)= spectra_x(ipart, irap, ipt, iphi) + (num(1)*f)
	  !spectra_y(ipart, irap, ipt, iphi)= spectra_y(ipart, irap, ipt, iphi) + (num(2)*f)
	  !spectra_eta(ipart, irap, ipt, iphi)= spectra_eta(ipart, irap, ipt, iphi) + (num(3)*f)
	  end do ! angle
	end do ! transverse momentum
      end do! particle rapidity
      !$OMP END DO
      !$OMP END PARALLEL
    end do  !particle specie

    !!! TO FIX: these equations need to be corrected
    ! enthalpy
    !w=Efo+prex+bulk
    ! u^mu dSigma_mu
    !umudV=(ufo(0)*dV(0) + ufo(1)*dV(1) + ufo(2)*dV(2) + ufo(3)*dV(3)*tausq)
    !pi^{0,mu} dSigma_mu
    !pi0mudV=(pitt*dV(0) + pitx*dV(1) + pity*dV(2) + pitz*dV(3)*tausq)
    ! T^{0,mu} dSigma_mu
    !Energy_int= Energy_int + w*ufo(0)*umudV + pi0mudV - (prex+bulk)*dV(0)
  return
  end subroutine work_thermal_viscous_3D  
  
! !------------------------------------------------------------------  
! !------------------------------------------------------------------  
! !------------------------------------------------------------------  
 subroutine work_thermal_visco_3D_vort_on(hyper_i, hyper_v, bulk, fracden, deriv_u, deriv_T) 
    use common, only : Tfo,Efo, npart, m, mu, g, bf, pdg_number, name, Energy_int
    use common_thermal, only : nrap, nphi, npt, &
    & phi, cosphi, senphi, &
    & cosphi2, sencosphi, senphi2, &
    & rapidity,pt, mt, &
    & maxpar, &
    & width, charge, isospin, baryon, strange, charme, bottom, &
    & n_decays, &
    & spectra, spectra_tau, spectra_x, spectra_y, spectra_eta, &
    & exter, dtbx, dtby, dtbz, dxby, dxbz, dybz, pola, pola_boost, &
!	debugging
    & s_0_x, s_0_y, s_0_tau, s_0_eta, spectra_0, pola_0, pola_boost_0
    implicit none 
    real, intent(in), dimension(1:13) :: hyper_i
    real, intent(in), dimension(1:10) :: hyper_v
    real, intent(in) :: bulk, fracden
    
    real, intent(in), dimension(0:3,0:3):: deriv_u
    real, intent(in), dimension(0:3):: deriv_T
!         
    real, dimension(0:3):: subpola
    real, dimension(0:3):: subpola_cart, mom
    
    real betasq, tausq, tau_inv, w, umudV, pi0mudV
    integer i, ivort
    
    real one_ov_epT2
    real pixx, piyy, pizz, pitt, pitx, pity, pitz, pixy, pixz, piyz
    real sh, sh2, ch, ch2, chsh, pt2, mt2, ptmt 
    real f, f0, delta_f, pimunupmupnu
    
    real c_t_xy, c_t_xz, c_t_yz, c_x_ty, c_x_tz, c_x_yz
    real c_y_tx, c_y_tz, c_y_xz, c_z_tx, c_z_ty, c_z_xy
    real w_t_xy, w_t_xz, w_t_yz, w_x_ty, w_x_tz, w_x_yz
    real w_y_tx, w_y_tz, w_y_xz, w_z_tx, w_z_ty, w_z_xy
    
    real contr_Pip, eem_inv, dSpf, dSpf0

    real tst(0:3)
    tst=0.0
    num=0.0
    cnum=0.0
    !shear stress components
    pixy=hyper_v(1)
    pixz=hyper_v(2)
    piyz=hyper_v(3)
    pixx=hyper_v(4)
    piyy=hyper_v(5)
    pizz=hyper_v(6)    
    pitt=hyper_v(7)
    pitx=hyper_v(8)
    pity=hyper_v(9)
    pitz=hyper_v(10)
    
    
    ! tau  or t 
    xfo(0)=hyper_i(1)
    ! x
    xfo(1)=hyper_i(2)
    ! y
    xfo(2)=hyper_i(3)
    ! eta or z
    xfo(3)=hyper_i(4)
    
    tausq=xfo(0)*xfo(0)
    tau_inv=1.0/xfo(0)
!     baryon density
    rho=hyper_i(5)
    
!     velocities
    vfo(1)=hyper_i(6)
    vfo(2)=hyper_i(7)
    vfo(3)=hyper_i(8)    
    
!     pressure
    prex=hyper_i(9)
    
!     Volume element
    dV(0)=hyper_i(10)	!fm^3
    dV(1)=hyper_i(11)	!fm^3
    dV(2)=hyper_i(12)	!fm^3
!     YOU MUST CHEK in ECHO IF THIS LAST COMPONENT IS 
! 	tau dtau dx dy      (= fm ^ 4)
! 	OR
! 	1/tau dtau dx dy       (= fm ^ 2)

! 	currently it is fm^4 so:
    dV(3)=hyper_i(13)/xfo(0)	!fm^3
!       otherwise     
!     dV(3)=hyper_i(13)	!fm^2
    
    betasq=vfo(1)*vfo(1)+vfo(2)*vfo(2)+vfo(3)*vfo(3)*tausq
    ufo(0)=1.0/sqrt(1-betasq)
    do i=1, 3
      ufo(i)=ufo(0)*vfo(i)
    end do
    exter=0.0
    call der_beta(deriv_u, deriv_T, Tfo, ufo, xfo(0))
!     print *, xfo,exter
    
! ! ---    
    one_ov_epT2=0.5/((Tfo*Tfo)*(efo+prex))
    ! fm^3 GeV^-3
    ! while p^mu p^nu pi^{mu nu} = GeV^2 GeV fm^{-3}
    ! hence delta_f is adimensional:		correct
! ! ---
    
!   writing the coefficients for the subpola. 
!   the first is the component, the last two are referred to "exter"
!   I set the signs here, so that I can sum them up later with all + 
    c_t_xz= tau_inv*exter(dxbz)
    c_t_yz=-tau_inv*exter(dybz)
!     c_t_xy 
    c_x_tz=-tau_inv*exter(dtbz)
!     c_x_ty
    c_y_tz= tau_inv*exter(dtbz)
!     c_y_tx
    c_z_ty=-tau_inv*exter(dtby)
    c_z_tx= tau_inv*exter(dtbx)
 
!   u^3 p^3 g_33 = u^3 mt/tau sinh(y-eta) tau*tau = u^3 mt sinh(y-eta) tau
!   = (u^3 * tau) mt sinh(y-eta) = pre_ex3 * mt sinh(y-eta)
!   = [pre_ex3 * sinh(y-eta)] mt
!   = ex3 * mt
    pre_ex3=ufo(3)*xfo(0)
    do ipart=1,npart
      !$OMP PARALLEL
      !$OMP DO PRIVATE(irap,ch,sh,ch2,sh2,chsh,ex0,ex3,cnum,c_t_xy,c_x_yz,c_x_ty,c_y_xz,c_y_tx,c_z_xy,ipt,mt2,pt2,ptmt,ex1,ex2)&
      !$OMP& PRIVATE(expnophi,num,w_t_xz,w_t_yz,w_t_xy,w_x_yz,w_x_tz,w_x_ty,w_y_xz,w_y_tz,w_y_tx,w_z_xy,w_z_ty,w_z_tx,iphi,umupmu)&
      !$OMP& PRIVATE(esp,den,numtot,pimunupmupnu,f0,delta_f,f,dSpf,dSpf0,subpola,subpola_cart,mom,contr_Pip,eem_inv,ivort)

      do irap=1, nrap
      ch=cosh(rapidity(irap)-xfo(3))
      sh=sinh(rapidity(irap)-xfo(3))
      !
      ch2=ch*ch
      sh2=sh*sh
      chsh=ch*sh
      !      
      ex0=ufo(0)*ch
      ex3=pre_ex3*sh
      cnum(0)=ch*dV(0)
      cnum(3)=sh*dV(3)
      
      ! subpola
      c_t_xy=-sh*exter(dxby)
      c_x_yz=-tau_inv*ch*exter(dybz)      
      c_x_ty= sh*exter(dtby)
      c_y_xz= ch*tau_inv*exter(dxbz)
      c_y_tx=-sh*exter(dtbx)
      c_z_xy=-tau_inv*ch*exter(dxby)      
      ! from now on all subpola signs are already assigned
      
      do ipt=1,npt
	mt2=mt(ipart,ipt)*mt(ipart,ipt)
	pt2=pt(ipt)*pt(ipt)
	ptmt=mt(ipart,ipt)*pt(ipt)
	
	ex1=ufo(1)*pt(ipt)
	ex2=ufo(2)*pt(ipt)
	expnophi=(ex0-ex3)*mt(ipart, ipt)
	cnum(1)=pt(ipt)*dV(1)
	cnum(2)=pt(ipt)*dV(2)
	
	num(0)=cnum(0)*mt(ipart,ipt)
	num(3)=cnum(3)*mt(ipart,ipt)
	!subpola
	w_t_xz=pt(ipt)      *c_t_xz
	w_t_yz=pt(ipt)      *c_t_yz
	w_t_xy=mt(ipart,ipt)*c_t_xy
	w_x_yz=mt(ipart,ipt)*c_x_yz
	w_x_tz=pt(ipt)      *c_x_tz
	w_x_ty=mt(ipart,ipt)*c_x_ty
	w_y_xz=mt(ipart,ipt)*c_y_xz
	w_y_tz=pt(ipt)      *c_y_tz
	w_y_tx=mt(ipart,ipt)*c_y_tx
	w_z_xy=mt(ipart,ipt)*c_z_xy
	w_z_ty=pt(ipt)      *c_z_ty
	w_z_tx=pt(ipt)      *c_z_tx
	!----
	do iphi=1, nphi
	  umupmu= expnophi - ex1*cosphi(iphi) - ex2*senphi(iphi)			  
	  esp=(umupmu-mu(ipart))/Tfo
! ! ! ! 	    den=exp(esp)			!! Maxwell
	  den=exp(esp)-(1.0*bf(ipart))  !! Bose
	  
	  num(1)=cnum(1)*cosphi(iphi)
	  num(2)=cnum(2)*senphi(iphi)	  
	  
	  numtot=num(0)+num(1)+num(2)+num(3)

	  ! p=particle momentum pi=shear stress tensor
	  ! this is 
	  !p^mu p^nu pi^alpha pi^beta g_{alpha mu}g_{beta nu}
	  pimunupmupnu=   pitt*ch2*mt2  &
		      & + pixx*pt2*cosphi2(iphi)  &
		      & + piyy*pt2*senphi2(iphi)  &
		      & + pizz*tausq*sh2*mt2      &
		      & -2.0*pitx*ch*ptmt*cosphi(iphi) &
		      & -2.0*pity*ch*ptmt*senphi(iphi) &
		      & -2.0*xfo(0)*pitz*ch*sh*mt2 &
		      & +2.0*pixy*pt2*sencosphi(iphi) &
		      & +2.0*xfo(0)*pixz*sh*ptmt*cosphi(iphi) &
		      & +2.0*xfo(0)*piyz*sh*ptmt*senphi(iphi)       
	  
	  f0=1.0/den

    ! 	  ------- debugging purpose
	  spectra_0(ipart, irap, ipt, iphi)= spectra_0(ipart, irap, ipt, iphi) + (numtot*f0)
	  !s_0_tau(ipart, irap, ipt, iphi)= s_0_tau(ipart, irap, ipt, iphi) + (num(0)*f0)
	  !s_0_x(ipart, irap, ipt, iphi)= s_0_x(ipart, irap, ipt, iphi) + (num(1)*f0)
	  !s_0_y(ipart, irap, ipt, iphi)= s_0_y(ipart, irap, ipt, iphi) + (num(2)*f0)
	  !s_0_eta(ipart, irap, ipt, iphi)= s_0_eta(ipart, irap, ipt, iphi) + (num(3)*f0)
! !     	  ------- 	  
	  
! 	  DOI: 10.1103/PhysRevC.78.034915 ---- eq.(19)
! 	       10.1103/PhysRevC.79.039903, 
! 	       10.1103/PhysRevC.78.034915, 
! 	       10.1103/PhysRevC.79.039903
          delta_f=pimunupmupnu*f0*(1.0+(f0*bf(ipart)))*one_ov_epT2    ! Quantum
!!!       delta_f=pimunupmupnu*f0*one_ov_epT2        ! Maxwell

	  f=f0+delta_f
	  
	  dSpf=(numtot*f)
	  dSpf0=(numtot*f0)
	  !---                                                                    |||||
! ! ! !   spectra(ipart, irap, ipt, iphi)= spectra(ipart, irap, ipt, iphi) + dSpf vvvvv	  
	  spectra(ipart, irap, ipt, iphi)= spectra(ipart, irap, ipt, iphi) + (numtot*f)
	  !spectra_tau(ipart, irap, ipt, iphi)= spectra_tau(ipart, irap, ipt, iphi) + (num(0)*f)
	  !spectra_x(ipart, irap, ipt, iphi)= spectra_x(ipart, irap, ipt, iphi) + (num(1)*f)
	  !spectra_y(ipart, irap, ipt, iphi)= spectra_y(ipart, irap, ipt, iphi) + (num(2)*f)
	  !spectra_eta(ipart, irap, ipt, iphi)= spectra_eta(ipart, irap, ipt, iphi) + (num(3)*f)
	  !---
	  !!subpola
	  subpola(0)= w_t_xz*senphi(iphi) + w_t_yz*cosphi(iphi) + w_t_xy
	  subpola(1)= w_x_yz              + w_x_tz*senphi(iphi) + w_x_ty
	  subpola(2)= w_y_xz              + w_y_tz*cosphi(iphi) + w_y_tx
	  subpola(3)= w_z_xy              + w_z_ty*cosphi(iphi) + w_z_tx*senphi(iphi)
	  
	  !cartesian coordinates for subpola --> subpola_cart
	  subpola_cart(0)=subpola(0)*cosh(xfo(3)) + subpola(3)*xfo(0)*sinh(xfo(3))
	  subpola_cart(1)=subpola(1)
	  subpola_cart(2)=subpola(2)
	  subpola_cart(3)=subpola(0)*sinh(xfo(3)) + subpola(3)*xfo(0)*cosh(xfo(3))
! 	  cartesian coordinates for the particle momentum
	  mom(0)=(cosh(xfo(3))*ch + sh*sinh(xfo(3)))*mt(ipart,ipt)
	  mom(1)=pt(ipt)*cosphi(iphi)
	  mom(2)=pt(ipt)*senphi(iphi)
	  mom(3)=(sinh(xfo(3))*ch + sh*cosh(xfo(3)))*mt(ipart,ipt)
	  
	  contr_Pip=mom(1)*subpola_cart(1) + &
	            mom(2)*subpola_cart(2) + &
	            mom(3)*subpola_cart(3) 
! 	  contr_Pip=contr_Pip/(mom(0)*(mom(0)+m(ipart)))
	  eem_inv=1.0/(mom(0)*(mom(0)+m(ipart)))
	  
	  ivort=0
! 	  --- f=f0+delta_f
	  pola(ivort, ipart, irap, ipt, iphi)=pola(ivort, ipart, irap, ipt, iphi) +    mom(0)*subpola_cart(ivort)*dSpf
	  pola_boost(ivort, ipart, irap, ipt, iphi)=pola_boost(ivort, ipart, irap, ipt, iphi) + contr_Pip*dSpf   
! 	  ---- f=f0
	  pola_0(ivort, ipart, irap, ipt, iphi)=pola_0(ivort, ipart, irap, ipt, iphi) +    mom(0)*subpola_cart(ivort)*dSpf0
	  pola_boost_0(ivort, ipart, irap, ipt, iphi)=pola_boost_0(ivort, ipart, irap, ipt, iphi) + contr_Pip*dSpf0
! 	  --------
! !     	polarization  ------- 
	  do ivort=1,3    
	    ! f=f0+ delta_f
	    pola(ivort, ipart, irap, ipt, iphi)=pola(ivort, ipart, irap, ipt, iphi) + subpola_cart(ivort)*dSpf    
	    pola_boost(ivort, ipart, irap, ipt, iphi)= &
	      & pola_boost(ivort, ipart, irap, ipt, iphi) + dSpf*contr_Pip*eem_inv*mom(ivort)
	      
	    ! f=f0
	    pola_0(ivort, ipart, irap, ipt, iphi)=pola_0(ivort, ipart, irap, ipt, iphi) + subpola_cart(ivort)*dSpf0  
	    pola_boost_0(ivort, ipart, irap, ipt, iphi)= &
	      & pola_boost_0(ivort, ipart, irap, ipt, iphi) + dSpf0*contr_Pip*eem_inv*mom(ivort)	 	    
	  end do
! !     	  ------- 	  
	  end do ! angle
	end do ! transverse momentum
      end do! particle rapidity
      !$OMP END DO
      !$OMP END PARALLEL
    end do  !particle specie

    ! TO FIX: these equations need to be corrected
    ! enthalpy
    !w=Efo+prex+bulk
    ! u^mu dSigma_mu
    !umudV=(ufo(0)*dV(0) + ufo(1)*dV(1) + ufo(2)*dV(2) + ufo(3)*dV(3)*tausq)
    !pi^{0,mu} dSigma_mu
    !pi0mudV=(pitt*dV(0) + pitx*dV(1) + pity*dV(2) + pitz*dV(3)*tausq)
    ! T^{0,mu} dSigma_mu
    !Energy_int= Energy_int + w*ufo(0)*umudV + pi0mudV - (prex+bulk)*dV(0)
  return
  end subroutine work_thermal_visco_3D_vort_on
  
  
! !------------------------------------------------------------------  
! !------------------------------------------------------------------  
  subroutine work_thermal_ideal_3D_vort_on(hyper_i, deriv_u, deriv_T)
    use common, only : Tfo,Efo, npart, m, mu, g, bf, pdg_number, name, Energy_int
    use common_thermal, only : nrap, nphi, npt, &
    & phi, cosphi, senphi, &
    & rapidity,pt, mt, &
    & maxpar, &
    & width, charge, isospin, baryon, strange, charme, bottom, &
    & n_decays, &
    & spectra, spectra_tau, spectra_x, spectra_y, spectra_eta, &
    & exter, dtbx, dtby, dtbz, dxby, dxbz, dybz, pola, pola_boost
    implicit none 
    real, intent(in), dimension(1:13) :: hyper_i
    real, intent(in), dimension(0:3,0:3):: deriv_u
    real, intent(in), dimension(0:3):: deriv_T
!         
    real, dimension(0:3):: subpola
    real, dimension(0:3):: subpola_cart, mom
    
    real betasq, tausq, tau_inv, w
    integer i, ivort
    
    real c_t_xy, c_t_xz, c_t_yz, c_x_ty, c_x_tz, c_x_yz
    real c_y_tx, c_y_tz, c_y_xz, c_z_tx, c_z_ty, c_z_xy
    real w_t_xy, w_t_xz, w_t_yz, w_x_ty, w_x_tz, w_x_yz
    real w_y_tx, w_y_tz, w_y_xz, w_z_tx, w_z_ty, w_z_xy
    
    real contr_Pip, eem_inv, dSpf

    real tst(0:3)
    tst=0.0
    num=0.0
    cnum=0.0
    ! Vorticity
    
    
    ! tau  or t 
    xfo(0)=hyper_i(1)
    ! x
    xfo(1)=hyper_i(2)
    ! y
    xfo(2)=hyper_i(3)
    ! eta or z
    xfo(3)=hyper_i(4)
    
    tausq=xfo(0)*xfo(0)
    
!     baryon density
    rho=hyper_i(5)
    
!     velocities
    vfo(1)=hyper_i(6)
    vfo(2)=hyper_i(7)
    vfo(3)=hyper_i(8)    
    
!     pressure
    prex=hyper_i(9)
    
!     Volume element
    dV(0)=hyper_i(10)	!fm^3
    dV(1)=hyper_i(11)	!fm^3
    dV(2)=hyper_i(12)	!fm^3
!     YOU MUST CHEK in ECHO IF THIS LAST COMPONENT IS 
! 	tau dtau dx dy      (= fm ^ 4)
! 	OR
! 	1/tau dtau dx dy       (= fm ^ 2)


! 	currently it is fm^4 so:
    dV(3)=hyper_i(13)/xfo(0)	!fm^3
!       otherwise     
!     dV(3)=hyper_i(13)	!fm^2
    
    betasq=vfo(1)*vfo(1)+vfo(2)*vfo(2)+vfo(3)*vfo(3)*tausq
    ufo(0)=1.0/sqrt(1-betasq)
    do i=1, 3
      ufo(i)=ufo(0)*vfo(i)
    end do
    exter=0.0
    call der_beta(deriv_u, deriv_T, Tfo, ufo, xfo(0))
!     print *, xfo,exter
    
    tau_inv=1.0/xfo(0)
!   writing the coefficients for the subpola. 
!   the first is the component, the last two are referred to "exter"
!   I set the signs here, so that I can sum them up later with all + 
    c_t_xz= tau_inv*exter(dxbz)
    c_t_yz=-tau_inv*exter(dybz)
!     c_t_xy 

    c_x_tz=-tau_inv*exter(dtbz)
!     c_x_ty

    c_y_tz= tau_inv*exter(dtbz)
!     c_y_tx

    c_z_ty=-tau_inv*exter(dtby)
    c_z_tx= tau_inv*exter(dtbx)
 
!   u^3 p^3 g_33 = u^3 mt/tau sinh(y-eta) tau*tau = u^3 mt sinh(y-eta) tau
!   = (u^3 * tau) mt sinh(y-eta) = pre_ex3 * mt sinh(y-eta)
!   = [pre_ex3 * sinh(y-eta)] mt
!   = ex3 * mt
    pre_ex3=ufo(3)*xfo(0)
    do ipart=1,npart
      !$OMP PARALLEL
      !$OMP DO PRIVATE(irap,ch,sh,ex0,ex3,cnum,c_t_xy,c_x_yz,c_x_ty,c_y_xz,c_y_tx,c_z_xy,ipt,ex1,ex2,expnophi,num)&
      !$OMP& PRIVATE(w_t_xz, w_t_yz, w_t_xy, w_x_yz, w_x_tz, w_x_ty, w_y_xz, w_y_tz, w_y_tx, w_z_xy, w_z_ty, w_z_tx)&
      !$OMP& PRIVATE(iphi,umupmu,esp,den,numtot,dSpf,subpola, subpola_cart,mom,contr_Pip,eem_inv, ivort)
      
      do irap=1, nrap
      ch=cosh(rapidity(irap)-xfo(3))
      sh=sinh(rapidity(irap)-xfo(3))
      ex0=ufo(0)*ch
      ex3=pre_ex3*sh
      cnum(0)=ch*dV(0)
      cnum(3)=sh*dV(3)
      
      ! subpola
      c_t_xy=-sh*exter(dxby)
      c_x_yz=-tau_inv*ch*exter(dybz)      
      c_x_ty= sh*exter(dtby)
      c_y_xz= ch*tau_inv*exter(dxbz)
      c_y_tx=-sh*exter(dtbx)
      c_z_xy=-tau_inv*ch*exter(dxby)      
      ! from now on all subpola signs are already assigned
      
      do ipt=1,npt
	ex1=ufo(1)*pt(ipt)
	ex2=ufo(2)*pt(ipt)
	expnophi=(ex0-ex3)*mt(ipart, ipt)
	cnum(1)=pt(ipt)*dV(1)
	cnum(2)=pt(ipt)*dV(2)
	
	num(0)=cnum(0)*mt(ipart,ipt)
	num(3)=cnum(3)*mt(ipart,ipt)
	!subpola
	w_t_xz=pt(ipt)      *c_t_xz
	w_t_yz=pt(ipt)      *c_t_yz
	w_t_xy=mt(ipart,ipt)*c_t_xy
	w_x_yz=mt(ipart,ipt)*c_x_yz
	w_x_tz=pt(ipt)      *c_x_tz
	w_x_ty=mt(ipart,ipt)*c_x_ty
	w_y_xz=mt(ipart,ipt)*c_y_xz
	w_y_tz=pt(ipt)      *c_y_tz
	w_y_tx=mt(ipart,ipt)*c_y_tx
	w_z_xy=mt(ipart,ipt)*c_z_xy
	w_z_ty=pt(ipt)      *c_z_ty
	w_z_tx=pt(ipt)      *c_z_tx
	!----
	do iphi=1, nphi
	  umupmu= expnophi - ex1*cosphi(iphi) - ex2*senphi(iphi)			  
	  esp=(umupmu-mu(ipart))/Tfo
! ! ! ! 	    den=exp(esp)			!! Maxwell
	  den=exp(esp)-(1.0*bf(ipart))  !! Bose
	  
	  num(1)=cnum(1)*cosphi(iphi)
	  num(2)=cnum(2)*senphi(iphi)	  
	  
	  numtot=num(0)+num(1)+num(2)+num(3)
	  dSpf=(numtot/den)
	  spectra(ipart, irap, ipt, iphi)= spectra(ipart, irap, ipt, iphi) + dSpf
    ! 	  ------- debugging purpose
	  spectra_tau(ipart, irap, ipt, iphi)= spectra_tau(ipart, irap, ipt, iphi) + (num(0)/den)
	  spectra_x(ipart, irap, ipt, iphi)= spectra_x(ipart, irap, ipt, iphi) + (num(1)/den)
	  spectra_y(ipart, irap, ipt, iphi)= spectra_y(ipart, irap, ipt, iphi) + (num(2)/den)
	  spectra_eta(ipart, irap, ipt, iphi)= spectra_eta(ipart, irap, ipt, iphi) + (num(3)/den)
! ! !     	  ------- 
	  !!subpola
	  subpola(0)= w_t_xz*senphi(iphi) + w_t_yz*cosphi(iphi) + w_t_xy
	  subpola(1)= w_x_yz              + w_x_tz*senphi(iphi) + w_x_ty
	  subpola(2)= w_y_xz              + w_y_tz*cosphi(iphi) + w_y_tx
	  subpola(3)= w_z_xy              + w_z_ty*cosphi(iphi) + w_z_tx*senphi(iphi)
	  
	  !cartesian coordinates for subpola --> subpola_cart
	  subpola_cart(0)=subpola(0)*cosh(xfo(3)) + subpola(3)*xfo(0)*sinh(xfo(3))
	  subpola_cart(1)=subpola(1)
	  subpola_cart(2)=subpola(2)
	  subpola_cart(3)=subpola(0)*sinh(xfo(3)) + subpola(3)*xfo(0)*cosh(xfo(3))
! 	  cartesian coordinates for the particle momentum
	  mom(0)=(cosh(xfo(3))*ch + sh*sinh(xfo(3)))*mt(ipart,ipt)
	  mom(1)=pt(ipt)*cosphi(iphi)
	  mom(2)=pt(ipt)*senphi(iphi)
	  mom(3)=(sinh(xfo(3))*ch + sh*cosh(xfo(3)))*mt(ipart,ipt)
	  
	  contr_Pip=mom(1)*subpola_cart(1) + &
	            mom(2)*subpola_cart(2) + &
	            mom(3)*subpola_cart(3) 
! 	  contr_Pip=contr_Pip/(mom(0)*(mom(0)+m(ipart)))
	  eem_inv=1.0/(mom(0)*(mom(0)+m(ipart)))
	  
	  ivort=0
	  pola(ivort, ipart, irap, ipt, iphi)=pola(ivort, ipart, irap, ipt, iphi) +    mom(0)*subpola_cart(ivort)*dSpf
	  pola_boost(ivort, ipart, irap, ipt, iphi)=pola_boost(ivort, ipart, irap, ipt, iphi) + contr_Pip*dSpf   
	  
! !     	polarization  ------- 
	  do ivort=1,3    
	    pola(ivort, ipart, irap, ipt, iphi)=pola(ivort, ipart, irap, ipt, iphi) + subpola_cart(ivort)*dSpf    
	    
	    pola_boost(ivort, ipart, irap, ipt, iphi)=pola_boost(ivort, ipart, irap, ipt, iphi) + &
						    & dSpf*contr_Pip*eem_inv*mom(ivort)
	    
	  end do
! !     	  ------- 	  
	  end do ! angle
	end do ! transverse momentum
      end do! particle rapidity
      !$OMP END DO
      !$OMP END PARALLEL
    end do  !particle specie

    w=Efo+prex
    Energy_int=Energy_int+ w*ufo(0)*(ufo(0)*dV(0) + ufo(1)*dV(1) + ufo(2)*dV(2) + ufo(3)*dV(3)*tausq) - prex * dV(0)



! nc=nc+1
! call set_hysu(xfo,nc)
! call fill_histo_vel(ufo)
  return
  end subroutine work_thermal_ideal_3D_vort_on
! !------------------------------------------------------------------  
! !------------------------------------------------------------------  
  subroutine work_thermal_ideal_3D(hyper_i)
    use common, only : Tfo,Efo, npart, m, mu, g, bf, pdg_number, name, Energy_int
    use common_thermal, only : nrap, nphi, npt, &
    & phi, cosphi, senphi, &
    & rapidity,pt, mt, &
    & maxpar, &
    & width, charge, isospin, baryon, strange, charme, bottom, &
    & n_decays, &
    & spectra, spectra_tau, spectra_x, spectra_y, spectra_eta
    implicit none 
    real, intent(in), dimension(1:17) :: hyper_i
    real betasq, tau_inv, tausq, w, delta_rc
    real, parameter :: rc0=0.0955
    integer i
    
    ! tau  or t 
    xfo(0)=hyper_i(1)
    ! x
    xfo(1)=hyper_i(2)
    ! y
    xfo(2)=hyper_i(3)
    ! eta or z
    xfo(3)=hyper_i(4)
    
    tausq=xfo(0)*xfo(0)
    
!     baryon density
    rho=hyper_i(5)
    
!     velocities
    vfo(1)=hyper_i(6)
    vfo(2)=hyper_i(7)
    vfo(3)=hyper_i(8)    
    
!     pressure
    prex=hyper_i(9)
    
!     Volume element
    dV(0)=hyper_i(10)	!fm^3
    dV(1)=hyper_i(11)	!fm^3
    dV(2)=hyper_i(12)	!fm^3
!     YOU MUST CHEK in ECHO IF THIS LAST COMPONENT IS 
! 	tau dtau dx dy      (= fm ^ 4)
! 	OR
! 	1/tau dtau dx dy       (= fm ^ 2)

! 	currently it is fm^4 so:
    dV(3)=hyper_i(13)/xfo(0)	!fm^3   
!       otherwise     
!     dV(3)=hyper_i(13)	!fm^2
    
    betasq=vfo(1)*vfo(1)+vfo(2)*vfo(2)+vfo(3)*vfo(3)*tausq
    ufo(0)=1.0/sqrt(1-betasq)
    do i=1, 3
      ufo(i)=ufo(0)*vfo(i)
    end do

    delta_rc=hyper_i(17)/(2*rc0)
 
!   u^3 p^3 g_33 = u^3 mt/tau sinh(y-eta) tau*tau = u^3 mt sinh(y-eta) tau
!   = (u^3 * tau) mt sinh(y-eta) = pre_ex3 * mt sinh(y-eta)
!   = [pre_ex3 * sinh(y-eta)] mt
!   = ex3 * mt
    pre_ex3=ufo(3)*xfo(0)
    do ipart=1,npart
      !$OMP PARALLEL
      !$OMP DO PRIVATE(irap,ch,sh,ex0,ex3,cnum,ipt,ex1,ex2,expnophi,num,iphi,umupmu,esp,den,numtot)
      do irap=1, nrap
      ch=cosh(rapidity(irap)-xfo(3))
      sh=sinh(rapidity(irap)-xfo(3))
      ex0=ufo(0)*ch
      ex3=pre_ex3*sh
      cnum(0)=ch*dV(0)
      cnum(3)=sh*dV(3)
      do ipt=1,npt
	ex1=ufo(1)*pt(ipt)
	ex2=ufo(2)*pt(ipt)
	expnophi=(ex0-ex3)*mt(ipart, ipt)
	cnum(1)=pt(ipt)*dV(1)
	cnum(2)=pt(ipt)*dV(2)
	
	num(0)=cnum(0)*mt(ipart,ipt)
	num(3)=cnum(3)*mt(ipart,ipt)
	do iphi=1, nphi
	  umupmu= expnophi - ex1*cosphi(iphi) - ex2*senphi(iphi)			  
	  esp=(umupmu-mu(ipart))/Tfo
!!!          den=exp(esp)			!! Maxwell
	  den=exp(esp)-(1.0*bf(ipart))  !! Bose
	  
	  num(1)=cnum(1)*cosphi(iphi)
	  num(2)=cnum(2)*senphi(iphi)	  
	  
	  numtot=num(0)+num(1)+num(2)+num(3)
	  spectra(ipart, irap, ipt, iphi)= spectra(ipart, irap, ipt, iphi) + (numtot/den)*(1+charge(ipart)*delta_rc)
    ! 	  ------- debugging purpose
    !	  spectra_tau(ipart, irap, ipt, iphi)= spectra_tau(ipart, irap, ipt, iphi) + (num(0)/den)
    !	  spectra_x(ipart, irap, ipt, iphi)= spectra_x(ipart, irap, ipt, iphi) + (num(1)/den)
    !	  spectra_y(ipart, irap, ipt, iphi)= spectra_y(ipart, irap, ipt, iphi) + (num(2)/den)
    !	  spectra_eta(ipart, irap, ipt, iphi)= spectra_eta(ipart, irap, ipt, iphi) + (num(3)/den)
! !     	  ------- 
	  end do ! angle
	end do ! transverse momentum
      end do! particle rapidity
      !$OMP END DO
      !$OMP END PARALLEL
    end do  !particle specie

        !w=Efo+prex
        !Energy_int=Energy_int+ w*ufo(0)*(ufo(0)*dV(0) + ufo(1)*dV(1) + ufo(2)*dV(2) + ufo(3)*dV(3)*tausq) - prex * dV(0)



! nc=nc+1
! call set_hysu(xfo,nc)
! call fill_histo_vel(ufo)
  return
  end subroutine work_thermal_ideal_3D
! ! !------------------------------------------------------------------  
! ! !------------------------------------------------------------------  
  subroutine work_thermal_ideal_2D(hyper_i)
  use bessel
  use common, only : Tfo, Efo, npart, m, mu, g, bf, pdg_number, name, &
    & Energy_int
    use common_thermal, only : nrap, nphi, npt, &
    & phi, cosphi, senphi, &
    & rapidity,pt, mt, &
    & maxpar, &
    & width, charge, isospin, baryon, strange, charme, bottom, &
    & n_decays, &
    & spectra, spectra_tau, spectra_x, spectra_y, spectra_eta
    implicit none 
    real, intent(in), dimension(1:13) :: hyper_i
    real betasq, w, besarg
    integer i


!     
    ! tau  or t 
    xfo(0)=hyper_i(1)
    ! x
    xfo(1)=hyper_i(2)
    ! y
    xfo(2)=hyper_i(3)

!     
! !     baryon density
    rho=hyper_i(5)
!     
! !     velocities
    vfo(1)=hyper_i(6)
    vfo(2)=hyper_i(7)
    vfo(3)=hyper_i(8)    
!     
! !     pressure
    prex=hyper_i(9)
!     
! !     Volume element
    dV(0)=hyper_i(10)
    dV(1)=hyper_i(11)
    dV(2)=hyper_i(12)
! !     dV(3)=hyper_i(13)
!     
    betasq=vfo(1)*vfo(1)+vfo(2)*vfo(2)
! ! ! !     +vfo(3)*vfo(3)*xfo(0)*xfo(0)
    ufo(0)=1.0/sqrt(1-betasq)
    do i=1, 2
      ufo(i)=ufo(0)*vfo(i)
    end do
! 
    do ipart=1,npart
        !$OMP PARALLEL
        !$OMP DO PRIVATE(besarg, num,cnum,ipt,esp,den,numtot,iphi)
	do ipt=1,npt
	  besarg=ufo(0)*mt(ipart, ipt)/Tfo
	  num(0)=2.0*bessk1(besarg)*mt(ipart, ipt)*dV(0)
	  
	  cnum(1)= 2.0*bessk0(besarg)*pt(ipt) 
	  ! |^|  per non ripetere due volte questa operazione e non dichiarare una nuova variabile
	  cnum(2)= cnum(1)*dV(2)
	  cnum(1)= cnum(1)*dV(1)
! 	  
	  do iphi=1, nphi
	    num(1)=cnum(1)*cos(phi(iphi))
	    num(2)=cnum(2)*sin(phi(iphi))
	    esp=(ufo(1)*cos(phi(iphi)) + ufo(2)*sin(phi(iphi))) * pt(ipt) + mu(ipart)
	    den=exp(-esp/Tfo)
! 
! 	    
	    numtot=num(0)+num(1)+num(2)
	    spectra(ipart, 1, ipt, iphi)= spectra(ipart, 1, ipt, iphi) + (numtot/den)
!       ! 	  ------- debugging purpose
	    spectra_tau(ipart, 1, ipt, iphi)= spectra_tau(ipart, 1, ipt, iphi) + (num(0)/den)
	    spectra_x(ipart, 1, ipt, iphi)= spectra_x(ipart, 1, ipt, iphi) + (num(1)/den)
	    spectra_y(ipart, 1, ipt, iphi)= spectra_y(ipart, 1, ipt, iphi) + (num(2)/den)
!   ! !     	  ------- 
	  end do ! angle
	end do ! transverse momentum
       !$OMP END DO
       !$OMP END PARALLEL
    end do  !particle specie

  w=Efo+prex
  Energy_int=Energy_int+ w*ufo(0)*(ufo(0)*dV(0) + ufo(1)*dV(1) + ufo(2)*dV(2)) - prex * dV(0)
!   
!   return
  end subroutine work_thermal_ideal_2D
! !------------------------------------------------------------------  
  subroutine work_thermal_visco_2D(hyper_i, hyper_v, bulk, denomin)
    use bessel
  use common, only : Tfo, Efo, npart, m, mu, g, bf, pdg_number, name, &
    & Energy_int
    use common_thermal, only : nrap, nphi, npt, &
    & phi, cosphi, senphi, &
    & rapidity,pt, mt, &
    & maxpar, &
    & width, charge, isospin, baryon, strange, charme, bottom, &
    & n_decays, &
    & spectra, spectra_tau, spectra_x, spectra_y, spectra_eta
    implicit none 
    real, intent(in), dimension(1:13) :: hyper_i
    real, intent(in), dimension(1:10) :: hyper_v
    real, intent(in):: bulk, denomin
    
    real betasq, tau_inv, w, besarg
    integer i

    real one_ov_epT2
    real pixx, piyy, pizz, pitt, pitx, pity, pitz, pixy, pixz, piyz  
    real sh, sh2, ch, ch2, chsh, pt2, mt2, ptmt, mt3, ptmt2, pt2mt  
    real bek2, bek3
    real viscorr0, viscorr1, viscorr2, viscorr3, viscorrtot
    real cphi, sphi, cossindV
    real pittzz, pitxty, pixxzz

! contravariant components of the shear viscosity tensor 
    pixy=hyper_v(1)
    pixz=hyper_v(2)
    piyz=hyper_v(3)
    pixx=hyper_v(4)
    piyy=hyper_v(5)
    pizz=hyper_v(6)    
    pitt=hyper_v(7)
    pitx=hyper_v(8)
    pity=hyper_v(9)
    pitz=hyper_v(10)

    
    ! tau  or t 
    xfo(0)=hyper_i(1)
    ! x
    xfo(1)=hyper_i(2)
    ! y
    xfo(2)=hyper_i(3)
    ! eta or z
!     xfo(3)=hyper_i(4)
    
!     baryon density
    rho=hyper_i(5)
    
!     velocities
    vfo(1)=hyper_i(6)
    vfo(2)=hyper_i(7)
    vfo(3)=hyper_i(8)    
    
!     pressure
    prex=hyper_i(9)
! calcolo di coefficiente che moltiplica correzioni viscose
    one_ov_epT2=0.5/((Tfo*Tfo)*(efo+prex))
    
!     Volume element
    dV(0)=hyper_i(10)
    dV(1)=hyper_i(11)
    dV(2)=hyper_i(12)
!     dV(3)=hyper_i(13)
    
    betasq=(vfo(1)*vfo(1))+(vfo(2)*vfo(2))
! ! !     +vfo(3)*vfo(3)*xfo(0)*xfo(0)
    ufo(0)=1.0/sqrt(1.0-betasq)
    do i=1, 2
      ufo(i)=ufo(0)*vfo(i)
    end do

    do ipart=1,npart
      !$OMP PARALLEL
      !$OMP DO PRIVATE(ipt,mt2,pt2,ptmt2,pt2mt,bek2,bek3,besarg,num,cnum,pittzz,viscorr3,iphi,cphi,sphi,cossindV,pitxty,pixxzz)&
      !$OMP& PRIVATE(esp,den,numtot,viscorr2,viscorr1,viscorr0,viscorrtot)
      do irap=1, nrap
	do ipt=1,npt
	  mt2=mt(ipart,ipt)*mt(ipart,ipt)
	  pt2=pt(ipt)*pt(ipt)
          mt3=mt2*mt(ipart,ipt)
          ptmt2=pt(ipt)*mt2
          pt2mt=pt2*mt(ipart,ipt)         

	  besarg=ufo(0)*mt(ipart, ipt)/Tfo
          bek2=2.0*bessk1(besarg)/besarg+bessk0(besarg)
          bek3=4.0*bek2/besarg+bessk1(besarg)
	  num(0)=2.0*bessk1(besarg)*mt(ipart, ipt)*dV(0)
	  cnum(1)= 2.0*bessk0(besarg)*pt(ipt)
	  cnum(2)= cnum(1)

      ! we evaluate the correction that does not depend on phi (i.e. that with Bessel K3)
          pittzz=pitt+xfo(0)*xfo(0)*pizz
          viscorr3=mt3*bek3*pittzz*dV(0)	  

	  do iphi=1, nphi
             cphi=cos(phi(iphi))
             sphi=sin(phi(iphi))
             cossindV=cphi*dV(1)+sphi*dV(2)
             pitxty=cphi*pitx+sphi*pity
             pixxzz=pt2*(pixx*cphi**2+piyy*sphi**2+2.0*pixy*sphi*cphi)-mt2*xfo(0)*xfo(0)*pizz

	     num(1)=cnum(1)*cos(phi(iphi))*dV(1)
	     num(2)=cnum(2)*sin(phi(iphi))*dV(2)
	     esp=(ufo(1)*cos(phi(iphi)) + ufo(2)*sin(phi(iphi))) * pt(ipt) + mu(ipart)
	     den=exp(-esp/Tfo)    
	     numtot=num(0)+num(1)+num(2)

             viscorr2=ptmt2*bek2*(-2.0*pitxty*dV(0)+pittzz*cossindV)
             viscorr1=bessk1(besarg)*(mt(ipart,ipt)*pixxzz*dV(0)-2.0*pt2mt*pitxty*cossindV)
             viscorr0=bessk0(besarg)*pt(ipt)*pixxzz*cossindV

             viscorrtot=viscorr0+viscorr1+viscorr2+viscorr3
	     spectra(ipart, irap, ipt, iphi)= spectra(ipart, irap, ipt, iphi) + (numtot/den)+2.0*one_ov_epT2*(viscorrtot/den)
      ! 	  ------- debugging purpose
	  !   spectra_tau(ipart, irap, ipt, iphi)= spectra_tau(ipart, irap, ipt, iphi) + (num(0)/den)
	  !   spectra_x(ipart, irap, ipt, iphi)= spectra_x(ipart, irap, ipt, iphi) + (num(1)/den)
	  !   spectra_y(ipart, irap, ipt, iphi)= spectra_y(ipart, irap, ipt, iphi) + (num(2)/den)
  ! !     	  ------- 
	  end do ! angle
	end do ! transverse momentum
      end do! particle rapidity
      !$OMP END DO
      !$OMP END PARALLEL
    end do  !particle specie

  w=Efo+prex
  Energy_int=Energy_int+ w*ufo(0)*(ufo(0)*dV(0) + ufo(1)*dV(1) + ufo(2)*dV(2)) - prex * dV(0)
  
  return
  end subroutine work_thermal_visco_2D

! !------------------------------------------------------------------  
! !------------------------------------------------------------------   
! !-------------------- VORTICITY ROUTINES HERE   -------------------  
! !------------------------------------------------------------------   
  
  subroutine  der_beta(d_i_uj, dTdx, T, u, tau)
    use common_thermal, only : exter, dtbx, dtby, dtbz, dxby,dxbz, dybz
    implicit none 
    real, intent(in),dimension(0:3,0:3)::  d_i_uj ! d_i u^j = du^j/dx^i
    real, intent(in),dimension(0:3)::  dTdx
    real, intent(in), dimension(0:3):: u
    real, intent(in) :: T, tau
    
    real,dimension(0:3,0:3):: der_i_beta_j

    real g(0:3)
    
    integer i, j
    real T_inv
    
    !! beta_rho = u_rho / T 
    !! d_nu beta_rho = (1/T)[ -u^gamma d_nu g_{gamma rho} +
    !!                         g_{gamma rho} d_nu u^gamma 
    !!                       - g_{alpha rho} (u^alpha/T) d_nu T  ]
    !! ECHO provides d_nu T = dT/dx^nu
    !!               u^mu 
    !!               d_nu u^rho = du^rho/dx^nu
    
    g(0)=1.0 	!! there's just 1 index because I'm using diagonal terms
    g(1)=-1.0
    g(2)=-1.0
    g(3)=-tau*tau    
    
    T_inv=1.0/T
    
    der_i_beta_j=0.0
    do i=0,3
      do j=0,3
	if (i /= j) then 
! 	  der_i_beta_j(i,j)=( d_i_u_j(i,j) - g(j)*u(j)*dTdx(i)  )/T
	  der_i_beta_j(i,j)= g(j)*( d_i_uj(i,j)- u(j)*dTdx(i)*T_inv )* T_inv
	endif 
      end do
    end do
    ! The formula is NOT complete, now we add the -2tau u^eta component
    ! i = 0  and j = 3
    der_i_beta_j(0,3)= der_i_beta_j(0,3) - (2.0*tau*u(3)*T_inv)
    
    
!     since in the calculation of the polarization only includes
!     antisym. combinations we can directly store:
  exter(dtbx)= der_i_beta_j(0,1)-der_i_beta_j(1,0)
  exter(dtby)= der_i_beta_j(0,2)-der_i_beta_j(2,0)
  exter(dtbz)= der_i_beta_j(0,3)-der_i_beta_j(3,0)
  exter(dxby)= der_i_beta_j(1,2)-der_i_beta_j(2,1)
  exter(dxbz)= der_i_beta_j(1,3)-der_i_beta_j(3,1)
  exter(dybz)= der_i_beta_j(2,3)-der_i_beta_j(3,2)

  return 
  end subroutine der_beta
  
end module work





! !------------------------------------------------------------------
! !	Valentina Nov 2013
! ! 	want to produce particles
! !-------------------------------------------------------------------

 program thermal
    use common
    use common_thermal
    use io
    use init
!     use output_thermal
    use work
    implicit none 
    
    character*2 dummy
    INTEGER lid, filerror, eof, AllocateStatus
    integer i,j,k, jv
    integer, dimension(:), allocatable ::  cellsintau
    integer cit
    
    real, dimension(:), allocatable :: hysu_visco
    real, dimension (1:17)::  hysu_ide
    real fracden, bulk
    
    
    call io_read_setup() 
    
    lid=index(inputdir,' ')-1
    
!     ihysu=1
!     open(unit=14,status='old',file=input1(1:li)//'.dat',form='unformatted', iostat=filerror, access='stream')

    open(unit=22,status='old',file=inputdir(1:lid)//'hypersurface.txt',form='formatted', iostat=filerror)
    call check_file(filerror, inputdir(1:lid)//'hypersurface.txt')
    
001 format (I3,ES14.7)    
    read (22, *, iostat=eof) freeze_type, Tfo, Efo
    print *, "freeze type:", freeze_type
    select case (freeze_type)    
    case(0)
      print *, "Isothermal hypersurface detected, freeze-out temperature:", Tfo
      print *, "Associated energy density", Efo
    case(1)
      print *, "Constant energy density hypersurface detected, energydensity:", Efo
      print *, "Associated Temperature", Tfo
    case default
      print *, "The freeze-out procedure can be applied only on isothermal or iso-energy hypersurfaces"
      print *, "I am sorry, I am going to quit!!"
      Tfo=freeze_value
    end select
    call io_read_particlelist()
    call init_allocate_spectra()
    
002 format (A2,I8)    
003 format (17(ES14.7, 2x))

012 format (A2,I8)    
013 format (25(ES14.7, 2x))
011 format (3(A24,2x,I12,3x))
023 format (20(ES14.7, 2x))  


    allocate(cellsintau(1000*nx(0)), stat=AllocateStatus)
    if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate cellsintau ***" 
    
    ! !	Reads the  grid of the system: the number of cells, 
! ! 	the physical boundaries (fm) and the interstitial space (fm)
	print*,'---------------------------------------------------'
	print*,'grid of the system '
	write(*,'(a,4(" | ",i5))') 'nx',nx
	write(*,'(a,4(" | ",f5.2))') 'dx',dx	
	write(*,'(a,4(" | ",f7.2,"-> ",f7.2))') 'xlim',xlim
	!dtau is not fixed, we compute it later
    
    call io_print_setup(0, 'thermal ')
    print *, ""    
    print *, "°°°°°°°°°°°°°°°° here we go! °°°°°°°°°°°°°°°°"
    k=1
    cit=0

!     FLAGS HIERARCHY
!     1) dimensionality (1D, 2D, 3D)
!     2) ECHO viscosity
!     3) vorticity
!     4) viscosity correction to the distribution function
! I reset the indentation for the sake of clearness
! =================================================================
! ====================           CORE           ===================
! =================================================================
select case (dimension_flag)
  case(1)
    print *, "1+1D still not implemented"
    call exit
  case(2)  !---------------------------------------------------------
    ! ! !     IDEAL 2D-----------------------------------------------
      print *, "2+1D IDEAL ", nx(0)
      if (vorticity_flag == 1) then 
	print *, "Sorry, Vorticity can be calculated only in 3+1D... I am going to quit"
	call exit
      endif 
      if (viscosity_echo == 0) then 
	DO WHILE (eof .eq. 0)
	  read(22, *, iostat=eof, end=100) cellsintau(k)
	  cit=cit+cellsintau(k)
	  write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	  do i=1, cellsintau(k)
	    read(22, 003, iostat=eof) (hysu_ide(j), j=1, 13)
	    call work_thermal_ideal_2D(hysu_ide)
	  end do
	  print *, "Energy up to now", Energy_int
	  k=k+1
	END DO     
	
      else if(viscosity_echo==1 )then    
	allocate (hysu_visco(1:10), stat=AllocateStatus)
	if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate hysu_visco ***" 
!   ! !     VISCOUS 2D WITH F=F0-----------------------------------------------	
	if (viscosity_spectra == 0) then
	  print *, "2+1D VISCOUS AND f = f0 "      
	  DO WHILE (eof .eq. 0)
	    read(22, *, iostat=eof, end=100) cellsintau(k)
	    cit=cit+cellsintau(k)
	    write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	    do i=1, cellsintau(k)
	      read(22, 013, iostat=eof) (hysu_ide(j), j=1, 13), (hysu_visco(jv), jv=1,10), fracden
	      call work_thermal_ideal_2D(hysu_ide)
	    end do
	    print *, "Energy up to now", Energy_int
	    k=k+1
	  END DO 	    
	else if (viscosity_spectra == 1) then
	  print *, "2+1D VISCOUS  AND f = f0 + delta f " 
	  DO WHILE (eof .eq. 0)
	    read(22, *, iostat=eof, end=100) cellsintau(k)
	    cit=cit+cellsintau(k)
	    write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	    do i=1, cellsintau(k)
	      read(22, 013, iostat=eof) (hysu_ide(j), j=1, 13), bulk, (hysu_visco(jv), jv=1,10), fracden
	      call work_thermal_visco_2D(hysu_ide, hysu_visco, bulk, fracden)
	    end do
	    print *, "Energy up to now", Energy_int
	    k=k+1
	  END DO 	    
	else !! viscosity_spectra
	  print *, "something's wrong in settings.txt, visco_spe=", viscosity_spectra
	  call exit 
	endif	  
      else !!  viscosity_echo
	print *, "something's wrong in settings.txt, visco_hyd=", viscosity_echo
	call exit 
      endif     
      
  case(3)  !-------------------------------------------------------------------------   
  
    if( viscosity_echo==0)then 
      ! ! !     IDEAL 3D     ---------------------------------------------------------  
      print *, "3+1D IDEAL "
      if (vorticity_flag == 0) then!vorticity_flag....................................
	DO WHILE (eof .eq. 0)
	  read(22, *, iostat=eof, end=100) cellsintau(k)
	  cit=cit+cellsintau(k)
	  write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	  do i=1, cellsintau(k)
	    read(22, 003, iostat=eof) (hysu_ide(j), j=1, 17)
	    call work_thermal_ideal_3D(hysu_ide)
	  end do
	  print *, "Energy up to now", Energy_int
	  k=k+1
	END DO
      else if (vorticity_flag == 1) then        !....................................
	open(unit=23,status='old',file=inputdir(1:lid)//'hypersurf_deriv.txt',form='formatted', iostat=filerror)
	call check_file(filerror, inputdir(1:lid)//'hypersurf_deriv.txt')
	DO WHILE (eof .eq. 0)
	  read(22, *, iostat=eof, end=100) cellsintau(k)
	  cit=cit+cellsintau(k)
	  write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	  do i=1, cellsintau(k)
	    read(22, 003, iostat=eof) (hysu_ide(j), j=1, 13)
	    read(23, 023) (u_der(j,0), j=1,3), (u_der(j,1), j=1,3), (u_der(j,2), j=1,3), (u_der(j,3), j=1,3),  (u_der(0,j), j=0,3),&
	    & (T_der(j), j=0,3)
	    call work_thermal_ideal_3D_vort_on(hysu_ide, u_der, t_der)
	  end do
	  print *, "Energy up to now", Energy_int
	  k=k+1
	END DO
      else 
	print *, "**** The vorticity_flag can be either 0 or 1."
	print *, "**** Your choice was", vorticity_flag
	print *, "**** please check settings.txt again and relaunch"
	print *, "**** I am quitting."      
	call exit
      endif !vorticity_flag ................................................................
      
    else if ( viscosity_echo==1)then !-----------------------------------------------------
      Allocate (hysu_visco(1:10), stat=AllocateStatus)
      if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate hysu_visco ***" 
      
      if (vorticity_flag == 0) then!vorticity_flag....................................
	if (viscosity_spectra == 0) then ! __________________________________________________
	  print *, "3+1D VISCOUS  AND f = f0  " 
	  DO WHILE (eof .eq. 0)
	    read(22, *, iostat=eof, end=100) cellsintau(k)
	    cit=cit+cellsintau(k)
	    write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	    do i=1, cellsintau(k)
	      read(22, 013, iostat=eof) (hysu_ide(j), j=1, 13), bulk, (hysu_visco(jv), jv=1, 10), fracden
	      call work_thermal_ideal_3D(hysu_ide)
	    end do
	    print *, "Energy up to now", Energy_int
	    k=k+1
	  END DO	  	  
	else if (viscosity_spectra == 1) then! ______________________________________________
	  print *, "3+1D VISCOUS  AND f = f0 + delta f " 
	  DO WHILE (eof .eq. 0)
	    read(22, *, iostat=eof, end=100) cellsintau(k)
	    cit=cit+cellsintau(k)
	    write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	    do i=1, cellsintau(k)
	      read(22, 013, iostat=eof) (hysu_ide(j), j=1, 13), bulk, (hysu_visco(jv), jv=1, 10), fracden	      
	      call work_thermal_viscous_3D(hysu_ide, hysu_visco, bulk, fracden)
	    end do
	    print *, "Energy up to now", Energy_int
	    k=k+1
	  END DO	  
	else                             ! __________________________________________________
	    print *, "**** The viscosity flag visco_spe can be either 0 or 1."
	    print *, "**** Your choice was", viscosity_spectra
	    print *, "**** please check settings.txt again and relaunch"
	    print *, "**** I am quitting."      
	    call exit
	endif                           ! __________________________________________________
      else if (vorticity_flag == 1) then!vorticity_flag....................................
	open(unit=23,status='old',file=inputdir(1:lid)//'hypersurf_deriv.txt',form='formatted', iostat=filerror)	
	call check_file(filerror, inputdir(1:lid)//'hypersurf_deriv.txt')
	if (viscosity_spectra == 0) then ! __________________________________________________
	  print *, "3+1D VISCOUS  AND f = f0 with the calculation of VORTICITY" 
	  DO WHILE (eof .eq. 0)
	    read(22, *, iostat=eof, end=100) cellsintau(k)
	    cit=cit+cellsintau(k)
! ! 	    cambiare qui 
	    write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	    do i=1, cellsintau(k)
	      read(22, 013, iostat=eof) (hysu_ide(j), j=1, 13), bulk, &
				      & (hysu_visco(jv), jv=1, 10), fracden
	      read(23, 023) (u_der(j,0), j=1,3), (u_der(j,1), j=1,3), &
			   & (u_der(j,2), j=1,3), (u_der(j,3), j=1,3), &
			   & (u_der(0,j), j=0,3), (T_der(j), j=0,3)
	      call work_thermal_ideal_3D_vort_on(hysu_ide, u_der, t_der)
	    end do
	    print *, "Energy up to now", Energy_int
! 	    print *, " " 
	    k=k+1
	  END DO	  	  
	else if (viscosity_spectra == 1) then! ______________________________________________
	  print *, "3+1D VISCOUS  AND f = f0 + delta f with the calculation of VORTICITY" 
	  DO WHILE (eof .eq. 0)
	    read(22, *, iostat=eof, end=100) cellsintau(k)
	    cit=cit+cellsintau(k)
	    write (*, 011) "reading time step number", k, "containing freezing cells:", cellsintau(k), "total frozen cells", cit
	    do i=1, cellsintau(k)
	      read(22, 013, iostat=eof) (hysu_ide(j), j=1, 13), bulk, &
				      & (hysu_visco(jv), jv=1, 10), fracden
	      read(23, 023) (u_der(j,0), j=1,3), (u_der(j,1), j=1,3), &
			   & (u_der(j,2), j=1,3), (u_der(j,3), j=1,3), &
			   & (u_der(0,j), j=0,3), (T_der(j), j=0,3)
	      call work_thermal_visco_3D_vort_on(hysu_ide, hysu_visco,bulk, fracden, u_der, t_der)
	    end do
	    print *, "Energy up to now", Energy_int
	    k=k+1
	  END DO	  	  	  
	else                             ! __________________________________________________
	    print *, "**** The viscosity flag visco_spe can be either 0 or 1."
	    print *, "**** Your choice was", viscosity_spectra
	    print *, "**** please check settings.txt again and relaunch"
	    print *, "**** I am quitting."      
	    call exit
	endif         
      else              !(vorticity_flag)              ....................................
      	print *, "**** The vorticity_flag can be either 0 or 1."
	print *, "**** Your choice was", vorticity_flag
	print *, "**** please check settings.txt again and relaunch"
	print *, "**** I am quitting."      
	call exit
      endif		!(vorticity_flag)              ....................................
    else ! ( viscosity_echo)     !--------------------------------------------------------- 
	print *, "**** The viscosity flag visco_hyd can be either 0 or 1."
	print *, "**** Your choice was", viscosity_echo
	print *, "**** please check settings.txt again and relaunch"
	print *, "**** I am quitting."      
	call exit
    endif !( viscosity_echo      )--------------------------------------------------------- 
  end select
    
    099 print *, "**** There is something wrong, I expected " ,cellsintau(k), "cells"
    print *, "**** but the file hypersurface.txt terminated earlier."
    print *, "**** I am quitting, please contact the developer."
    if (cellsintau(k)>0) call exit
100 continue
! =================================================================
  frozen_cells=cit
  
  call io_end_setup(0)
    open(unit=36,status='replace' ,file=(outdir(1:LID_out)//'utils_therm'//'.txt'),  form='formatted', iostat=filerror)
    call check_file(filerror, outdir(1:	LID_out)//'utils_therm'//'.txt')    
    write (36,'(A12, 2x,ES14.7)') "Energy (GeV)", Energy_int
    write (36,'(A12, 2x,I16)') "Frozen cells", cit
    close (36)
    
    !print *, "Total Energy ", Energy_int, "GeV"
    print *, "Frozen cells", cit
    call io_print_spectra()
    
    
 return   
 end program thermal
