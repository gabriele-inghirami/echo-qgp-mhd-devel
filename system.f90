! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Copyright (C) 2015-2019 The ECHO-QGP team                                * 
! *                                                                           *
! *  File: system.f90                                                         *
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
! *  Contributors: Vinod Chandra (vchandra@iitgn.ac.in)                       *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           *
! *****************************************************************************


module system_eqgp
!-- Equations for (RHD in light-cone coordinates)
!-- See Del Zanna et al. 2007, A&A 473, 11

  use eos
  use viscous

  use common, only: nv,nu,nov,nou,krh,kvx,kvy,kvz,kpr,kpibu,kpizz,kpixz,kpixy,kpiyz,kpixx,kpiyy,ibc
  use common, only: mhd,kbx,kby,kbz,kglm,krc,glm_alpha,glm_ch,divclean,algo_solver,print_rho_comov,out_sel
  use common, only: rmhd,kex,key,kez,sigma_el_cond,sigma_chiral_cond,ratio_chir_el,imex_alpha_coeff,dt_int
  use common, only: coordinates,g_cov,gp,gm,g_cov0,t,init_type,pitx,pity,pitz,pitt,recal,nx,ny,nz,ngc,pizz,obtained
  use common, only: MINKOWSKI, BJORKEN,viscosity,init_type

  implicit none

  procedure(), pointer :: system_metric

  real(8), dimension(5)   :: solver_parameters !array containing the parameters to pass to the non lin. eq. solver for the mhd case
  real(8) :: w_algo_solver7 !for the third degree pol. - rtsafe solver
  real(8) :: w_algo_solver8_e_3p !for the third degree pol. - rtsafe solver
  real(8), dimension(2)   :: X_guess, X_arr !array containing the unkwown for the inner 2x2 nleq system in the algo 4 solver
  real(8) :: xinner !used for the internal 2x2 nlew solver in the 4 algo solver
    
contains

! *****************************************************************************

! we use the ACM library
INCLUDE "libacm554.f90"

! *****************************************************************************
subroutine system_setup(coordinates, system_metric, viscosity, nv,nu,nov,nou)

  procedure(), pointer :: system_metric
  integer :: coordinates, nv,nu,nov,nou
  logical :: viscosity

  select case (coordinates)
    case (MINKOWSKI)
     system_metric=>system_metric_minkowski
    case (BJORKEN)
     system_metric=>system_metric_bjorken
  end select

  if (viscosity) then
     nov=5
     nou=5
     nv=12
  else
     nov=5
     nou=5
     nv=5
  end if

  ! if we run mhd simulations, then we have got other 3 variables (three components for B)
  ! we need also to associate the index labels of B components to numbers
  if(mhd) then
    if(rmhd) then
      nov=nov+6
      nou=nou+6
      nv=nv+6
      kbx=nv-5
      kby=nv-4
      kbz=nv-3
      kex=nv-2
      key=nv-1
      kez=nv
    else
      nov=nov+3
      nou=nou+3
      nv=nv+3
      kbx=nv-2
      kby=nv-1
      kbz=nv
    end if
    if(divclean) then
      nov=nov+1
      nou=nou+1
      nv=nv+1
      kglm=nv
    end if 
    if((print_rho_comov) .or. (out_sel%rc .eq. 1))then
      nov=nov+1
      nou=nou+1
      nv=nv+1
      krc=nv
    end if 
    
  end if
  
  nu=nv

end subroutine system_setup

! *****************************************************************************

subroutine system_metric_minkowski()
!-- Define metric terms

  g_cov(1)=1.
  g_cov(2)=1.
  g_cov(3)=1.

  gp=1.
  gm=1.
  
  g_cov0=-1. !this component is outside the array because it has been added later

end subroutine system_metric_minkowski

! *****************************************************************************

subroutine system_metric_bjorken()
!-- Define metric terms

  g_cov(1)=1.
  g_cov(2)=1.
  g_cov(3)=t**2.
  
  gp=t
  gm=1./gp
  
  g_cov0=-1. !this component is outside the array because it has been added later

end subroutine system_metric_bjorken

! *****************************************************************************

subroutine system_prim2cons(v,u,errcode)
!-- Calculate conservative variables
  implicit none
  real(8),intent(inout ),dimension(nv) :: v
  real(8),intent(out),dimension(nv) :: u

  real(8) :: rh,vx,vy,vz,pr,v2,glf,en,w,ww
  real(8) :: pitx, pity,pitz,pitt
  real(8) :: pizz
  real(8) :: B2,E2,EEM,factor
  real(8), dimension(1:3) :: EXB,E_field_cova
  integer, intent(out) :: errcode

  errcode=0
  
  rh=v(krh)
  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)
  pr=v(kpr)
  
  
  call eos_energy(rh,en,pr,errcode)
  if (errcode .gt. 0) then
     write(*,*) "An error occurred inside subroutine prim2cons when computing energy density."
     run_crashed=.true.
     return
  end if
  
    v2=g_cov(1)*vx*vx+g_cov(2)*vy*vy+g_cov(3)*vz*vz
    glf=1./sqrt(1.-v2)
  
  w=en+pr
  ww=w*glf**2.

  u(krh)=rh*glf

  if(mhd) then
    if(rmhd) then
      E_field_cova(1:2)=v(kex:key)
      E_field_cova(3)=v(kez)*g_cov(3)
      E2=E_field_cova(1)**2+E_field_cova(2)**2+E_field_cova(3)*v(kez)
    else
      call system_compute_E_in_ideal_MHD(v(kvx:kvz),v(kbx:kbz),E_field_cova)
      E2=E_field_cova(1)**2+E_field_cova(2)**2+E_field_cova(3)**2/g_cov(3)
    end if
    call system_EXB(E_field_cova,v(kbx:kbz),EXB)
    B2=v(kbx)**2+v(kby)**2+v(kbz)**2*g_cov(3)
    EEM=0.5*(E2+B2)
  end if
  
  if(viscosity) then
   w=w+v(kpibu)
   ww=w*glf**2.
   
   if(obtained .eq. 'no') then
     call get_derived_pi(v,pitt,pitx,pity,pitz) 
   else
     call get_derived_pi_zz(v,pitt,pitx,pity,pitz,pizz)
     v(kpizz)=pizz
   end if
   

   u(kvx)=(ww*vx+pitx)
   u(kvy)=(ww*vy+pity)
   u(kvz)=(ww*vz+pitz)*g_cov(3)
   u(kpr)=ww-pr-v(kpibu)-pitt*g_cov0 

   u(kpibu)=v(kpibu)*u(krh)

   u(kpizz)=v(kpizz)*u(krh)
   u(kpixz)=v(kpixz)*u(krh)
   u(kpixy)=v(kpixy)*u(krh)
   u(kpiyz)=v(kpiyz)*u(krh)
   u(kpixx)=v(kpixx)*u(krh)
   u(kpiyy)=v(kpiyy)*u(krh)

  else 
   if(mhd) then
  
     u(kvx)=ww*vx+EXB(1)
     u(kvy)=ww*vy+EXB(2)
     u(kvz)=ww*vz*g_cov(3)+EXB(3)
     u(kpr)=ww-pr+EEM
     if(rmhd) then
       u(kbx:kez)=v(kbx:kez)
     else
       u(kbx:kbz)=v(kbx:kbz)
     end if
     if(divclean) u(kglm)=v(kglm)/gp !we add this division just to counter-balance the final multiplication by gp
   else
     u(kvx)=ww*vx
     u(kvy)=ww*vy
     u(kvz)=ww*vz*g_cov(3)
     u(kpr)=ww-pr
   end if
  end if
  
  u(1:nu)=u(1:nu)*gp
end subroutine system_prim2cons

! *****************************************************************************

subroutine system_prim2cons_0(v,u,errcode)
!-- Calculate conservative variables at initialization stage (only when viscosity is taken into account)
  implicit none
  real(8),intent(inout ),dimension(nv) :: v
  real(8),intent(out),dimension(nv) :: u

  real(8) :: rh,vx,vy,vz,pr,v2,glf,en,w,ww
  real(8) :: pitt,pitx,pity,pitz
  real(8) :: pizz
  integer, intent(out) :: errcode

  errcode=0
  
  rh=v(krh)
  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)
  pr=v(kpr)
  
  call eos_energy(rh,en,pr,errcode)
  if (errcode .gt. 0) then
     write(*,*) "An error occurred into system_prim2cons_0 when computing energy density"
     run_crashed=.true.
     return
  end if
  
   glf=1./sqrt(1.-vx*vx*g_cov(1)-vy*vy*g_cov(2)-vz*vz*g_cov(3))
   w=en+pr+v(kpibu)
   ww=w*glf**2.
  
   u(krh)=rh*glf
   
   if(obtained .eq. 'no') then
     call get_derived_pi(v,pitt,pitx,pity,pitz)
   else
     call get_derived_pi_zz(v,pitt,pitx,pity,pitz,pizz)
     v(kpizz)=pizz
   end if
   
   u(kvx)=ww*vx+pitx
   u(kvy)=ww*vy+pity
   u(kvz)=(ww*vz+pitz)*g_cov(3)
   u(kpr)=ww-pr-v(kpibu)-pitt*g_cov0 

   u(kpibu)=v(kpibu)*u(krh)

   u(kpizz)=v(kpizz)*u(krh)
   u(kpixz)=v(kpixz)*u(krh)
   u(kpixy)=v(kpixy)*u(krh)
   u(kpiyz)=v(kpiyz)*u(krh)
   u(kpixx)=v(kpixx)*u(krh)
   u(kpiyy)=v(kpiyy)*u(krh)

  u(1:nu)=u(1:nu)*gp

end subroutine system_prim2cons_0

! *****************************************************************************

subroutine system_cons2prim0(u,v,ierr,errcode)
!-- Calculate primitive variables in the ideal HD case

  implicit none
  real(8),parameter :: tol=1.d-12, eps=1.d-12

  real(8),intent(inout),dimension(1:5) :: u,v
  real(8), dimension(nv) :: uin

  integer,intent(out) :: ierr

  integer :: iter,iter_max=100

  real(8) :: d,et,sx,sy,sz,pr,ww
  real(8) :: s2,w,v2,rh,en,f,drhdpr,dendpr,pr_eos,dprdrh,dprden,dpr,df,sv
  real(8) :: vx, vy, vz, glf, glfg
 
  integer :: errcode
  
  errcode=0

  d =u(krh)
  sx=u(kvx)
  sy=u(kvy)
  sz=u(kvz)
  et=u(kpr)

  s2=g_cov(1)*sx*sx+g_cov(2)*sy*sy+g_cov(3)*sz*sz

  ierr=0

!-- Old value of pressure provided as initial guess
  pr=u(kpr)/3.
  w=et+pr
  v2=s2/w**2
  !we reduce the guessed initial pressure until we don't have superluminal velocities anymore
  do while (v2 .ge. 1)
     pr=pr*0.8
     w=et+pr
     v2=s2/w**2
     if(pr .lt. 1.e-10) then
       v(krh)=1.0
       v(kvx:kpr)=0.
       return
     end if
   end do 

!-- Newton-Raphson iteration to solve f(p)=eos(rh,en)-p=0
  do iter=1,iter_max

    w=et+pr
    v2=s2/w**2
    if(v2 .ge. 1.) then
      write(*,*) 'Error in ideal cons2prim, superluminal velocity: v2=',v2
      errcode=777
      write(*,*) 'Error:',errcode
      run_crashed=.true.
      return
    end if
    glf=1./sqrt(1.-v2)
    rh=d/glf
    en=w*(1.-v2)-pr

    drhdpr=rh*v2*glf/(en+pr)
    dendpr=v2
    call eos_pressure(rh,en,pr_eos,dprdrh,dprden,errcode)
    f=pr_eos-pr
    df=dprdrh*drhdpr+dprden*dendpr-1.
    dpr=-f/df
    

    if (abs(dpr)<tol*pr) exit
    pr=pr+dpr

  end do

  v(krh)=rh
  v(kvx)=sx/w
  v(kvy)=sy/w
  v(kvz)=sz/w
  v(kpr)=pr
  
  if (iter>=iter_max) ierr=iter
 
  if(pr<=0) then
    write(*,*) "Problems in cons2prim0:"
    write(*,*) "Pressure is negative: ",pr
    call exit(3)
  end if


end subroutine system_cons2prim0


! *****************************************************************************


subroutine system_cons2prim(u,v,ierr,errcode)
!-- Calculate primitive variables

  implicit none
  real(8),parameter :: tol=1.d-12, eps=1.d-12

  real(8),intent(inout),dimension(nv) :: u,v
  real(8), dimension(nv) :: uin

  integer,intent(out) :: ierr

  integer :: iter,iter_max=100

  real(8) :: d,et,sx,sy,sz,pr,ww
  real(8) :: s2,w,v2,rh,en,f,drhdpr,dendpr,pr_eos,dprdrh,dprden,dpr,df,sv
  real(8) :: vx, vy, vz, glf, glfg
  real(8) :: pitt,pitx,pity,pitz
  real(8) :: pizz
 
  integer :: errcode
  integer :: i,k,iv,ip !just counters
  
  real(8) :: sx_cov, sy_cov, sz_cov,en1
  real(8), parameter :: tolf=1.d-12
  real(8), parameter :: tolx=1.d-12
  real(8), dimension(1:nv) :: vprim,vprim1
  real(8), parameter :: epsp=1.d-6

  integer, parameter :: ivmax=100, ipmax=100

  real(8) :: B2,E2,EEM,SB,SS2,SB2
  real(8), dimension(1:3) :: E_field, EXB
 
  real(8) :: bx,by,bz,ex,ey,ez

  real(8) :: xs, ys, DD, UU
  ! variables needed for the Brentm solver
  real(8), dimension(2) :: X_eq, F_eq
  real(8), dimension(3) :: X3_eq, F3_eq
  real(8), dimension(3,3) :: fjac3
  integer :: info_solver, nfev
  real(8), dimension(2) :: sigma_brentm2,wa1_brentm2,wa2_brentm2
  real(8), dimension(3) :: sigma_brentm3,wa1_brentm3,wa2_brentm3
  

  real(8), dimension(18) :: workarr !working array for cernlib, with dimensions at least N*(N+3) (so N=3->3*(3+3)=18)

  integer, parameter :: LWA=(3*(3+13))
  real(8), dimension(LWA) :: WA

  real(8) :: ps, v2old

  !for NR with backtracking
  real(8), parameter :: smallvalue=4.d-16, accuracy=1.d-14
  real(8), dimension(1:3,1:3) :: matL ! L is the triangular lower matrix
  real(8), dimension(1:3) :: DX
  integer, dimension(1:3) :: PER !PER is the permutation matrix for pivoting
  real(8) :: Fnorm, Fnorm_lambda
  integer :: iterations
  integer, parameter :: maxiterations=300 !maximum number of iterations before leaving
  real(8) :: lambda
  logical :: keepon

  !for third degree pol. - NR approach
  real(8) :: gam,v2_old,pg,vb

  !for solver 9
  real(8) :: b9,c9,d9,e9,p9,q9,sbig9,qbig9,dnull9,deinz9,droot9,discr9,drad9,phi9
  real(8) :: sol1,sol2,sol3,sol4,rad129,rad349,root129,root349,rapp9,sol9,dist9,d19,d29,d39,d49,oldsol9,v9
  logical :: exact9

  !for algo_solver 7 (IMEX-SSP resistive MHD)
  real(8),dimension(1:3,1:3) :: MM
  real(8) :: den_M, a_M, s_M, AG, g3, AG2, AG3, gl3, vx2, vy2, vz2, g3m, g3z, vdotB
  real(8), dimension(1:3) :: vel_arr, s_arr, E_arr, vXB, Stiff_arr, k_arr

  errcode=0

  d =gm*u(krh)
  sx=gm*u(kvx)
  sy=gm*u(kvy)
  sz=gm*u(kvz)/g_cov(3)
  et=gm*u(kpr)


  
  s2=g_cov(1)*sx*sx+g_cov(2)*sy*sy+g_cov(3)*sz*sz

  if ((.not. viscosity) .and. (.not. mhd)) then 

  ierr=0

!-- Old value of pressure provided as initial guess
  pr=v(kpr)

!-- Newton-Raphson iteration to solve f(p)=eos(rh,en)-p=0
  do iter=1,iter_max

    w=et+pr
    v2=s2/w**2.
    if(v2 .ge. 1.) then
      write(*,*) 'Error in ideal cons2prim, superluminal velocity: v2=',v2
      errcode=777
      write(*,*) 'Error:',errcode
      run_crashed=.true.
      return
    end if
    glf=1./sqrt(1.-v2)
    rh=d/glf
    en=w*(1.-v2)-pr

    drhdpr=rh*v2*glf/(en+pr)
    dendpr=v2
    call eos_pressure(rh,en,pr_eos,dprdrh,dprden,errcode)
    f=pr_eos-pr
    df=dprdrh*drhdpr+dprden*dendpr-1.
    dpr=-f/df
    

    if (abs(dpr)<tol*pr) exit
    pr=pr+dpr

  end do

  v(krh)=rh
  v(kvx)=sx/w
  v(kvy)=sy/w
  v(kvz)=sz/w
  v(kpr)=pr
  
  if (iter>=iter_max) ierr=iter
 
  if(pr<=0) then
    write(*,*) "Problems in cons2prim:"
    write(*,*) "Pressure is negative: ",pr
    errcode=777
    write(*,*) 'Error:',errcode
    run_crashed=.true.
  end if


  else if(viscosity) then
  ! CASE WITH VISCOSITY

  ierr=0
!  
  uin(:)=gm*u(:)
  
  d=uin(krh) 
!  
  v(kpibu:kpizz)=uin(kpibu:kpizz)/uin(krh)

!-- v^i as the old values

  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)

  v2=g_cov(1)*vx**2+g_cov(2)*vy**2+g_cov(3)*vz**2

  vprim=v !to be used with get_derived_pi

!-- External cycle on v

  do iv=1, ivmax

    vprim(kvx)=vx
    vprim(kvy)=vy
    vprim(kvz)=vz
    
    if(obtained .eq. 'no') then
      call get_derived_pi(vprim,pitt,pitx,pity,pitz)
    else
      call get_derived_pi_zz(vprim,pitt,pitx,pity,pitz,pizz)
      vprim(kpizz)=pizz
    end if
     

    en1=et+v(kpibu)+pitt*g_cov0
     
    sx_cov=uin(kvx)-pitx
    sy_cov=uin(kvy)-pity
    sz_cov=uin(kvz)-pitz*g_cov(3)

    sx=sx_cov
    sy=sy_cov
    sz=sz_cov/g_cov(3)

    s2=sx*sx_cov+sy*sy_cov+sz*sz_cov

    pr=v(kpr)

    
    do ip=1,ipmax
      
      v2=s2/(en1+pr)**2.
      if(v2 .ge. 1.) then
        write(*,*) 'Error in viscous cons2prim, superluminal velocity: v2=',v2
        errcode=107
        write(*,*) 'Error:',errcode
        run_crashed=.true.
        return
      end if


      glf=1./sqrt(1. - v2) 

      rh=d/glf
     
      en=(en1+pr)*(1.-v2)-pr-vprim(kpibu)
      drhdpr=rh*v2*glf/(en+pr)
      dendpr=v2
      call eos_pressure(rh,en,pr_eos,dprdrh,dprden,errcode)

      f=pr_eos-pr
      df=dprdrh*drhdpr+dprden*dendpr-1.
      dpr=-f/df
      if (abs(dpr)<tol*pr) exit
      pr=pr+dpr

    end do
    
    vx=sx/(en1+pr); if (abs(vx)<eps) vx=0.
    vy=sy/(en1+pr); if (abs(vy)<eps) vy=0.
    vz=sz/(en1+pr); if (abs(vz)<eps) vz=0.
    
    vprim(kvx)=vx
    vprim(kvy)=vy
    vprim(kvz)=vz

    if (abs(g_cov(1)*vx**2+g_cov(2)*vy**2+g_cov(3)*vz**2-v2)<eps) exit


  end do

  if (iv>=ivmax) ierr=5

  v(krh:kpr)=(/ rh,vx,vy,vz,pr /)

  if(obtained .eq. 'zz') then
    call get_derived_pi_zz(v,pitt,pitx,pity,pitz,pizz)
    v(kpizz)=pizz
  end if
 
  else if(mhd) then
  if(algo_solver .eq. 1) then
  ierr=0
!  
  uin(:)=gm*u(:)
  
  DD=uin(krh)
  sx=uin(kvx)
  sy=uin(kvy)
  sz=uin(kvz)
  UU=uin(kpr)
  bx=uin(kbx) 
  by=uin(kby) 
  bz=uin(kbz) 

  v(kbx)=uin(kbx)
  v(kby)=uin(kby)
  v(kbz)=uin(kbz)


  B2=bx*bx+by*by+bz*bz*g_cov(3)
  SB=sx*bx+sy*by+sz*bz
  SB2=SB*SB
  SS2=sx*sx+sy*sy+sz*sz/g_cov(3)
!-- v^i and pr as the old values
 
  rh=v(krh)
  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)
  pr=v(kpr)
  call eos_energy(rh,en,pr,errcode)

  v2=g_cov(1)*vx**2+g_cov(2)*vy**2+g_cov(3)*vz**2

  xs=v2
  ys=(en+pr)/(1-v2)

  X3_eq(1)=xs
  X3_eq(2)=ys
  X3_eq(3)=pr

  solver_parameters(1)=B2
  solver_parameters(2)=SB2
  solver_parameters(3)=SS2
  solver_parameters(4)=UU
  solver_parameters(5)=DD

  !call DSNLEQ(3,X3_eq,F3_eq,1.e-14,1.e-14,300,0,info_solver,sist3,workarr)
  call BRENTM(sist3,3,X3_eq,F3_eq,real(1.d-14,8),real(1.d-14,8),300,3,info_solver,nfev,workarr,3,sigma_brentm3,wa1_brentm3,&
             &wa2_brentm3)
  xs=X3_eq(1)
  ys=X3_eq(2)
  pr=X3_eq(3)

  !insert a check to avoid xs>1
  rh=DD*sqrt(1-xs)
  vx=(sx+SB*bx/ys)/(ys+B2) 
  vy=(sy+SB*by/ys)/(ys+B2) 
  vz=(sz/g_cov(3)+SB*bz/ys)/(ys+B2) 

  v(krh:kpr)=(/ rh,vx,vy,vz,pr /)

  if(mhd .and. divclean) v(kglm)=u(kglm) !the artificial variable remains the same, both as primitive or conservative 
   
  else if(algo_solver .eq. 2) then
 
  ierr=0
!  
  uin(:)=gm*u(:)
  
  DD=uin(krh)
  sx=uin(kvx)
  sy=uin(kvy)
  sz=uin(kvz)
  UU=uin(kpr)
  bx=uin(kbx) 
  by=uin(kby) 
  bz=uin(kbz) 

  v(kbx)=uin(kbx)
  v(kby)=uin(kby)
  v(kbz)=uin(kbz)


  B2=bx*bx+by*by+bz*bz*g_cov(3)
  SB=sx*bx+sy*by+sz*bz
  SB2=SB*SB
  SS2=sx*sx+sy*sy+sz*sz/g_cov(3)
  
  solver_parameters(1)=B2
  solver_parameters(2)=SB2
  solver_parameters(3)=SS2
  solver_parameters(4)=UU
  solver_parameters(5)=DD
  
!-- v^i and pr as the old values
 
  rh=v(krh)
  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)
  pr=v(kpr)
  call eos_energy(rh,en,pr,errcode)
  v2old=g_cov(1)*vx**2+g_cov(2)*vy**2+g_cov(3)*vz**2
  
  !initial guess for the internal 2x2 system
  X_eq(1)=(en+pr)/(1-v2old)
  X_eq(2)=pr

  !call DSNLEQ(2,X_eq,F_eq,1.e-14,1.e-14,300,0,info_solver,sist4,workarr)
  call BRENTM(sist4,2,X_eq,F_eq,real(1.d-14,8),real(1.d-14,8),300,3,info_solver,nfev,workarr,2,sigma_brentm2,wa1_brentm2,&
             &wa2_brentm2)

  ys=X_eq(1)
  pr=X_eq(2)
  v2=xinner

  !insert a check to avoid xs>1
  rh=DD*sqrt(1-v2)
  vx=(sx+SB*bx/ys)/(ys+B2) 
  vy=(sy+SB*by/ys)/(ys+B2) 
  vz=(sz/g_cov(3)+SB*bz/ys)/(ys+B2) 
  
  v(krh:kpr)=(/ rh,vx,vy,vz,pr /)

  if(mhd .and. divclean) v(kglm)=u(kglm) !the artificial variable remains the same, both as primitive or conservative 

else if(algo_solver .eq. 3) then

  ierr=0
!  
  uin(:)=gm*u(:)

  DD=uin(krh)
  sx=uin(kvx)
  sy=uin(kvy)
  sz=uin(kvz)
  UU=uin(kpr)
  bx=uin(kbx)
  by=uin(kby)
  bz=uin(kbz)

  v(kbx)=uin(kbx)
  v(kby)=uin(kby)
  v(kbz)=uin(kbz)


  B2=bx*bx+by*by+bz*bz*g_cov(3)
  SB=sx*bx+sy*by+sz*bz
  SB2=SB*SB
  SS2=sx*sx+sy*sy+sz*sz/g_cov(3)
!-- v^i and pr as the old values

  rh=v(krh)
  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)
  pr=v(kpr)
  call eos_energy(rh,en,pr,errcode)

  v2=g_cov(1)*vx**2+g_cov(2)*vy**2+g_cov(3)*vz**2

  xs=v2
  ys=(en+pr)/(1-v2)

  X3_eq(1)=xs
  X3_eq(2)=ys
  X3_eq(3)=pr

  solver_parameters(1)=B2
  solver_parameters(2)=SB2
  solver_parameters(3)=SS2
  solver_parameters(4)=UU
  solver_parameters(5)=DD

  call getfun(X3_eq,F3_eq,fjac3,Fnorm,3,1) !this subroutine returns F(X), its Jacobian and its norm
      
  iterations=0

  do while(Fnorm >= accuracy)
            
     call ludecompose(fjac3,3,matL,PER,smallvalue)

     call solve(matL,fjac3,PER,DX,-F3_eq,3)
         
     !backtracking process
     lambda=1.
     call getfun(X3_eq+lambda*DX,F3_eq,fjac3,Fnorm_lambda,3,1)
     if((Fnorm_lambda > (1-lambda/2)*Fnorm) .and. (lambda > 1./64)) then
       keepon=.True.
       do while(keepon)
          lambda=lambda/2.
          call getfun(X3_eq+lambda*DX,F3_eq,fjac3,Fnorm_lambda,3,0)
          if(.not.((Fnorm_lambda > (1-lambda/2)*Fnorm) .and. (lambda > 1./64))) then
            keepon=.False.
            X3_eq=X3_eq+lambda*DX
            Fnorm=Fnorm_lambda
          end if
       end do 
     else
       X3_eq=X3_eq+lambda*DX
       Fnorm=Fnorm_lambda
     end if 
 
     iterations=iterations+1
     if(iterations .gt. maxiterations) then
       write(*,*) "Solver in cons2prim failed, maximum number of iterations reached... Quitting..."
       call exit(8)
     end if           
  end do
 
  xs=X3_eq(1)
  ys=X3_eq(2)
  pr=X3_eq(3)

  !insert a check to avoid xs>1
  rh=DD*sqrt(1-xs)
  vx=(sx+SB*bx/ys)/(ys+B2)
  vy=(sy+SB*by/ys)/(ys+B2)
  vz=(sz/g_cov(3)+SB*bz/ys)/(ys+B2)

  v(krh:kpr)=(/ rh,vx,vy,vz,pr /)

  if(mhd .and. divclean) v(kglm)=u(kglm) !the artificial variable remains the same, both as primitive or conservative 

else if(algo_solver .eq. 4) then
! solver for ideal gas eos
  ierr=0
  uin(:)=gm*u(:)

  gam=(gamma_adindex-1)/gamma_adindex

  D =uin(krh)
  sx=uin(kvx)
  sy=uin(kvy)
  sz=uin(kvz)
  UU=uin(kpr)
  bx=uin(kbx)
  by=uin(kby)
  bz=uin(kbz)



  S2=sx*sx+sy*sy+sz*sz/g_cov(3)
  B2=bx*bx+by*by+bz*bz*g_cov(3)
  SB=(sx*bx+sy*by+sz*bz)
  SB2=SB**2

  solver_parameters(1)=B2
  solver_parameters(2)=SB2
  solver_parameters(3)=S2
  solver_parameters(4)=UU
  solver_parameters(5)=D


  v2_old=v(kvx)**2+v(kvy)**2+v(kvz)**2
  v2=rtsafe(solver7,v2_old,0.d0,1.d0-1.d-11,1.d-12,ierr)

  if (ierr>0) return

  v(1:nv)=uin(1:nv)
  
  if(mhd .and. divclean) v(kglm)=u(kglm) !the artificial variable remains the same, both as primitive or conservative 

  w=w_algo_solver7
  rh=D*sqrt(1.-v2)
  pg=gam*(w*(1.-v2)-rh)

  v(krh)=rh
  v(kpr)=max(pg,eps) !it is important to keep it to prevent crashes for isolated issues

  vb=SB/w

  v(kvx)=(sx+vb*bx)/(w+B2)!; if (abs(vx)<eps) vx=0.
  v(kvy)=(sy+vb*by)/(w+B2)!; if (abs(vy)<eps) vy=0.
  v(kvz)=(sz/g_cov(3)+vb*bz)/(w+B2)!; if (abs(vz)<eps) vz=0.

else if(algo_solver .eq. 5) then

  ierr=0
  uin(:)=gm*u(:)

  D =uin(krh)
  sx=uin(kvx)
  sy=uin(kvy)
  sz=uin(kvz)
  UU=uin(kpr)
  bx=uin(kbx)
  by=uin(kby)
  bz=uin(kbz)

  S2=sx*sx+sy*sy+sz*sz/g_cov(3)
  B2=bx*bx+by*by+bz*bz*g_cov(3)
  SB=(sx*bx+sy*by+sz*bz)
  SB2=SB**2

  solver_parameters(1)=B2
  solver_parameters(2)=SB2
  solver_parameters(3)=S2
  solver_parameters(4)=UU
  solver_parameters(5)=D


  v2_old=v(kvx)**2+v(kvy)**2+v(kvz)**2*g_cov(3)
  v2=rtsafe(solver8,v2_old,0.d0,1.d0-1.d-11,1.d-12,ierr)
  if((v2 .ge. 1) .or. (ierr>0)) then
    !write(*,*) abs(v(kvx)**2+v(kvy)**2+v(kvz)**2*g_cov(3)-v2)
    write(*,*) "Problems in cons2prim:"
    write(*,*) "v2 from rtsafe is larger than 1 or the routine failed: ",v2
    write(*,*) "Setting error flag to 778"
    errcode=778
    write(*,*) 'Error:',errcode
    run_crashed=.true.
    return  
  end if

  v(1:nv)=uin(1:nv)

  if(mhd .and. divclean) v(kglm)=u(kglm) !the artificial variable remains the same, both as primitive or conservative 
  
  w=w_algo_solver8_e_3p
  rh=D*sqrt(1.-v2)
  v(kpr)=0.25*(1-v2)*w

  v(krh)=rh

  vb=SB/w

  v(kvx)=(sx+vb*bx)/(w+B2)
  v(kvy)=(sy+vb*by)/(w+B2)
  v(kvz)=(sz/g_cov(3)+vb*bz)/(w+B2)

else if(algo_solver .eq. 6) then

  ierr=0
  exact9=.true.
 
  v2_old=v(kvx)**2+v(kvy)**2+v(kvz)**2

  oldsol9=4*v(kpr)/(1-v2_old)

  uin(:)=gm*u(:)

  DD=uin(krh)
  sx=uin(kvx)
  sy=uin(kvy)
  sz=uin(kvz)
  UU=uin(kpr)
  bx=uin(kbx)
  by=uin(kby)
  bz=uin(kbz)


  S2=sx*sx+sy*sy+sz*sz
  B2=bx*bx+by*by+bz*bz
  SB=(sx*bx+sy*by+sz*bz)
  SB2=SB**2

  !filled for the cases when the approximate solver is needed
  solver_parameters(1)=B2
  solver_parameters(2)=SB2
  solver_parameters(3)=SS2
  solver_parameters(4)=UU
  solver_parameters(5)=DD

  !we always assume a9=3
  b9=8*B2-4*UU
  c9=7*B2**2-8*B2*UU+S2
  d9=2*B2**3-4*B2**2*UU+2*B2*S2
  e9=SB2*B2

  p9=(24*c9-3*b9**2)/72
  q9=(b9**3-12*b9*c9+72*d9)/216

  dnull9=c9**2-3*b9*d9+36*e9
  deinz9=2*c9**3-9*b9*c9*d9+27*b9**2*e9+81*d9**2-216*c9*e9
  discr9=deinz9**2-4*dnull9**3

  if(discr9 .gt. 0) then
    droot9=sqrt(discr9)

    qbig9=(0.5*(deinz9+droot9))**(1./3.)
    

    drad9=-2*p9/3+(qbig9+dnull9/qbig9)/9

    if(drad9 .le. 0) then
      write(*,*) "Error in system_cons2prim, algo_solver9: drad9<0,",drad9,"no real solutions. I quit"
      call exit(12)
    end if
  
    sbig9=0.5*sqrt(drad9)
  else if(dnull9 .gt. 0) then
    phi9=acos(deinz9/(2*sqrt(dnull9**3)))
    sbig9=0.5*sqrt(-2.*p9/3+2./9.*sqrt(dnull9)*cos(phi9/3.)) 
  else 
    call solver10_approximate(b9,c9,d9,oldsol9,sol9) !this subroutine already checks that sol9>0
    exact9=.false.
    v2=get_x9(sol9,S2,B2,SB2)
    if((v2 .lt. 0) .or. (v2 .gt. 1) .or. (sol9 .lt. 0)) then
      write(*,*) "Sorry, but I got v^2 and (e+p)gamma^2 equal to:",v2,sol9," so the algo_solver9 failed. I quit with shame..."
      call exit(12)
    end if
  end if

  if(exact9) then

    rad129=-4*sbig9**2-2*p9+q9/sbig9
    rad349=-4*sbig9**2-2*p9-q9/sbig9

    rapp9=-b9/12

    if(rad129 .ge. 0) then
      root129=0.5*sqrt(rad129)
      sol1=rapp9-sbig9-root129
      sol2=rapp9-sbig9+root129
    end if

    if(rad349 .ge. 0) then
      root349=0.5*sqrt(rad349)
      sol3=rapp9+sbig9-root349
      sol4=rapp9+sbig9+root349
    end if

!    write(*,*) "Solutions:",sol1,sol2,sol3,sol4,oldsol9

    sol9=-1.
    dist9=1.d100

    if(sol2 .gt. 0) then
      v9=get_x9(sol2,S2,B2,SB2)
      if((v9 .ge. 0) .and. (v9 .lt. 1)) then
        sol9=sol2
        v2=v9
        dist9=abs(oldsol9-sol9)
        if(sol1 .gt. 0) then
          v9=get_x9(sol1,S2,B2,SB2)
          if((v9 .ge. 0) .and. (v9 .lt. 1)) then
            d19=abs(oldsol9-sol1)
            if(d19 .le. dist9) then
              sol9=sol1
              dist9=d19
              v2=v9
            end if
          end if
        end if
      end if
    end if

    if(sol4 .gt. 0) then
      v9=get_x9(sol4,S2,B2,SB2)
      if((v9 .ge. 0) .and. (v9 .lt. 1)) then
        d49=abs(oldsol9-sol4)
        if(d49 .le. dist9) then
          sol9=sol4
          dist9=d49
          v2=v9
        end if
      end if
      if(sol3 .gt. 0) then
        v9=get_x9(sol3,S2,B2,SB2)
        if((v9 .ge. 0) .and. (v9 .lt. 1)) then
          d39=abs(oldsol9-sol3)
          if(d39 .le. dist9) then
            sol9=sol3
            dist9=d39
            v2=v9
          end if
        end if
      end if
    end if

    if(sol9 .lt. 0) then
      write(*,*) "Sorry, in system_cons2prim the algo_solver9 failed... Quitting..."
      call exit(12)
    end if
  end if !end of the part using exact solutions

  v(1:nv)=uin(1:nv)
  
  v(krh)=DD*sqrt(1.-v2)

  v(kpr)=sol9/4*(1-v2)

  vb=SB/sol9

  v(kvx)=(sx+vb*bx)/(sol9+B2)
  v(kvy)=(sy+vb*by)/(sol9+B2)
  v(kvz)=(sz+vb*bz)/(sol9+B2)

  if(mhd .and. divclean) v(kglm)=u(kglm) !the artificial variable remains the same, both as primitive or conservative 

  else !no match with any solver
    write(*,*) "Sorry, but it not clear which solver you wish to use in the cons2prim subroutine. Exiting..."
    call exit(11)
  end if !end if on solver
  end if !end if on mhd case

end subroutine system_cons2prim

! *****************************************************************************

subroutine sist(N,X,F,K)
  implicit none
  real(8), dimension(*) :: X,F
  integer :: K, N
  real(8) :: B2, SB2, S2, U, D
  real(8) :: xx,yy,p

  xx=X(1)
  yy=X(2)

  B2=solver_parameters(1)
  SB2=solver_parameters(2)
  S2=solver_parameters(3)
  U=solver_parameters(4)
  D=solver_parameters(5)
  
  if(K .eq. 1) then
    F(1)=(yy+B2)**2*xx-SB2/(yy**2)*(2*yy+B2)-S2
  end if

  if(K .eq. 2) then
    p=(gamma_adindex-1)/gamma_adindex*((1.-xx)*yy-D*sqrt(1-xx))
    F(2)=yy-p+0.5*(1.+xx)*B2-0.5*SB2/(yy**2)-U
  end if

end subroutine sist

! *****************************************************************************

subroutine sist3(N,X,F,K)
  implicit none
  real(8), dimension(*) :: X,F
  integer :: K, N
  real(8) :: B2, SB2, S2, U, D
  real(8) :: xx,yy,pp,en,rh
  integer :: errcode !just for compatibility, it will not be used to avoid slowing down

  xx=X(1)
  yy=X(2)
  pp=X(3)

  B2=solver_parameters(1)
  SB2=solver_parameters(2)
  S2=solver_parameters(3)
  U=solver_parameters(4)
  D=solver_parameters(5)
  
  if(K .eq. 1) then
    F(1)=(yy+B2)**2*xx-SB2/(yy**2)*(2*yy+B2)-S2
  end if

  if(K .eq. 2) then
    F(2)=yy-pp+0.5*(1.+xx)*B2-0.5*SB2/(yy**2)-U
  end if

  if(K .eq. 3) then
    rh=D*sqrt(1-xx)
    call eos_energy(rh,en,pp,errcode)
    F(3)=yy-(pp+en)/(1-xx)
  end if

end subroutine sist3

! *****************************************************************************

subroutine sist4(N,X,F,K)
  implicit none
  real(8), dimension(*) :: X,F
  integer :: K, N
  real(8) :: B2, SB2, S2, U, D
  real(8) :: xx,yy,pp,en,rh
  integer :: errcode !just for compatibility, it will not be used to avoid slowing down

  yy=X(1)
  pp=X(2)

  B2=solver_parameters(1)
  SB2=solver_parameters(2)
  S2=solver_parameters(3)
  U=solver_parameters(4)
  D=solver_parameters(5)

  xx=(S2+SB2*(2.*yy+B2)/(yy**2))/((yy+B2)**2)
  
  xinner=xx
  
  if(K .eq. 2) then
    rh=D*sqrt(1-xx)
    call eos_energy(rh,en,pp,errcode)
    F(2)=yy-(pp+en)/(1-xx)
  end if

  if(K .eq. 1) then
    F(1)=yy-pp+0.5*(1.+xx)*B2-0.5*SB2/(yy**2)-U
  end if

end subroutine sist4

! *****************************************************************************

subroutine sist5(N,X,F,fjac3,ldjac,iflag)
  implicit none
  real(8), dimension(3) :: X,F
  real(8), dimension(3,3) :: fjac3
  integer :: iflag, ldjac, N
  real(8) :: B2, SB2, S2, U, D
  real(8) :: xx,yy,pp,en,rh,dendpr,dendrh
  integer :: errcode !just for compatibility, it will not be used to avoid slowing down

  !this procedure assume en=en(pr), there are some other terms if en=en(rho,pr)

  xx=X(1)
  yy=X(2)
  pp=X(3)

  B2=solver_parameters(1)
  SB2=solver_parameters(2)
  S2=solver_parameters(3)
  U=solver_parameters(4)
  D=solver_parameters(5)

  if(iflag .eq. 1) then
    F(1)=(yy+B2)**2*xx-SB2/(yy**2)*(2*yy+B2)-S2
    F(2)=yy-pp+0.5*(1.+xx)*B2-0.5*SB2/(yy**2)-U
    rh=D*sqrt(1-xx)
    call eos_energy(rh,en,pp,errcode)
    F(3)=yy-(pp+en)/(1-xx)
  end if
  if(iflag .eq. 2) then
    rh=D*sqrt(1-xx)
    call eos_energy_with_ders(rh,en,pp,dendrh,dendpr,errcode)
    fjac3(1,1)=(yy+B2)**2
    fjac3(1,2)=2*(yy+B2)*xx+(2*SB2*(2*yy+B2))/(yy**3)-2*SB2/(yy**2)
    fjac3(1,3)=0.
    fjac3(2,1)=B2/2.
    fjac3(2,2)=1.+SB2/(yy**3)
    fjac3(2,3)=-1.
    fjac3(3,1)=-(pp+en)/((1-xx)**2)
    fjac3(3,2)=1.
    fjac3(3,3)=-(1+dendpr)/(1-xx)
  end if

end subroutine sist5

! *****************************************************************************

subroutine system_flux(v,u,f,vf,dir,errcode)
!-- Calculate fluid fluxes and characteristic velocities
  use common, only: viscosity
  
  implicit none

  real(8),dimension(nv), intent(inout) :: v 
  real(8),intent(out),dimension(nv) :: u,f
  real(8),intent(out),dimension(2 ) :: vf
  integer,intent(in) :: dir
  real(8) :: gcov_dir  !value of the component of metric tensor needed for computing characteristic speeds
  real(8) :: B2,E2,EEM
  real(8), dimension(1:3) :: E_field_cova, EXB 
  real(8) :: rh,vx,vy,vz,pr,en,v2,glf,w,ww,wwdir,cs2,a2,den,vf1,vf2,vdir,ca2,bb2
  integer, intent(out) :: errcode

  real(8) :: tvx, tey, tbz, tez, tby, tvz,tvy, tbx, tex


  errcode=0

  rh=v(krh)
  vx=v(kvx)
  vy=v(kvy)
  vz=v(kvz)
  pr=v(kpr)
 
  call eos_energy(rh,en,pr,errcode)
  if (errcode .gt. 0) then
     write(*,*) "An error occurred into system_flux when computing energy density"
     run_crashed=.true.
     return
  end if

  v2=g_cov(1)*vx*vx+g_cov(2)*vy*vy+g_cov(3)*vz*vz
  if (v2 .gt. 1) then
     write(*,*) "Superluminal velocity found in system_flux"
     write(*,*) "Value:", v2
     errcode=333
     run_crashed=.true.
     return
  end if
  glf=1./sqrt(1.-v2)
  if(.not. viscosity) then
    w=en+pr
    if(mhd) then
      ! remind that here we use E covariant and B contravariant
      if(rmhd) then
        E_field_cova(1:2)=v(kex:key)
        E_field_cova(3)=v(kez)
        E2=E_field_cova(1)*v(kex)+E_field_cova(2)*v(key)+E_field_cova(3)*v(kez)
      else
        call system_compute_E_in_ideal_MHD(v(kvx:kvz),v(kbx:kbz),E_field_cova)
        E2=E_field_cova(1)**2+E_field_cova(2)**2+E_field_cova(3)**2/g_cov(3)
      end if
      call system_EXB(E_field_cova,v(kbx:kbz),EXB)
      B2=v(kbx)**2+v(kby)**2+v(kbz)**2*g_cov(3)
      EEM=0.5d0*(E2+B2)

    end if
  else
    w=en+pr+v(kpibu)
  end if

  ww=w*glf**2.

   u(krh)=rh*glf

  if(viscosity) then
 
   if(obtained .eq. 'no') then
     call get_derived_pi(v,pitt,pitx,pity,pitz)
   else
     call get_derived_pi_zz(v,pitt,pitx,pity,pitz,pizz)
     v(kpizz)=pizz
   end if
   u(kpibu:)=v(kpibu:)*u(krh)


   u(kvx)=ww*vx+pitx
   u(kvy)=ww*vy+pity
   u(kvz)=(ww*vz+pitz)*g_cov(3)
   u(kpr)=ww-pr-v(kpibu)-pitt*g_cov0

  else 
   if(mhd) then
       
     u(kvx)=ww*vx+EXB(1)
     u(kvy)=ww*vy+EXB(2)
     u(kvz)=ww*vz*g_cov(3)+EXB(3)
     u(kpr)=ww-pr+EEM
     if(rmhd) then
       u(kbx:kez)=v(kbx:kez)
     else
       u(kbx:kbz)=v(kbx:kbz)
     end if
     if(divclean) u(kglm)=v(kglm)/gp !to counterbalance the multiplication by gp at the end

   else
     u(kvx)=ww*vx
     u(kvy)=ww*vy
     u(kvz)=ww*vz*g_cov(3)
     u(kpr)=ww-pr
   end if
  end if

  
!-- Fluxes

select case(dir)
case(1)

  wwdir=ww*vx
  vdir=vx
  gcov_dir=g_cov(1)
  
  f(krh)=rh*glf*vx
  if(viscosity) then
    f(kvx)=wwdir*vx*g_cov(1)+pr+v(kpibu)+v(kpixx)*g_cov(1)
    f(kvy)=wwdir*vy*g_cov(2)+v(kpixy)*g_cov(2)
    f(kvz)=wwdir*vz*g_cov(3)+v(kpixz)*g_cov(3)
    f(kpr)=wwdir+pitx
    f(nou+1:nv)=u(nou+1:nv)*vx
  else
    if(mhd) then
      f(kvx)=wwdir*vx-E_field_cova(1)**2-v(kbx)**2+pr+EEM !(e+p)g^2*v^x*v_x-E_x*E^x-B_x*B^x+p+EEM
      f(kvy)=wwdir*vy-E_field_cova(1)*E_field_cova(2)-v(kbx)*v(kby) !(e+p)g^2*v^x*v_y-E_y*E^x-B_y*B^x
      f(kvz)=wwdir*vz*g_cov(3)-E_field_cova(1)*E_field_cova(3)-v(kbx)*v(kbz)*g_cov(3) !(e+p)g^2*v^x*v_z-E_z*E^x-B_z*B^x
      f(kpr)=wwdir+EXB(1)

      f(kbx)=0.
      f(kby)=-E_field_cova(3)*gm
      f(kbz)=+E_field_cova(2)*gm
      if(rmhd) then
        f(kex)=0.
        f(key)=+v(kbz)*g_cov(3)*gm
        f(kez)=-v(kby)*g_cov(2)*gm
      end if
      if(divclean) then
        f(kbx)=v(kglm)
        f(kglm)=v(kbx)/gp !for the non-physical variable we do not use the sqrt(-g) factor, so we counter-balance the final *gp
      end if

    else
      f(kvx)=wwdir*vx*g_cov(1)+pr
      f(kvy)=wwdir*vy*g_cov(2)
      f(kvz)=wwdir*vz*g_cov(3) 
      f(kpr)=wwdir
    end if
  end if

case(2)

  wwdir=ww*vy
  vdir=vy
  gcov_dir=g_cov(2)
  
  f(krh)=rh*glf*vy
  if(viscosity) then
    f(kvx)=wwdir*vx*g_cov(1)+v(kpixy)*g_cov(1)
    f(kvy)=wwdir*vy*g_cov(2)+pr+v(kpibu)+v(kpiyy)*g_cov(2)
    f(kvz)=wwdir*vz*g_cov(3)+v(kpiyz)*g_cov(3)
    f(kpr)=wwdir+pity
    f(nou+1:nv)=u(nou+1:nv)*vy
  else
   if(mhd) then
      f(kvx)=wwdir*vx-E_field_cova(2)*E_field_cova(1)-v(kby)*v(kbx)!(e+p)g^2*v^y*v_x-E_x*E^y-B_x*B^y
      f(kvy)=wwdir*vy-E_field_cova(2)**2-v(kby)**2+pr+EEM !(e+p)g^2*v^y*v_y-E_y*E^y-B_y*B^y
      f(kvz)=wwdir*vz*g_cov(3)-E_field_cova(2)*E_field_cova(3)-v(kby)*v(kbz)*g_cov(3) !(e+p)g^2*v^y*v_z-E_z*E^y-B_z*B^y
      f(kpr)=wwdir+EXB(2)
      
      f(kbx)=E_field_cova(3)*gm
      f(kby)=0.
      f(kbz)=-E_field_cova(1)*gm
      if(rmhd) then
        f(kex)=-v(kbz)*g_cov(3)*gm
        f(key)=0.
        f(kez)=+v(kbx)*g_cov(1)*gm
      end if
      if(divclean) then
        f(kby)=v(kglm)
        f(kglm)=v(kby)/gp !for the non-physical variable we do not use the sqrt(-g) factor, so we counter-balance the final *gp
      end if

    else
      f(kvx)=wwdir*vx*g_cov(1)
      f(kvy)=wwdir*vy*g_cov(2)+pr
      f(kvz)=wwdir*vz*g_cov(3) 
      f(kpr)=wwdir
    end if
  end if

case(3)

  wwdir=ww*vz
  vdir=vz
  gcov_dir=g_cov(3)
  
  f(krh)=rh*glf*vz
  if(viscosity) then
    f(kvx)=wwdir*vx*g_cov(1)+v(kpixz)*g_cov(1)
    f(kvy)=wwdir*vy*g_cov(2)+v(kpiyz)*g_cov(2)
    f(kvz)=wwdir*vz*g_cov(3)+pr+v(kpibu)+v(kpizz)*g_cov(3)
    f(kpr)=wwdir+pitz
    f(nou+1:nv)=u(nou+1:nv)*vz
  else
    if(mhd) then
      f(kvx)=wwdir*vx-E_field_cova(3)/g_cov(3)*E_field_cova(1)-v(kbz)*v(kbx) !(e+p)g^2*v^z*v_x-E_x*E^z-B_x*B^z
      f(kvy)=wwdir*vy-E_field_cova(3)/g_cov(3)*E_field_cova(2)-v(kbz)*v(kby) !(e+p)g^2*v^z*v_y-E_y*E^z-B_y*B^z
      f(kvz)=wwdir*vz*g_cov(3)-E_field_cova(3)**2/g_cov(3)-v(kbz)**2*g_cov(3)+pr+EEM !(e+p)g^2*v^z*v_z-E_z*E^z-B_z*B^z
      f(kpr)=wwdir+EXB(3)/g_cov(3)

      f(kbx)=-E_field_cova(2)*gm
      f(kby)=E_field_cova(1)*gm
      f(kbz)=0.
      if(rmhd) then
        f(kex)=+v(kby)*g_cov(2)*gm
        f(key)=-v(kbx)*g_cov(1)*gm
        f(kez)=0.
      end if
      if(divclean) then
        f(kbz)=v(kglm)/g_cov(3)
        f(kglm)=v(kbz)/gp !for the non-physical variable we do not use the sqrt(-g) factor, so we counter-balance the final *gp
      end if


    else
      f(kvx)=wwdir*vx*g_cov(1)
      f(kvy)=wwdir*vy*g_cov(2)
      f(kvz)=wwdir*vz*g_cov(3)+pr 
      f(kpr)=wwdir
    end if
  end if
end select


  u(1:nu)=u(1:nu)*gp
  f(1:nu)=f(1:nu)*gp

!-- Fast speeds

  if(mhd) then
    bb2=B2-E2
    ca2=bb2/(bb2+(en+pr))
  else
    ca2=0.
  endif

  call eos_sound(rh,en,pr,cs2,errcode)
  if(errcode .gt. 0) then
    write(*,*) "An error occurred when computing the speed of sound into the system_flux subroutine."
    write(*,*) "Error code:", errcode
    run_crashed=.true.
    return
  end if

    a2=cs2+ca2-cs2*ca2

    den=1./(1.-v2*a2)

    vf1=den*(1.-a2)*vdir
    vf2=den*sqrt(a2*(1.-v2)*((1.-v2*a2)/gcov_dir-(1.-a2)*vdir**2))

    vf(1)=vf1+vf2
    vf(2)=vf1-vf2


end subroutine system_flux

! *****************************************************************************

subroutine system_source(ix,iy,iz,vars,v,src,src_stiff,errcode)
!-- Calculate source terms
  use common, only: nv,nu,nov,nou,krh,kvx,kvy,kvz,kpr,kpibu,kpizz,kpixz,kpixy,kpiyz,kpixx,kpiyy
  use common, only: coordinates,g_cov,gp,gm,g_cov0,viscosity,z,x,y,pitz,pitt,pitx,pity
  use common, only: MINKOWSKI, BJORKEN
  use viscous

  implicit none
  integer,intent(in) :: ix, iy, iz
  real(8),intent(inout ),dimension(nv) :: v
  real(8),intent(out),dimension(nv) :: src,src_stiff
  real(8),intent(in),dimension(:,:,:,:), allocatable :: vars

  real(8) :: rh,vx,vy,vz,pr,en,v2,glf,w,ww
  real(8) :: dutdt, dutdx, dutdy, dutdz, duxdx, duydy, duzdz, duzdt, duzdx, duzdy, duxdt, duxdy, duxdz, duydt, duydx, duydz
  real(8) :: sum_der ! sum of derivatives
  real(8) :: sum_der_t ! sum of derivatives + ut/t
  real(8) :: ut
  real(8) :: sum_der_0, sum_der_x, sum_der_y, sum_der_z
  real(8) :: w0x, w0y, w0z, wxy, wxz, wyz, wyx, wzx, wzy, wx0, wy0, wz0
  real(8), parameter :: lambda0=1. !this controls the -4/3 terms
  real(8), parameter :: lambda1=1., lambda2=0.
  
  real(8) :: IV00, IV0x, IV0y, IV0z, IVxx, IVyy, IVzz, IVxy, IVxz, IVyz
  
  real(8) :: B2,E2,EEM, rhoC, vdotE, vdotB
  real(8), dimension(1:3) :: E_field_cova, EXB, vXB, vXE
  
  integer :: errcode
  

  src(1:nv)=0.
  src_stiff(1:nv)=0.
 
 
    rh=v(krh)
    vx=v(kvx)
    vy=v(kvy)
    vz=v(kvz)
    pr=v(kpr)

    call eos_energy(rh,en,pr,errcode)
    if(errcode .gt. 0) then
      write(*,*) "Error into the out_vars subroutine when trying to get derived quantities (temp, entr_dens, en_dens)"
      write(*,*) "Error code:", errcode
      write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
      run_crashed=.true.
      return
    end if


    v2=g_cov(1)*vx*vx+g_cov(2)*vy*vy+g_cov(3)*vz*vz
    glf=1./sqrt(1.-v2)
    ut=glf

  if (coordinates .eq. BJORKEN) then !BJORKEN COORDINATES
!-- Only surviving term: -0.5*gp*W^(xx)*d(h_x^2)/dtau=-tau^2*W^(xx)=-W^x_x
    if(.not. viscosity) then
      w=en+pr
      ww=w*glf**2
      if(mhd) then
        call system_compute_E_in_ideal_MHD(v(kvx:kvz),v(kbx:kbz),E_field_cova)
        B2=v(kbx)**2+v(kby)**2+v(kbz)**2*g_cov(3)
        E2=E_field_cova(1)**2+E_field_cova(2)**2+E_field_cova(3)**2/g_cov(3)
        EEM=0.5*(E2+B2)
        src(kpr)=-t*(ww*vz*vz-(E_field_cova(3)/g_cov(3))**2-v(kbz)*v(kbz)+(EEM+pr)/g_cov(3)) 
        if(rmhd) then
          call system_compute_divE(ix,iy,iz,rhoC)
          call system_crossprod(v(kvx:kvz),v(kbx:kbz),vXB)
          call system_crossprod(v(kvx:kvz),v(kex:kez),vXE)
          vdotB=v(kvx)*v(kbx)+v(kvy)*v(kby)+v(kvz)*v(kbz)*g_cov(3)
          vdotE=v(kvx)*v(kex)+v(kvy)*v(key)+v(kvz)*v(kez)*g_cov(3)
          src(kex:kez)=-rhoC*v(kvx:kvz)
          src_stiff(kex:kez)=-glf*(sigma_el_cond*(v(kex:kez)-vdotE*v(kvx:kvz)+vXB(:))+sigma_chiral_cond*(v(kbx:kbz)-vdotB*&
                            &v(kvx:kvz)-vXE(:)))
        end if
      else
        src(kpr)=-t*(ww*vz*vz+pr/g_cov(3)) 
      endif
    else
      call get_sigma(ix,iy,iz,vars,sigma,dutdt, dutdx, dutdy, dutdz, duxdx, duydy, duzdz, duzdt, duzdx, duzdy,&
                  &duxdt, duxdy, duxdz, duydt, duydx, duydz)
      call get_viscous_parameters(eta_vis,zeta_vis,tau_pi,tau_pi_big,v(krh),v(kpr),errcode)
      if(errcode .gt. 0) then
        write(*,*) "An error occurred into the system_source subroutine after trying to get viscous parameters."
        write(*,*) "Error code:", errcode
        write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
        run_crashed=.true.
        return
      end if
    
      if(obtained .eq. 'no') then
        call get_derived_pi(v,pitt,pitx,pity,pitz)
      else
        call get_derived_pi_zz(v,pitt,pitx,pity,pitz,pizz)
        v(kpizz)=pizz
      end if
    
      sum_der = dutdt+duxdx+duydy+duzdz
      sum_der_t = sum_der + ut/t
    
      sum_der_0=g_cov0*(glf*dutdt+glf*v(kvx)*dutdx+glf*v(kvy)*dutdy+glf*v(kvz)*dutdz)
      sum_der_x=g_cov(1)*(glf*duxdt+glf*v(kvx)*duxdx+glf*v(kvy)*duxdy+glf*v(kvz)*duxdz)
      sum_der_y=g_cov(2)*(glf*duydt+glf*v(kvx)*duydx+glf*v(kvy)*duydy+glf*v(kvz)*duydz)
      sum_der_z=glf*(t**2.*duzdt+2.*t*v(kvz))+(t**2.)*(glf*v(kvx)*duzdx+glf*v(kvy)*duzdy+glf*v(kvz)*duzdz)

      w=en+pr+v(kpibu)
      ww=w*glf**2.
    
      src(kpr)=-t*(ww*vz*vz+(pr+v(kpibu))/g_cov(3)+v(kpizz)) !the t term depends on the christoffel symbol
      src(kpibu)=-1./tau_pi_big*(v(kpibu)+zeta_vis*sum_der_t)-lambda0*4./3.*v(kpibu)*sum_der_t!DDDAAA 4/3 to be checked 
     
      w0x=-0.5*(dutdx+duxdt)+0.5*glf*(sum_der_x-v(kvx)*g_cov(1)*sum_der_0)
      w0y=-0.5*(dutdy+duydt)+0.5*glf*(sum_der_y-v(kvy)*g_cov(2)*sum_der_0)
      w0z=-0.5*(dutdz+duzdt)+0.5*glf*(sum_der_z-v(kvz)*g_cov(3)*sum_der_0)
      wxy=0.5*(g_cov(2)*duydx-duxdy)+0.5*glf*(g_cov(1)*v(kvx)*sum_der_y-g_cov(2)*v(kvy)*sum_der_x)
      wxz=0.5*(g_cov(3)*duzdx-duxdz)+0.5*glf*(g_cov(1)*v(kvx)*sum_der_z-g_cov(3)*v(kvz)*sum_der_x)
      wyz=0.5*(g_cov(3)*duzdy-duydz)+0.5*glf*(g_cov(2)*v(kvy)*sum_der_z-g_cov(3)*v(kvz)*sum_der_y)
      wx0=w0x
      wy0=w0y
      wz0=w0z
      wyx=-wxy
      wzx=-wxz
      wzy=-wyz
      
      IV00=2.*(pitx*w0x+pity*w0y+pitz*w0z)
      IV0x=(pitt*wx0+pity*wxy+v(kpixx)*w0x+v(kpixy)*w0y+v(kpixz)*w0z+pitz*wxz)
      IV0y=(pitt*wy0+pitx*wyx+v(kpiyy)*w0y+v(kpixy)*w0x+v(kpiyz)*w0z+pitz*wyz)
      IV0z=(pitt*wz0+pitx*wzx+v(kpizz)*w0z+v(kpixz)*w0x+v(kpiyz)*w0y+pity*wzy)
      IVxx=2.*(pitx*wx0+v(kpixy)*wxy+v(kpixz)*wxz)
      IVyy=2.*(pity*wy0+v(kpixy)*wyx+v(kpiyz)*wyz)
      IVzz=2.*(pitz*wz0+v(kpixz)*wzx+v(kpiyz)*wzy)
      IVxy=v(kpixx)*wyx+v(kpiyy)*wxy+v(kpixz)*wyz+v(kpiyz)*wxz+pitx*wy0+pity*wx0
      IVxz=v(kpixx)*wzx+v(kpizz)*wxz+v(kpixy)*wzy+v(kpiyz)*wxy+pitx*wz0+pitz*wx0
      IVyz=v(kpiyy)*wzy+v(kpizz)*wyz+v(kpixy)*wzx+v(kpixz)*wyx+pity*wz0+pitz*wy0
               
      src(kpizz)=-1./tau_pi*(v(kpizz)+2.*eta_vis*sigma(szz))-2.*glf*v(kpizz)/t-&
              &2.*(glf*vz)*pitz/t-lambda0*4./3.*v(kpizz)*sum_der_t+lambda1*glf*(2.*v(kvz)*(pitz*sum_der_0+&
              &v(kpixz)*sum_der_x+v(kpiyz)*sum_der_y+v(kpizz)*sum_der_z))-lambda2*IVzz

      src(kpixz)=-1./tau_pi*(v(kpixz)+2.*eta_vis*sigma(sxz))+lambda0*(-4./3.*v(kpixz)*sum_der_t)-glf/t*v(kpixz)-glf*vz*pitx/t+&
              &lambda1*glf*((pitx*v(kvz)+pitz*v(kvx))*sum_der_0+(v(kpixx)*v(kvz)+v(kpixz)*v(kvx))*sum_der_x+(v(kpixy)*v(kvz)+&
              &v(kpiyz)*v(kvx))*sum_der_y+(v(kpixz)*v(kvz)+v(kpizz)*v(kvx))*sum_der_z)-lambda2*IVxz

      src(kpiyz)=-1./tau_pi*(v(kpiyz)+2.*eta_vis*sigma(syz))+lambda0*(-4./3.*v(kpiyz)*sum_der_t)-glf/t*v(kpiyz)-glf*vz*pity/t+&
              &lambda1*glf*((pity*v(kvz)+pitz*v(kvy))*sum_der_0+(v(kpixy)*v(kvz)+v(kpixz)*v(kvy))*sum_der_x+(v(kpiyy)*v(kvz)+&
              &v(kpiyz)*v(kvy))*sum_der_y+(v(kpiyz)*v(kvz)+v(kpizz)*v(kvy))*sum_der_z)-lambda2*IVyz

      src(kpixy)=-1./tau_pi*(v(kpixy)+2.*eta_vis*sigma(sxy))+lambda0*(-4./3.*v(kpixy)*sum_der_t)+lambda1*glf*((pitx*&
              &v(kvy)+pity*v(kvx))*sum_der_0+(v(kpixx)*v(kvy)+v(kpixy)*v(kvx))*sum_der_x+(v(kpixy)*v(kvy)+v(kpiyy)*v(kvx))*&
              &sum_der_y+(v(kpixz)*v(kvy)+v(kpiyz)*v(kvx))*sum_der_z)-lambda2*IVxy

      src(kpixx)=-1./tau_pi*(v(kpixx)+2.*eta_vis*sigma(sxx))+lambda0*(-4./3.*v(kpixx)*sum_der_t)+lambda1*glf*(2.*v(kvx)*(pitx*&
               &sum_der_0+v(kpixx)*sum_der_x+v(kpixy)*sum_der_y+v(kpixz)*sum_der_z))-lambda2*IVxx
    
      src(kpiyy)=-1./tau_pi*(v(kpiyy)+2.*eta_vis*sigma(syy))+lambda0*(-4./3.*v(kpiyy)*sum_der_t)+lambda1*glf*(2.*v(kvy)*(pity*&
               &sum_der_0+v(kpixy)*sum_der_x+v(kpiyy)*sum_der_y+v(kpiyz)*sum_der_z))-lambda2*IVyy

    end if
  else
    ! MINKOWSKI COORDINATES 
    if(rmhd) then
       call system_compute_divE(ix,iy,iz,rhoC)
       call system_crossprod(v(kvx:kvz),v(kbx:kbz),vXB)
       call system_crossprod(v(kvx:kvz),v(kex:kez),vXE)
       vdotE=v(kvx)*v(kbx)+v(kvy)*v(kby)+v(kvz)*v(kbz)*g_cov(3)
       vdotE=v(kvx)*v(kex)+v(kvy)*v(key)+v(kvz)*v(kez)*g_cov(3)
       src(kex:kez)=-rhoC*v(kvx:kvz)
       src_stiff(kex:kez)=-glf*(sigma_el_cond*(v(kex:kez)-vdotE*v(kvx:kvz)+vXB(:))+sigma_chiral_cond*(v(kbx:kbz)-vdotB*&
                         &v(kvx:kvz)-vXE(:)))
    end if
      
    if (viscosity) then
      call get_sigma(ix,iy,iz,vars,sigma,dutdt, dutdx, dutdy, dutdz, duxdx,duydy, duzdz, duzdt, duzdx, duzdy,&
                 &duxdt, duxdy, duxdz, duydt, duydx, duydz)
      call get_viscous_parameters(eta_vis,zeta_vis,tau_pi,tau_pi_big,v(krh),v(kpr),errcode)
      if(errcode .gt. 0) then
       write(*,*) "An error occurred into the system_source subroutine after trying to get viscous parameters."
       write(*,*) "Error code:", errcode
       write(*,*) "Position on the grid: ix=",ix,"iy=",iy,"iz=",iz
       run_crashed=.true.
       return
      end if

      if(obtained .eq. 'no') then
        call get_derived_pi(v,pitt,pitx,pity,pitz)
      else
        call get_derived_pi_zz(v,pitt,pitx,pity,pitz,pizz)
        v(kpizz)=pizz
      end if
      sum_der = dutdt+duxdx+duydy+duzdz
      !AAA think about the possibility of removing the metric tensor component, since we are assuming to know what they are
      sum_der_0=g_cov0*(glf*dutdt+v(kvx)*dutdx+v(kvy)*dutdy+v(kvz)*dutdz)
      sum_der_x=g_cov(1)*(glf*duxdt+v(kvx)*duxdx+v(kvy)*duxdy+v(kvz)*duxdz)
      sum_der_y=g_cov(2)*(glf*duydt+v(kvx)*duydx+v(kvy)*duydy+v(kvz)*duydz)
      sum_der_z=g_cov(3)*(glf*duzdt+v(kvx)*duzdx+v(kvy)*duzdy+v(kvz)*duzdz)



    
!   still Minkowski coordinates
    src(kpibu)=-1./tau_pi_big*(v(kpibu)+zeta_vis*sum_der)+lambda0*(-1./3.*v(kpibu)*sum_der)
     
    w0x=-0.5*(dutdx+duxdt)+0.5*glf*(sum_der_x-v(kvx)*g_cov(1)*sum_der_0)
    w0y=-0.5*(dutdy+duydt)+0.5*glf*(sum_der_y-v(kvy)*g_cov(2)*sum_der_0)
    w0z=-0.5*(dutdz+duzdt)+0.5*glf*(sum_der_z-v(kvz)*g_cov(3)*sum_der_0)
    wxy=0.5*(g_cov(2)*duydx-duxdy)+0.5*glf*(g_cov(1)*v(kvx)*sum_der_y-g_cov(2)*v(kvy)*sum_der_x)
    wxz=0.5*(g_cov(3)*duzdx-duxdz)+0.5*glf*(g_cov(1)*v(kvx)*sum_der_z-g_cov(3)*v(kvz)*sum_der_x)
    wyz=0.5*(g_cov(3)*duzdy-duydz)+0.5*glf*(g_cov(2)*v(kvy)*sum_der_z-g_cov(3)*v(kvz)*sum_der_y)
    wx0=w0x
    wy0=w0y
    wz0=w0z
    wyx=-wxy
    wzx=-wxz
    wzy=-wyz
     
    IV00=2.*(pitx*w0x+pity*w0y+pitz*w0z)
    IV0x=(pitt*wx0+pity*wxy+v(kpixx)*w0x+v(kpixy)*w0y+v(kpixz)*w0z+pitz*wxz)
    IV0y=(pitt*wy0+pitx*wyx+v(kpiyy)*w0y+v(kpixy)*w0x+v(kpiyz)*w0z+pitz*wyz)
    IV0z=(pitt*wz0+pitx*wzx+v(kpizz)*w0z+v(kpixz)*w0x+v(kpiyz)*w0y+pity*wzy)
    IVxx=2.*(pitx*wx0+v(kpixy)*wxy+v(kpixz)*wxz)
    IVyy=2.*(pity*wy0+v(kpixy)*wyx+v(kpiyz)*wyz)
    IVzz=2.*(pitz*wz0+v(kpixz)*wzx+v(kpiyz)*wzy)
    IVxy=v(kpixx)*wyx+v(kpiyy)*wxy+v(kpixz)*wyz+v(kpiyz)*wxz+pitx*wy0+pity*wx0
    IVxz=v(kpixx)*wzx+v(kpizz)*wxz+v(kpixy)*wzy+v(kpiyz)*wxy+pitx*wz0+pitz*wx0
    IVyz=v(kpiyy)*wzy+v(kpizz)*wyz+v(kpixy)*wzx+v(kpixz)*wyx+pity*wz0+pitz*wy0
     
    src(kpizz)=-1./tau_pi*(v(kpizz)+2.*eta_vis*sigma(szz))+lambda0*(-4./3.*v(kpizz)*sum_der)+lambda1*(2.*v(kvz)*(pitz*sum_der_0&
              &+v(kpixz)*sum_der_x+v(kpiyz)*sum_der_y+v(kpizz)*sum_der_z))-lambda2*IVzz

    src(kpixz)=-1./tau_pi*(v(kpixz)+2.*eta_vis*sigma(sxz))+lambda0*(-4./3.*v(kpixz)*sum_der)+lambda1*((pitx*v(kvz)+pitz*v(kvx))*&
              &sum_der_0+(v(kpixx)*v(kvz)+v(kpixz)*v(kvx))*sum_der_x+(v(kpixy)*v(kvz)+v(kpiyz)*v(kvx))*sum_der_y+(v(kpixz)*v(kvz)&
              &+v(kpizz)*v(kvx))*sum_der_z)-lambda2*IVxz

    src(kpixy)=-1./tau_pi*(v(kpixy)+2.*eta_vis*sigma(sxy))+lambda0*(-4./3.*v(kpixy)*sum_der)+lambda1*((pitx*v(kvy)+pity*v(kvx))*&
              &sum_der_0+(v(kpixx)*v(kvy)+v(kpixy)*v(kvx))*sum_der_x+(v(kpixy)*v(kvy)+v(kpiyy)*v(kvx))*sum_der_y+(v(kpixz)*v(kvy)&
               &+v(kpiyz)*v(kvx))*sum_der_z)-lambda2*IVxy

    src(kpiyz)=-1./tau_pi*(v(kpiyz)+2.*eta_vis*sigma(syz))+lambda0*(-4./3.*v(kpiyz)*sum_der)+lambda1*((pity*v(kvz)+pitz*v(kvy))*&
              &sum_der_0+(v(kpixy)*v(kvz)+v(kpixz)*v(kvy))*sum_der_x+(v(kpiyy)*v(kvz)+v(kpiyz)*v(kvy))*sum_der_y+(v(kpiyz)*v(kvz)&
               &+v(kpizz)*v(kvy))*sum_der_z)-lambda2*IVyz

    src(kpixx)=-1./tau_pi*(v(kpixx)+2.*eta_vis*sigma(sxx))+lambda0*(-4./3.*v(kpixx)*sum_der)+lambda1*(2.*v(kvx)*(pitx*sum_der_0+&
              &v(kpixx)*sum_der_x+v(kpixy)*sum_der_y+v(kpixz)*sum_der_z))-lambda2*IVxx

    src(kpiyy)=-1./tau_pi*(v(kpiyy)+2.*eta_vis*sigma(syy))+lambda0*(-4./3.*v(kpiyy)*sum_der)+lambda1*(2.*v(kvy)*(pity*sum_der_0+&
              &v(kpixy)*sum_der_x+v(kpiyy)*sum_der_y+v(kpiyz)*sum_der_z))-lambda2*IVyy
   
   end if !end if viscosity 
  end if !end if Minkwoski coordinates
  if(viscosity) then 
    src(kpibu:)=src(kpibu:)*v(krh)
  end if
  src=src*gp

end subroutine system_source

! *****************************************************************************

subroutine system_EXB(arrE, arrB, arrEXB)
!it computes the cross product E x B
!E should be given covariant, B contravariant

implicit none
real(8), intent(in), dimension(1:3) :: arrE, arrB
real(8), intent(out), dimension(1:3) :: arrEXB
real(8) :: E_contra(3)

E_contra(1:2)=arrE(1:2)
E_contra(3)=arrE(3)/g_cov(3)

arrEXB(1)=gp*(E_contra(2)*arrB(3)-E_contra(3)*arrB(2))
arrEXB(2)=gp*(E_contra(3)*arrB(1)-E_contra(1)*arrB(3))   
arrEXB(3)=gp*(E_contra(1)*arrB(2)-E_contra(2)*arrB(1))

end subroutine system_EXB


! *****************************************************************************

subroutine system_crossprod(arrX, arrY, cross)
!it computes the cross product X x Y
!both X and Y should be covariant

implicit none
real(8), intent(in), dimension(1:3) :: arrX, arrY
real(8), intent(out), dimension(1:3) :: cross


cross(1)=gm*(arrX(2)*arrY(3)-arrX(3)*arrY(2))
cross(2)=gm*(arrX(3)*arrY(1)-arrX(1)*arrY(3))   
cross(3)=gm*(arrX(1)*arrY(2)-arrX(2)*arrY(1))

end subroutine system_crossprod


! *****************************************************************************

subroutine system_compute_E_in_ideal_MHD(arrV, arrB, arrE)
!it computes the electric field arrE given the arrays arrV (velocities) and arrB (magnetic fields)
!it returns the covariant E field and it takes in input the contrav. V and B 

implicit none
real(8), intent(in),dimension(1:3) :: arrV, arrB
real(8), intent(out),dimension(1:3) :: arrE

arrE(1)=-gp*(arrV(2)*arrB(3)-arrV(3)*arrB(2))
arrE(2)=-gp*(arrV(3)*arrB(1)-arrV(1)*arrB(3))
arrE(3)=-gp*(arrV(1)*arrB(2)-arrV(2)*arrB(1))

end subroutine system_compute_E_in_ideal_MHD

! *****************************************************************************

subroutine system_compute_divE(i, j, k, divEcontra)
!it computes the divergence of E contravariant

use common, only: w,ddx,ddy,ddz
implicit none
integer, intent(in) :: i,j,k
real(8), intent(out) :: divEcontra
real(8) :: dx,dy,dz

!we assume that the grid is equally spaced
!in this second order derivative, the delta_x is two cells wide
dx=2.d0/ddx(i)
dy=2.d0/ddy(j)
dz=2.d0/ddz(k)

divEcontra=(w(i+1,j,k,kex)-w(i-1,j,k,kex))/dx+(w(i,j+1,k,key)-w(i,j-1,k,key))/dy+(w(i,j,k+1,kez)-w(i,j,k-1,kez))/dz

end subroutine system_compute_divE

! *****************************************************************************

! personal routines for 3x3 sistem of equations solving

      subroutine ludecompose(A, n, L, P, smallvalue)
      !A is destroyed and transformed into an upper triangular matrix
      implicit none
      real(8), dimension(1:n,1:n) :: A
      real(8), dimension(1:n,1:n) :: L
      real(8), dimension(1:n) :: tmprow
      integer, dimension(1:n) :: P
      integer, intent(in) :: n
      integer :: i,j,k, indP, tmpindx
      real(8) :: m, piv, smallvalue, tmpreal
      integer :: allocate_result
  
      L=0.
      do i=1,n
         P(i)=i
      end do

      do i=1,n
         L(i,i)=1.d0 !we set the diagonal equal to 1
      end do 

      !write(*,*) "Initial matrix L is:"
      !call printmat(L,n)

      do i=1,n-1
         !write(*,*) "External loop:",i
         piv=A(i,i)
         indP=0
         !scan along rows for partial pivoting
         do j=i,n 
            !write(*,*) "Row pivoting scan:", j
            if(abs(A(j,i))>piv) then
              piv=A(j,i)
              indP=j
              !write(*,*) "Pivot index is:", indP
            end if
         end do
         if(indP >0) then
           !we swap the lines of the permutation matrix (well, not exactly, but we keep track of the exchanges of the rows)
           tmpindx=P(i)
           P(i)=P(indP)
           P(indP)=tmpindx
           !we swap the rows of the A (U) matrix
           tmprow(i:n)=A(i,i:n) !we don't need to copy all the row, the first part is 0
           A(i,i:n)=A(indP,i:n)
           A(indP,i:n)=tmprow(i:n)
           !we swap the rows of the L matrix
           tmprow(1:i-1)=L(i,1:i-1)
           L(i,1:i-1)=L(indP,1:i-1)
           L(indP,1:i-1)=tmprow(1:i-1)
         end if
         do j=i+1,n
            m=A(j,i)/piv
            L(j,i)=m !A becomes the lower matrix
            do k=i,n !now we reduce the remaining terms of the original matrix
               !write(*,*) "Subtracting m*U(",i,k,")=",m*U(i,k)," from U(",j,k,")=",U(j,k)
               A(j,k)=A(j,k)-m*A(i,k)
               if(abs(A(j,k)) .le. smallvalue) A(j,k)=0.
            end do
         end do
      end do       

      end subroutine ludecompose

      subroutine printmat(inmat,n)
      implicit none

      real(8), allocatable, dimension(:,:),intent(in) :: inmat
      integer, intent(in) :: n
      integer :: i

      do i=1,n
         write(*,*) inmat(i,1:n)
      end do
      write(*,*)

      end subroutine printmat 
                              
      subroutine printrealarr(inarr,n)
      implicit none

      real(8), allocatable, dimension(:),intent(in) :: inarr
      integer, intent(in) :: n
      integer :: i

      write(*,*) inarr(1:n)
      write(*,*)

      end subroutine printrealarr 

      subroutine printintarr(inarr,n)
      implicit none

      integer, allocatable, dimension(:),intent(in) :: inarr
      integer, intent(in) :: n
      integer :: i

      write(*,*) inarr(1:n)
      write(*,*)

      end subroutine printintarr 

      subroutine solve(L,U,P,X,B,n)
      implicit none
      real(8), dimension(1:n,1:n),intent(in) :: L,U
      real(8), dimension(1:n), intent(in) :: B
      real(8), dimension(1:n) :: X, Z
      integer, dimension(1:n),intent(in) :: P
      
      real(8) :: rowsum 
      integer, intent(in) :: n
      integer :: j, i, k, allocate_result
       
      !solve LZ=PB with forward substitution
      Z(1)=B(P(1))
      do i=2,n
         rowsum=B(P(i))
         do k=1,i-1
            rowsum=rowsum-L(i,k)*Z(k)
         end do
         Z(i)=rowsum
      end do
       
      !solve UX=Z  
      X(n)=Z(n)/U(n,n)
      do j=1,n-1
         i=n-j !we use this trick to count backward from n to 1
         rowsum=Z(i)
         do k=i+1,n
            rowsum=rowsum-U(i,k)*X(k)
         end do
         X(i)=rowsum/U(i,i)
      end do

      end subroutine solve
      
      subroutine getfun(Xin,F,fjac3,Fnorm,n,switch)
      implicit none
      integer, intent(in) :: n
      real(8), dimension(1:n) :: Xin, F
      real(8), dimension(1:n,1:n) :: fjac3
      real(8) :: xx,yy,pp
      real(8), intent(out) :: Fnorm
      integer :: switch !if 0 it does not compute the Jacobian
      real(8) :: B2,SB2,S2,U,D,en,rh,dendrh,dendpr
      integer :: errcode
            
      xx=Xin(1)
      yy=Xin(2)
      pp=Xin(3)

      B2=solver_parameters(1)
      SB2=solver_parameters(2)
      S2=solver_parameters(3)
      U=solver_parameters(4)
      D=solver_parameters(5)

      F(1)=(yy+B2)**2*xx-SB2/(yy**2)*(2*yy+B2)-S2
      F(2)=yy-pp+0.5*(1.+xx)*B2-0.5*SB2/(yy**2)-U
      rh=D*sqrt(1-xx)
      call eos_energy_with_ders(rh,en,pp,dendrh,dendpr,errcode)
      F(3)=yy-(pp+en)/(1-xx)
      if(switch .eq. 0) return
      fjac3(1,1)=(yy+B2)**2
      fjac3(1,2)=2*(yy+B2)*xx+(2*SB2*(2*yy+B2))/(yy**3)-2*SB2/(yy**2)
      fjac3(1,3)=dendrh*(-D/(2.*sqrt(1.-xx)))
      fjac3(2,1)=B2/2.
      fjac3(2,2)=1.+SB2/(yy**3)
      fjac3(2,3)=-1.
      fjac3(3,1)=-(pp+en)/((1-xx)**2)
      fjac3(3,2)=1.
      fjac3(3,3)=-(1+dendpr)/(1-xx)

      end subroutine getfun
                              

! *****************************************************************************

subroutine solver7(x,f,df)
  implicit none
  real(8),intent(in ) :: x
  real(8),intent(out) :: f,df

  real(8),parameter :: third=1.d0/3.d0

  integer :: iter
  real(8)    :: v2,rho,c0,c2,c3,dw,dlogw,dc2,dc3,wb,vb2

  real(8) gam,w, B2,SB2,S2,D,UU

  gam=(gamma_adindex-1)/gamma_adindex

  B2=solver_parameters(1)
  SB2=solver_parameters(2)
  S2=solver_parameters(3)
  UU=solver_parameters(4)
  D=solver_parameters(5)
  
  v2=x
  rho=D*sqrt(1.d0-v2)

  c3=1.d0-gam*(1.d0-v2)
  c2=gam*rho+.5d0*B2*(1.d0+v2)-UU
  c0=-0.5d0*SB2

! For every x=v2, we solve a cubic in W of the form: 
! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)

  w=max(-c2/c3,(-c0/c3)**third)
  do iter=1,50
    dw=-((c3*w+c2)*w**2+c0)/((3*c3*w+2*c2)*w)
    if (abs(dw/w)<1.e-10) exit
    w=w+dw
  end do

  dc3=gam
  dc2=.5*(B2-gam*rho/(1.-v2))
  dlogw=-(dc3*w+dc2)/(3*c3*w+2*c2)
  wb=w+B2
  vb2=SB2/w**2
  f=wb**2*v2-(2*w+B2)*vb2-S2
  df=wb*(wb+2*dlogw*(w*v2+vb2))
  w_algo_solver7=w

end subroutine solver7

! *****************************************************************************

subroutine solver8(x,f,df) !echo 2007 solver, adapted for p=e/3
  implicit none
  real(8),intent(in ) :: x
  real(8),intent(out) :: f,df

  real(8),parameter :: third=1.d0/3.d0

  integer :: iter
  real(8)    :: v2,rho,c0,c2,c3,dw,dlogw,dc2,dc3,wb,vb2

  real(8) w, B2,SB2,S2,D,UU
  real(8), parameter :: gam=0.25d0


  B2=solver_parameters(1)
  SB2=solver_parameters(2)
  S2=solver_parameters(3)
  UU=solver_parameters(4)
  D=solver_parameters(5)
  
  v2=x

  c3=1.d0-gam*(1.d0-v2)
  c2=.5d0*B2*(1.d0+v2)-UU !the last term is essentially U
  c0=-0.5d0*SB2

! For every x=v2, we solve a cubic in W of the form: 
! c3*W^3+c2*W^2+c0=0 (c3>0, c0<=0)

  if (c0==0.d0) then
    w=-c2/c3
  else
    w=max(-c2/c3,(-c0/c3)**third)
    do iter=1,100
      dw=-((c3*w+c2)*w**2+c0)/((3*c3*w+2*c2)*w)
      if (abs(dw/w)<1.d-12) exit
      w=w+dw
    end do
  end if

  dc3=gam
  dc2=.5d0*B2
  dlogw=-(dc3*w+dc2)/(3*c3*w+2*c2)
  wb=w+B2
  vb2=SB2/w**2
  f=wb**2*v2-(2*w+B2)*vb2-S2
  df=wb*(wb+2*dlogw*(w*v2+vb2))
  w_algo_solver8_e_3p=w

end subroutine solver8

! *****************************************************************************

subroutine solver10_approximate(b,c,d,oldsol,x) !approximate solver
  implicit none
  real(8),intent(in ) :: b,c,d,oldsol

  integer, parameter :: maxiterations=100
  integer :: iter
  real(8)    :: x,dx

  x=oldsol
  do iter=1,maxiterations
    dx=-(12*x**3+3*b*x**2+2*c*x+d)
    if (abs(dx/x)<1.d-12) exit
    x=x+dx
  end do
  if(iter==maxiterations) then
    write(*,*) "Sorry, the solver10_approximate subroutine failed... Leaving."
    call exit(12)
  end if
end subroutine solver10_approximate

! *****************************************************************************

function rtsafe(func,x0,x1,x2,xacc,ierr)

  interface
    subroutine func(x,f,df)
      real(8),intent(in ) :: x
      real(8),intent(out) :: f,df
    end subroutine func
  end interface
  real(8)                :: rtsafe
  real(8)   ,intent(in)  :: x0,x1,x2,xacc
  integer,intent(out) :: ierr

  integer,parameter  :: maxit=100

  integer :: j
  real(8)    :: f,fl,fh,df,xl,xh,swap,dx,dxold,temp
 
  ierr=0
  xl=x1; xh=x2
  rtsafe=x0 !0.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call func(rtsafe,f,df)
  do j=1,maxit
    if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f)>0. .or. &
      abs(2.*f)>abs(dxold*df)) then
      dxold=dx
      dx=0.5*(xh-xl)
      rtsafe=xl+dx
      if (xl==rtsafe) then
         !we call rtsafe again because the value of rtsafe=w has been modified
         call func(rtsafe,f,df)
         return
      end if
    else
      dxold=dx
      dx=f/df
      temp=rtsafe
      rtsafe=rtsafe-dx
      if (temp==rtsafe) then
         !we call rtsafe again because the value of rtsafe=w has been modified
         call func(rtsafe,f,df)
         return
      end if
    end if
    if (abs(dx)<xacc) return
    call func(rtsafe,f,df)
    if (f<0.) then
      xl=rtsafe
      fl=f
    else
      xh=rtsafe
      fh=f
    end if
  end do
! stop 'RTSAFE exceeding maximum iterations'
  ierr=1
  return

end function rtsafe

! *****************************************************************************

real(8) function get_x9(y,s2,b2,sb2)
implicit none

real(8) y,s2,b2,sb2

get_x9=(s2*y**2+b2*sb2+2*sb2*y)/(y**2*(b2**2+2*b2*y+y**2))

end function get_x9

! *****************************************************************************
end module system_eqgp

! *****************************************************************************
