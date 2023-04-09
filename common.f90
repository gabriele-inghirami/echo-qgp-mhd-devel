! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *         
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2016,2018 The ECHO-QGP team                           * 
! *                                                                           *
! *  File: common.f90                                                         *
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
! *           Gabriele Inghirami (inghirami@fias.uni-frankfurt.de)            *
! *                                                                           *
! *  Contributors:                                                            *
! *                                                                           *
! *  Acknowledgments:                                                         *
! *                                                                           * 
! *****************************************************************************

module common
!-- Main settings to be accessed throughout the code
  use holib
  use parallel
  use eos,only: eos_tab_name, eqstate, przero, enezero,numerically
  implicit none
  integer :: init_type

  real(8) :: t
  real(8) :: tout, thyp
  real(8) :: tmin, tmax
  real(8) :: t_collect_d_a
  integer :: errcode
 
  real(8), parameter :: hbar=0.197326 !actually, it is hbar*c
  real(8), parameter :: sqrh=0.444214 !it is sqrt(hbar*c)=sqrt(0.197326)
  real(8), parameter :: pi=acos(-1.)
  real(8), parameter :: am=0.93891897   ! nucleon mass
  real(8), parameter :: Qbig=0.30282212 !elementary electric charge, it is sqrt(4*Pi*alpha)
  real(8) :: maxspeed !value to reset velocity if it becomes superluminal
 
  real(8) :: B_amplify_factor !the initial magnetic field is multiplied by this factor before introducing any cuts


  integer :: nv,nu,nov,nou

  integer,parameter :: krh=1,kvx=2,kvy=3,kvz=4,kpr=5,kpibu=6,kpixy=7, kpixz=8,kpiyz=9,kpixx=10,kpiyy=11,kpizz=12

  integer :: kbx, kby, kbz, kex, key, kez !their values will be assigned later depending on the type of run (viscid/unviscid)
  integer :: krc !charge density in the comoving frame
  integer :: kglm !index for the General Lagrange Multiplier, i.e. for divergence cleaning
  integer :: kvold_end !index at which vold should end

  integer, parameter :: dtx=1, dty=2, dtz=3, dxx=4, dxy=5, dxz=6, dyx=7, dyy=8, dyz=9, dzx=10, dzy=11, dzz=12
  integer, parameter :: dtt=13, dxt=14, dyt=15, dzt=16
  integer, parameter :: dtet=17, dtex=18, dtey=19, dtez=20

  real(8), allocatable, dimension(:,:,:,:) :: derivatives_all

  integer,parameter :: nderivatives=20

  real(8) :: pitt, pitx, pity, pitz ! components of shear viscous tensor that are computed from the indipendent ones
  real(8) :: pizz
  character(len=2) :: obtained !it obtains another shear viscous tensor component from the others using null trace condition

  integer :: coordinates !flag that select the coordinates system: 1=Minkowski, 2=Bjorken
  integer, parameter :: MINKOWSKI=1, BJORKEN=2 

! if .true. then we'll evolve also electric and magnetic fields
  logical :: mhd, rmhd
  logical :: print_rho_comov !decides whether to print or not the charge density in the comoving frame on the freezeout hypersurface
  real(8) :: time_comp_rho_comov !time at which the electric charge density in the comoving frame is computed

! interval between output printing, interval between log printing, maximum timestep
  real(8) :: dtout, dtlog, maxdt

! temperature at which simulations ends
  real(8) :: temp_end

! restart type: 0=no restart possibilities, 1=restart possible only from last

  integer :: restart_type

  logical :: restarted_run=.false.

  logical :: viscosity,bulkvis
 
  real(8) :: eta_over_s

! parameters to tune the relaxation times
  real(8) :: tau_pi_coeff
  real(8) :: zeta_vis_coeff

  real(8) :: g_cov(3),gp,gm, g_cov0, g_cov_3_old

!-- Grid points (x main direction, if parallel nx/npe points for each pe)

!
  real(8),dimension(3)   :: d_ini
  real(8),dimension(2,3) :: xlim
!
!-- Grid extension
  real(8) :: x1,x2,y1,y2,z1,z2

!-- Selection of printed variables into the output files
  type out_type
       integer density
       integer vx
       integer vy
       integer vz
       integer pressure
       integer energy_density
       integer temperature
       integer entropy_density
       integer bulk
       integer pitt
       integer pitx
       integer pity
       integer pitz
       integer pixy
       integer pixz
       integer piyz
       integer pixx
       integer piyy
       integer pizz
       integer v0
       integer bx
       integer by
       integer bz
       integer ex
       integer ey
       integer ez
       integer glm
       integer rc
  end type out_type

  character(len=9)  :: dir_out 
  character(len=:), allocatable :: prefix_dir
  character(len=120) :: param_file_cmdline
  character(len=:), allocatable :: param_file
  integer :: commandline_outdir_length


  integer,allocatable,dimension(:) :: recvcounts,displs

  type(out_type) out_sel

  integer :: output_precision

!-- Possible initialization types 
  integer, parameter :: GLAUBER_GEOMETRIC=0, SHOCK_TUBE_2D=1, SHEAR_VISCOUS_1D=2, GLAUBER_MONTECARLO=3, GUBSER_IDEAL=4,&
  &GUBSER_VISCOUS=5, TABULATED_INIT=6, ALFVEN2D=7, BLAST3D=8, MHD1DBJ=9, SODMHD=10, LYU=11, BERHAD=12, ROTOR=13, OT=14,&
  &DBG=15, EXT_CONS=16, EXT_GLISSANDO=17
       
!-- Numerical settings (used in evolve.f90 and work.f90)

  real(8) :: cfl

  logical,parameter :: ho=.false.

  character(3),parameter :: solver='HLL'

  !reconstruction algorithm
  character(5) :: recal

  character(10),parameter :: der='    DER-E6'

  integer,parameter :: ngc=3

  integer, parameter :: RK2=2,RK3=3,SSP=4

  character(4) :: nrk_string

  integer :: nrk

!-- Declarations

  integer,allocatable,dimension(:) :: kv,ku

  integer,allocatable,dimension(:,:,:) :: ibc
  integer :: boundary_cond

  real(8),allocatable,dimension(:) :: x,ddx,x0
  real(8),allocatable,dimension(:) :: y,ddy,y0
  real(8),allocatable,dimension(:) :: z,ddz,z0

  real(8),dimension(:,:,:,:),allocatable :: v,u
  real(8),dimension(:,:,:,:),allocatable :: vold, uold, vnewstep, unewstep
  real(8),dimension(:,:,:,:),allocatable :: deriv 

  real(8), dimension(:,:,:,:),allocatable :: vall, uall, vallold, uallold, vderivatives 
  !DDD they will be allocated only by first processor and it will store global grid data

  real(8) :: timeold !it holds the value of time before last timestep
  real(8) :: timeinterval !it holds the value of time interval for time derivatives

  real(8),dimension(:,:,:,:),allocatable :: w,u0,ax,ay,az,wsend,wrecv
  real(8),dimension(:,:,:,:,:),allocatable :: du,du_stiff !it is used to store the intermediate steps in IMEX-SSP schemes

  !parameters specific for Glauber-MonteCarlo simulations
  integer :: nconf         !number of nuclear configurations for generating events
  integer :: events_start  ! number of the event from which to start
  integer :: events_stop   ! number of the event at which to stop
  integer :: nbcoll      ! number of impact parameters per config.
  real(8) :: kappa          ! model parameters (taken from Eskola et al., PRC83, 034901)
  real(8) :: sig_mc         ! smearing parameter 
  real(8) :: ah             !initial hardness or collision weight (0=pure participant dependence, 1=pure collision dependence)
  integer :: kind_of_collision !1=AA, 2=dA, 3=pA
  logical :: fixed_b        !fixed impact parameter in GMC initial conditions
  real    :: eprot          !energy of the proton beam in pA collisions
  integer :: min_participants !minimum number of_participants to register the event for the AA collisions
  integer :: avg_events  !instead of performing event by ev. simulations, it does just 1 run with averaged in. cond. (0=off,1=on)
  
  !id of the run
  integer :: id_of_the_run

  !output directory name that can be given from command line
  character(len=120) :: custom_output_directory
  integer :: len_of_cmd_line_out_dir

  logical :: derivatives_out !flag to enable/disable the print of the derivatives into the output
  logical :: flows_out     !flag to enable/disable the print of the flows (directed, elliptic and eccentricity) into the output

  logical :: enable_smooth_viscosity !flag to enable/disable the smoothing of viscous tensor components below a certain temper.
  real(8) :: smooth_temp    !temperature under which the values viscosity tensor components are reduced

  logical :: run_crashed !flag enabled when something goes wrong in multiple runs (eg Glauber-MC) to skip to the next simulation

  character(5) :: name_of_nucleus
  
  character(5)          :: nucleus
  real(8)                  :: tauzero
  real(8)                  :: deta,sigeta,yb,projmass,radius,delta,roze,rads,sigma_in,bimpact,ecenter,rhocenter
  integer               :: zelectrons

  integer               :: ienentr         ! initial conditions flag: 0 (energy), 1 (entropy)

  !freezout and hypersurface computation section
  real(8), dimension(3) :: d_val(3)
  logical :: out_freeze   !.true. means that we print freeze-out hypersurfaces, .false. means not
  integer :: freeze_type  !0 for freezout based on temperature, 1 based on energy density
  real(8)  :: freeze_value   !freeze out value in GeV
  real(8)  :: freeze_time_interval !time interval for hypersurface computation

  !initialization with external files
  character(len=120)  :: input_edf !input file for entropy density distribution initialization
  character(len=120)  :: ext_cons_file !external file for initialization with conservative variables
  character(len=120)  :: ext_glissando_file !external file for initialization with Glissando

  !transverse plane analysis section
  real(8), allocatable, dimension(:) :: xcm, ycm, eccentricity, elliptic_flow, directed_flow
  logical :: tp_anal !it enables the computationo of the fluid transverse flow

  !parameters used for producing "tilted" initial energy density profiles
  real(8) :: ueta_coeff !for initializations with u^eta!=0
  logical :: tilting=.true. !initial energy density tilting as in http://arxiv.org/pdf/1501.04468v2.pdf
  real(8) :: etam !eta_m to produce initial en. dens. tilting as in http://arxiv.org/pdf/1501.04468v2.pdf, disabled if <0

  !parameters related to MHD simulations
  real(8) :: sigma_el_cond !electrical conductivity of QGP
  real(8) :: sigma_chiral_cond !chiral conductivity of QGP
  real(8) :: ratio_chir_el !ratio between sigma_chiral_cond and sigma_el_cond
  integer :: magfield_type !type of magnetic field (1=classical+chiral,2=classical only,3=chiral only)
  character(len=120)     :: input_B_field !input file for magnetic field initialization
  logical :: divclean !it enables divergence cleaning
  real(8) :: glm_alpha !look at: Journal of Computational Physics 229 (2010) 2117–2138
  real(8) :: glm_ch, glm_ch_old=0.5 !look at: Journal of Computational Physics 229 (2010) 2117–2138
  real(8), dimension(3) :: aflux !array in which to store the maximum speed along each direction
  logical :: dump_initial_B !if true it will dumpen B when the energy density is below a certain treshold
  real(8) :: edens_for_B_dump !value for the energy density, in GeV/fm^3, under which to dump the B field
  real(8) :: B_th_press_ratio !ratio between thermal and magnetic pressure in the dumping zone
  real(8) :: pw !weight of participants in computing the initial B field with GMC initial conditions
  integer :: algo_solver !solver for the cons2prim subroutine in the MHD case: 1 3x3 system cernlib, 2 rootfunc, 3 3x3 Newton w bracketing,
                         ! 4 solver assuming ideal gas eos, 5 solver assuming e=3p, 9 exact solver

  !these variables are just used to be known at the system_cons2prim when dealing with IMEX-SSP schemes without passing them
  real(8) :: imex_alpha_coeff, dt_int
  
  !for debugging purposes: NOT IMPLEMENTED, the idea was to use this variables, global but different for each process
  !to keep track of the cell currently evolved
  !integer :: ielab, jelab, kelab

  !variables and parameters for initialization:
  real(8)                  :: b1
  real(8)                  :: beta,weightbin_inv,weightwound_inv,weightbin,weightwound
  real(8)                  :: amassa
  integer               :: nmax    !number of elements for the thickness vector (named thick)
  real(8)                  :: gamma_col, beta_col
  real(8)                  :: thminus, thplus
  logical               :: search_in_cond !flag used for initialization type 1 to proceed with variable initialization
  real(8),parameter        :: dz1=0.001      ! dz for thickness function integration (fm)
  real(8),parameter        :: dr=0.001       ! dr step for thickness function calculation (fm)
  !real(8),parameter        :: dz1=0.01        ! dz for thickness function integration (fm)
  !real(8),parameter        :: dr=0.01         ! dr step for thickness function calculation (fm)
  real(8),parameter        :: acc=1.e-07     ! accuracy of thickness function integration

  !to keep track of where the parameter has been set: 0=default(common.f90), 1=parameter file, 2=command line

  integer :: init_type_option, coordinates_option, viscosity_flag_option, mhd_flag_option, bulkvis_flag_option
  integer :: smooth_temp_option, nx_option, ny_option, nz_option
  integer :: x1_option, x2_option, y1_option, y2_option, z1_option, z2_option, tmin_option, tmax_option, temp_end_option
  integer :: dtlog_option, dtout_option, output_precision_option, maxdt_option, restart_type_option, cfl_option
  integer :: nucleus_option, rads_option, sigma_in_option, bimpact_option, ienentr_option, ah_option, ecenter_option
  integer :: enezero_option, rhocenter_option, deta_option, sigeta_option, eta_over_s_option, tau_pi_coeff_option
  integer :: obtained_option, eqstate_option, eos_tab_name_option, numerically_option, nconf_option, nbcoll_option
  integer :: events_start_option, events_stop_option, kappa_option, sig_mc_option, kind_of_collision_option
  integer :: out_freeze_flag_option, freeze_type_option, freeze_time_interval_option, input_edf_option, etam_option
  integer :: ueta_coeff_option, input_B_field_option, density_option, vx_option, vy_option, vz_option, pressure_option
  integer :: energy_density_option, temperature_option, entropy_density_option, bulk_option, pitt_option, pitx_option
  integer :: pity_option, pitz_option, pixy_option, pixz_option, piyz_option, pixx_option, piyy_option, pizz_option
  integer :: v0_option, bx_option, by_option, bz_option, ex_option, ey_option, ez_option, glm_option, freeze_value_option
  integer :: derivatives_flag_option, flows_flag_option, nuclei_data_file_option, custom_output_directory_option
  integer :: sigma_el_cond_option, B_amplify_factor_option, divclean_flag_option, glm_alpha_option, dump_initial_B_flag_option
  integer :: recal_option, edens_for_B_dump_option, B_th_press_ratio_option, pw_option, maxspeed_option, algo_solver_option
  integer :: ext_cons_file_option,boundary_cond_option, fixed_b_flag_option, eprot_option, min_participants_option
  integer :: sigma_chiral_cond_option, magfield_type_option, avg_events_option
  integer :: rc_option, print_rho_comov_flag, print_rho_comov_flag_option
  integer :: nrk_option,ext_glissando_file_option


  character(len=10) :: flagstring
  character(len=120) :: inputstring
  character(len=120) :: nuclei_data_file


contains

! *****************************************************************************

subroutine common_alloc
!-- Allocate main arrays

  integer :: i
  integer :: allocate_result


  allocate(kv(nv), ku(nu), STAT=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate kv or ku into common_alloc, contained into common.90"
    write(*,*)  "(source file common.f90)"
    call exit(1)
  end if
  kv=(/(i,i=1,nv)/)
  ku=(/(i,i=1,nu)/)

  allocate(v(ix1:ix2,iy1:iy2,iz1:iz2,1:nv),u(ix1:ix2,iy1:iy2,iz1:iz2,1:nu), stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate v or u multidim. array into common_alloc, contained into common.f90"
    write(*,*)  "(source file common.f90)"
    call exit(1)
  end if
  v=0.
  u=0.
  
  allocate(deriv(ix1:ix2,iy1:iy2,iz1:iz2,1:nderivatives), stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate deriv multidim. array into common_alloc, contained into common.f90"
    write(*,*)  "(source file common.f90)"
    call exit(1)
  end if
  deriv=0.
 
  if(mhd) then
    kvold_end=kbz 
    if(rmhd) then
      kvold_end=kbz 
    end if
  else
    kvold_end=kpr
  end if 
  allocate(vold(ix1:ix2,iy1:iy2,iz1:iz2,krh:kvold_end),uold(ix1:ix2,iy1:iy2,iz1:iz2,krh:kvold_end),&
  &vnewstep(ix1:ix2,iy1:iy2,iz1:iz2,krh:kvold_end),unewstep(ix1:ix2,iy1:iy2,iz1:iz2,krh:kvold_end),stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate vold, uold, vnewstep or unewstep into common_alloc"
    write(*,*)  "(source file common.f90)"
    call exit(1)
  end if
  vold=0.
  uold=0.
  vnewstep=0.
  unewstep=0.

  allocate(x(1:nx),ddx(1:nx),x0(0:nx), &
           y(1:ny),ddy(1:ny),y0(0:ny), &
           z(1:nz),ddz(1:nz),z0(0:nz), stat=allocate_result)
  if(allocate_result /=0) then 
    write(*,*)  "Error, I can't allocate x,dd,x0,y,ddy,y0,z,ddz or z0 into common_alloc"
    write(*,*)  "(source file common.f90)"
    call exit(1)
  end if
  x=0.; ddx=0.; x0=0.; y=0.; ddy=0.; y0=0.; z=0.; ddz=0.; z0=0.

  allocate(ibc(nv,2,3),stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*)  "Error, I can't allocate ibc into common_alloc, contained into common.f90"
    write(*,*)  "(source file common.f90)"
    call exit(1)
  end if

  if (ipe .eq. 0) then 
     allocate(vall(1:nx,1:ny,1:nz,1:nv), vderivatives(1:nx,1:ny,1:nz,4), stat=allocate_result)
     if(allocate_result /=0) then
       write(*,*)  "Proc.0 - Error, I can't allocate 'vall' or/and 'vderivatives' arrays into common_alloc"
       write(*,*)  "(source file common.f90)"
       call exit(1)
     end if
     vall=0.
     vderivatives=0.
  end if

  if ((ipe .eq. 0) .and. tp_anal) then !this can be improved if each processor computes its own grid section
                                       !and then we combine results using mpi functions
     allocate(xcm(1:nz),ycm(1:nz),eccentricity(1:nz),directed_flow(1:nz),elliptic_flow(1:nz), stat=allocate_result)
     if(allocate_result /=0) then
       write(*,*)  "Proc.0 - Error, I can't allocate arrays for transverse plane analysis into common_alloc"
       write(*,*)  "(source file common.f90)"
       call exit(1)
     end if
  end if

end subroutine common_alloc

! *****************************************************************************

subroutine common_grid(x1,x2,y1,y2,z1,z2)
!-- Define uniform grids

  real(8),intent(in) :: x1,x2,y1,y2,z1,z2

  integer :: i
  real(8)    :: dx,dy,dz

  dx=(x2-x1)/nx
  dy=(y2-y1)/ny
  dz=(z2-z1)/nz

  ddx(1:nx)=0.; if (nx>1) ddx(1:nx)=1./dx
  ddy(1:ny)=0.; if (ny>1) ddy(1:ny)=1./dy
  ddz(1:nz)=0.; if (nz>1) ddz(1:nz)=1./dz

  x=(/(x1+(i-0.5)*dx,i=1,nx)/); x0=(/(x1+i*dx,i=0,nx)/)
  y=(/(y1+(i-0.5)*dy,i=1,ny)/); y0=(/(y1+i*dy,i=0,ny)/)
  z=(/(z1+(i-0.5)*dz,i=1,nz)/); z0=(/(z1+i*dz,i=0,nz)/)
!
	d_ini(1)=dx
	d_ini(2)=dy
	d_ini(3)=dz
!
	xlim(1,1)=x1
	xlim(2,1)=x2
	xlim(1,2)=y1
	xlim(2,2)=y2
	xlim(1,3)=z1
	xlim(2,3)=z2
!
  d_val(1)=d_ini(1)
  d_val(2)=d_ini(2)
  if (nz>1) then 
    d_val(3)=d_ini(3)
  else 
    d_val(3)=1.0
  endif


end subroutine common_grid

! *****************************************************************************

subroutine common_get_all_parameters()
       implicit none
       integer bulkvis_flag,viscosity_flag,navierstokes_flag,derivatives_flag,flows_flag,out_freeze_flag,mhd_flag
       integer dump_initial_B_flag, divclean_flag, fixed_b_flag
       integer filerror, i
       integer :: read_status
       integer :: num_arguments !number of arguments from command line
       character(len=32) :: arg
       integer :: ia !just a counter
       integer :: skip=0 !a flag to control the behaviour of the parameters reading
       logical :: param_file_cmdline_set=.false.

       if(pe0) then
  
       
014     format(17x,f8.0,9x,f7.0,13x,f8.0,28x,f8.0,46x,i3)

075     format(a10,a120)

        !default values for parameters:
        init_type=0
        init_type_option=0

        coordinates=2
        coordinates_option=0

        viscosity_flag=0
        viscosity_flag_option=0

        bulkvis_flag=0
        bulkvis_flag_option=0
         
        mhd_flag=0
        mhd_flag_option=0

        divclean_flag=0
        divclean_flag_option=0

        glm_alpha=5.
        glm_alpha_option=0

        dump_initial_B_flag=1
        dump_initial_B_flag_option=0

        edens_for_B_dump=0.05
        edens_for_B_dump_option=0

        B_th_press_ratio=0.01
        B_th_press_ratio_option=0

        algo_solver=1
        algo_solver_option=0

        smooth_temp=0.08
        smooth_temp_option=0
        
        nx=120
        nx_option=0
 
        ny=120
        ny_option=0
       
        nz=120
        nz_option=0

        x1=-12.
        x1_option=0

        x2=12.
        x2_option=0

        y1=-12.
        y1_option=0
    
        y2=12.
        y2_option=0

        z1=-12.
        z1_option=0

        z2=12.
        z2_option=0

        tmin=1.
        tmin_option=0

        tmax=11.
        tmax_option=0

        temp_end=0.
        temp_end_option=0

        dtlog=0.01     
        dtlog_option=0

        dtout=1.
        dtout_option=0

        output_precision=8
        output_precision_option=0

        maxdt=0.002         
        maxdt_option=0

        restart_type=0
        restart_type_option=0

        cfl=0.4
        cfl_option=0

        recal="MPE5"
        recal_option=0

        nrk_string="RK2"
        nrk_option=0

        maxspeed=0.995
        maxspeed_option=0

        nucleus="Au"
        nucleus_option=0

        rads=200. 
        rads_option=0

        sigma_in=42.
        sigma_in_option=0

        bimpact=5.
        bimpact_option=0

        ienentr=0
        ienentr_option=0

        ah=0.15
        ah_option=0

        ecenter=40.
        ecenter_option=0
       
        enezero=0.001
        enezero_option=0

        rhocenter=1.
        rhocenter_option=0

        deta=2.0
        deta_option=0
  
        sigeta=1.
        sigeta_option=0

        eta_over_s=0.08
        eta_over_s_option=0

        tau_pi_coeff=5.
        tau_pi_coeff_option=0

        obtained="zz"
        obtained_option=0

        eqstate=1
        eqstate_option=0

        eos_tab_name="qcdIEOS0.dat"
        eos_tab_name_option=0

        numerically=0
        numerically_option=0

        nconf=100
        nconf_option=0

        nbcoll=1
        nbcoll_option=0

        fixed_b_flag=0
        fixed_b_flag_option=0

        events_start=1
        events_start_option=0

        events_stop=5
        events_stop_option=0

        kappa=19.0
        kappa_option=0

        sig_mc=0.8
        sig_mc_option=0

        eprot=200.
        eprot_option=0
 
        min_participants=18
        min_participants_option=0

        kind_of_collision=1
        kind_of_collision_option=0

        avg_events=0
        avg_events_option=0

        out_freeze_flag=0
        out_freeze_flag_option=0
 
        print_rho_comov_flag=0
        print_rho_comov_flag_option=0
 
        freeze_type=1
        freeze_type_option=0

        freeze_value=0.5
        freeze_value_option=0
 
        freeze_time_interval=0.1
        freeze_time_interval_option=0

        input_edf="ed.dat"
        input_edf_option=0

        ext_cons_file="cons.dat"
        ext_cons_file_option=0

        ext_glissando_file="glissando.dat"
        ext_glissando_file_option=0

        etam=-1
        etam_option=0

        ueta_coeff=0.
        ueta_coeff_option=0

        input_B_field="initialB.dat"
        input_B_field_option=0

        boundary_cond=2
        boundary_cond_option=0

        out_sel%density=1 
        density_option=0

        out_sel%vx=1
        vx_option=0

        out_sel%vy=1
        vy_option=0

        out_sel%vz=1 
        vz_option=0

        out_sel%pressure=1
        pressure_option=0

        out_sel%energy_density=1
        energy_density_option=0

        out_sel%temperature=0
        temperature_option=0

        out_sel%entropy_density=0
        entropy_density_option=0

        out_sel%bulk=0
        bulk_option=0

        out_sel%pitt=0
        pitt_option=0

        out_sel%pitx=0
        pitx_option=0

        out_sel%pity=0
        pity_option=0

        out_sel%pitz=0
        pitz_option=0

        out_sel%pixy=0
        pixy_option=0

        out_sel%pixz=0
        pixz_option=0

        out_sel%piyz=0
        piyz_option=0

        out_sel%pixx=0
        pixx_option=0

        out_sel%piyy=0
        piyy_option=0

        out_sel%pizz=0
        pizz_option=0

        out_sel%v0=0
        v0_option=0

        out_sel%bx=1
        bx_option=0

        out_sel%by=1
        by_option=0

        out_sel%bz=1
        bz_option=0

        out_sel%ex=1
        ex_option=0

        out_sel%ey=1
        ey_option=0

        out_sel%ez=1    
        ez_option=0

        out_sel%glm=0    
        glm_option=0

        out_sel%rc=0    
        rc_option=0

        derivatives_flag=0
        derivatives_flag_option=0

        flows_flag=0
        flows_flag_option=0

        nuclei_data_file="nuclear_data.dat"       
        nuclei_data_file_option=0

        custom_output_directory=""
        custom_output_directory_option=0

        sigma_el_cond=0.0058
        sigma_el_cond_option=0

        sigma_chiral_cond=0.0015
        sigma_chiral_cond_option=0

        magfield_type=1
        magfield_type_option=0

        B_amplify_factor=1.
        B_amplify_factor_option=0

        pw=0.5
        pw_option=0

        !now reading command line arguments
        num_arguments=command_argument_count()

        !we check if we should use a non standard parameter file
        if (num_arguments .gt. 0) then
         do ia = 1, num_arguments
          if(skip .eq. 1) then
           skip=0
           cycle
          end if
          call get_command_argument(ia, arg)

          if (arg .eq. "-PARAM_FILE") then
           call get_command_argument(ia+1, arg)
           read(arg,"(a120)") param_file_cmdline
           param_file_cmdline_set=.true.
           exit
          end if
         end do
        end if 

        if(param_file_cmdline_set) then
          allocate(character(len=len_trim(adjustl(param_file_cmdline))) :: param_file)
          param_file=trim(adjustl(param_file_cmdline))
        else
          allocate(character(len=9) :: param_file)
          param_file='param.dat'
        end if


        !first, we check that parameter file exists
        open(unit=29,status='OLD',file=param_file, iostat=filerror, form='formatted')
        if (filerror .ne. 0) then
            write(*,*) "File "//param_file//" cannot be opened, so I'm forced to quit!"
            close(29)
            call exit(1)
        end if

 
        do 
         read(29,075,IOSTAT=filerror) flagstring,inputstring
         if(filerror .lt. 0) then
          exit
         else
          if((flagstring(1:1) .eq. "!") .or. (flagstring .eq. "") .or. (flagstring(1:1) .eq. "#") .or.&
            & (flagstring(1:2) .eq. "//")) then
            cycle
          else if(flagstring .eq. "INIT_TYPE=") then
            read(inputstring,"(i2)") init_type
            init_type_option=1
          else if(flagstring .eq. "COORD....=") then
            read(inputstring,"(i1)") coordinates 
            coordinates_option=1
          else if(flagstring .eq. "VISCOUS..=") then
            read(inputstring,"(i1)") viscosity_flag
            viscosity_flag_option=1
          else if(flagstring .eq. "BULK.....=") then
            read(inputstring,"(i1)") bulkvis_flag 
            bulkvis_flag_option=1
          else if(flagstring .eq. "MHD......=") then
            read(inputstring,"(i1)") mhd_flag 
            mhd_flag_option=1
          else if(flagstring .eq. "DIVCLEAN.=") then
            read(inputstring,"(i1)") divclean_flag 
            divclean_flag_option=1
          else if(flagstring .eq. "GLM_PARAM=") then
            read(inputstring,"(f14.0)") glm_alpha
            glm_alpha_option=1
          else if(flagstring .eq. "DUMP_IN_B=") then
            read(inputstring,"(i1)") dump_initial_B_flag 
            dump_initial_B_flag_option=1
          else if(flagstring .eq. "B_DUM_EN.=") then
            read(inputstring,"(f14.0)") edens_for_B_dump
            edens_for_B_dump_option=1
          else if(flagstring .eq. "Bp_ov_Tp.=") then
            read(inputstring,"(f14.0)") B_th_press_ratio
            B_th_press_ratio_option=1
          else if(flagstring .eq. "CUT_TEMP.=") then
            read(inputstring,"(f14.0)") smooth_temp
            smooth_temp_option=1
          else if(flagstring .eq. "NX.......=") then
            read(inputstring,"(i6)") nx
            nx_option=1
          else if(flagstring .eq. "NY.......=") then
            read(inputstring,"(i6)") ny
            ny_option=1
          else if(flagstring .eq. "NZ.......=") then
            read(inputstring,"(i6)") nz
            nz_option=1
          else if(flagstring .eq. "XMIN.....=") then
            read(inputstring,"(f14.0)") x1
            x1_option=1
          else if(flagstring .eq. "XMAX.....=") then
            read(inputstring,"(f14.0)") x2
            x2_option=1
          else if(flagstring .eq. "YMIN.....=") then
            read(inputstring,"(f14.0)") y1
            y1_option=1
          else if(flagstring .eq. "YMAX.....=") then
            read(inputstring,"(f14.0)") y2
            y2_option=1
          else if(flagstring .eq. "ZMIN.....=") then
            read(inputstring,"(f14.0)") z1
            z1_option=1
          else if(flagstring .eq. "ZMAX.....=") then
            read(inputstring,"(f14.0)") z2
            z2_option=1
          else if(flagstring .eq. "TSTART...=") then
            read(inputstring,"(f14.0)") tmin
            tmin_option=1
          else if(flagstring .eq. "TSTOP....=") then
            read(inputstring,"(f14.0)") tmax
            tmax_option=1
          else if(flagstring .eq. "TEMP_END.=") then
            read(inputstring,"(f14.0)") temp_end
            temp_end_option=1
          else if(flagstring .eq. "DTLOG....=") then
            read(inputstring,"(f14.0)") dtlog
            dtlog_option=1
          else if(flagstring .eq. "DTOUT....=") then
            read(inputstring,"(f14.0)") dtout
            dtout_option=1
          else if(flagstring .eq. "OUTP_PREC=") then
            read(inputstring,"(i1)") output_precision
            output_precision_option=1
          else if(flagstring .eq. "MAXDT....=") then
            read(inputstring,"(f14.0)") maxdt
            maxdt_option=1
          else if(flagstring .eq. "RESTART..=") then
            read(inputstring,"(i1)") restart_type
            restart_type_option=1
          else if(flagstring .eq. "CFL......=") then
            read(inputstring,"(f14.0)") cfl
            cfl_option=1
          else if(flagstring .eq. "REC_ALGO.=") then
            read(inputstring,"(a5)") recal
            recal_option=1
          else if(flagstring .eq. "INT_ALGO.=") then
            read(inputstring,"(a3)") nrk_string
            nrk_option=1
          else if(flagstring .eq. "MAXSPEED.=") then
            read(inputstring,"(f14.0)") maxspeed
            maxspeed_option=1
          else if(flagstring .eq. "NUCLEUS..=") then
            read(inputstring,"(a5)") nucleus
            nucleus_option=1
          else if(flagstring .eq. "RADS.....=") then
            read(inputstring,"(f14.0)") rads
            rads_option=1
          else if(flagstring .eq. "SIGMA_IN.=") then
            read(inputstring,"(f14.0)") sigma_in
            sigma_in_option=1
          else if(flagstring .eq. "B........=") then
            read(inputstring,"(f14.0)") bimpact
            bimpact_option=1
          else if(flagstring .eq. "IENENTR..=") then
            read(inputstring,"(i1)") ienentr
            ienentr_option=1
          else if(flagstring .eq. "AH.......=") then
            read(inputstring,"(f14.0)") ah
            ah_option=1
          else if(flagstring .eq. "ECENTER..=") then
            read(inputstring,"(f14.0)") ecenter
            ecenter_option=1
          else if(flagstring .eq. "ENEZERO..=") then
            read(inputstring,"(f14.0)") enezero
            enezero_option=1
          else if(flagstring .eq. "RHOCENTER=") then
            read(inputstring,"(f14.0)") rhocenter
            rhocenter_option=1
          else if(flagstring .eq. "DETA.....=") then
            read(inputstring,"(f14.0)") deta
            deta_option=1
          else if(flagstring .eq. "SIGETA...=") then
            read(inputstring,"(f14.0)") sigeta
            sigeta_option=1
          else if(flagstring .eq. "ETA_S....=") then
            read(inputstring,"(f14.0)") eta_over_s
            eta_over_s_option=1
          else if(flagstring .eq. "TAU_PI_C.=") then
            read(inputstring,"(f14.0)") tau_pi_coeff
            tau_pi_coeff_option=1
          else if(flagstring .eq. "TRACE_IMP=") then
            read(inputstring,"(a2)") obtained
            obtained_option=1
          else if(flagstring .eq. "EOS......=") then
            read(inputstring,"(i1)") eqstate
            eqstate_option=1
          else if(flagstring .eq. "EOS_FILE.=") then
            read(inputstring,"(a120)") eos_tab_name
            eos_tab_name_option=1
          else if(flagstring .eq. "NUM_DER..=") then
            read(inputstring,"(i1)") numerically
            numerically_option=1
          else if(flagstring .eq. "NCONF....=") then
            read(inputstring,"(i6)") nconf
            nconf_option=1
          else if(flagstring .eq. "NBCOLL...=") then
            read(inputstring,"(i6)") nbcoll
            nbcoll_option=1
          else if(flagstring .eq. "FIXED_B..=") then
            read(inputstring,"(i1)") fixed_b_flag
            fixed_b_flag_option=1
          else if(flagstring .eq. "EV_START.=") then
            read(inputstring,"(i6)") events_start
            events_start_option=1
          else if(flagstring .eq. "EV_STOP..=") then
            read(inputstring,"(i6)") events_stop
            events_stop_option=1
          else if(flagstring .eq. "KAPPA....=") then
            read(inputstring,"(f14.0)") kappa
            kappa_option=1
          else if(flagstring .eq. "SIG......=") then
            read(inputstring,"(f14.0)") sig_mc
            sig_mc_option=1
          else if(flagstring .eq. "EPROT....=") then
            read(inputstring,"(f14.0)") eprot
            eprot_option=1
          else if(flagstring .eq. "AVG_IC...=") then
            read(inputstring,"(i1)") avg_events
            avg_events_option=1
          else if(flagstring .eq. "MIN_PART.=") then
            read(inputstring,"(i3)") min_participants
            min_participants_option=1
          else if(flagstring .eq. "COLLISION=") then
            read(inputstring,"(i1)") kind_of_collision
            kind_of_collision_option=1
          else if(flagstring .eq. "HYP_COMPU=") then
            read(inputstring,"(i1)") out_freeze_flag
            out_freeze_flag_option=1
          else if(flagstring .eq. "HYP_RHOEL=") then
            read(inputstring,"(i1)") print_rho_comov_flag
            print_rho_comov_flag_option=1
          else if(flagstring .eq. "FREEZKIND=") then
            read(inputstring,"(i1)") freeze_type
            freeze_type_option=1
          else if(flagstring .eq. "FREEZEVAL=") then
            read(inputstring,"(f14.0)") freeze_value
            freeze_value_option=1
          else if(flagstring .eq. "HYPSURFTI=") then
            read(inputstring,"(f14.0)") freeze_time_interval
            freeze_time_interval_option=1
          else if(flagstring .eq. "IN_D_FILE=") then
            read(inputstring,"(a120)") input_edf
            input_edf_option=1
          else if(flagstring .eq. "EX_C_FILE=") then
            read(inputstring,"(a120)") ext_cons_file
            ext_cons_file_option=1
          else if(flagstring .eq. "EX_G_FILE=") then
            read(inputstring,"(a120)") ext_glissando_file
            ext_glissando_file_option=1
          else if(flagstring .eq. "ETAM_TILT=") then
            read(inputstring,"(f14.0)") etam
            etam_option=1
          else if(flagstring .eq. "UETA_COEF=") then
            read(inputstring,"(f14.0)") ueta_coeff
            ueta_coeff_option=1
          else if(flagstring .eq. "B_in_FILE=") then
            read(inputstring,"(a120)") input_B_field
            input_B_field_option=1
          else if(flagstring .eq. "EL_COND..=") then
            read(inputstring,"(f14.0)")  sigma_el_cond
            sigma_el_cond_option=1
          else if(flagstring .eq. "CHIR_COND=") then
            read(inputstring,"(f14.0)")  sigma_chiral_cond
            sigma_chiral_cond_option=1
          else if(flagstring .eq. "MAGFIELDT=") then
            read(inputstring,"(i1)")  magfield_type
            magfield_type_option=1
          else if(flagstring .eq. "B_amp_fac=") then
            read(inputstring,"(f14.0)") B_amplify_factor
            B_amplify_factor_option=1
          else if(flagstring .eq. "bound_con=") then
            read(inputstring,"(i1)") boundary_cond
            boundary_cond_option=1
          else if(flagstring .eq. "density..=") then
            read(inputstring,"(i1)") out_sel%density
            density_option=1
          else if(flagstring .eq. "vx.......=") then
            read(inputstring,"(i1)") out_sel%vx
            vx_option=1
          else if(flagstring .eq. "vy.......=") then
            read(inputstring,"(i1)") out_sel%vy
            vy_option=1
          else if(flagstring .eq. "vz.......=") then
            read(inputstring,"(i1)") out_sel%vz
            vz_option=1
          else if(flagstring .eq. "pressure.=") then
            read(inputstring,"(i1)") out_sel%pressure
            pressure_option=1
          else if(flagstring .eq. "ene_dens.=") then
            read(inputstring,"(i1)") out_sel%energy_density
            energy_density_option=1
          else if(flagstring .eq. "temper...=") then
            read(inputstring,"(i1)") out_sel%temperature
            temperature_option=1
          else if(flagstring .eq. "entr_dens=") then
            read(inputstring,"(i1)") out_sel%entropy_density
            entropy_density_option=1
          else if(flagstring .eq. "bulk_visc=") then
            read(inputstring,"(i1)") out_sel%bulk
            bulk_option=1
          else if(flagstring .eq. "pi^tt....=") then
            read(inputstring,"(i1)") out_sel%pitt
            pitt_option=1
          else if(flagstring .eq. "pi^tx....=") then
            read(inputstring,"(i1)") out_sel%pitx
            pitx_option=1
          else if(flagstring .eq. "pi^ty....=") then
            read(inputstring,"(i1)") out_sel%pity
            pity_option=1
          else if(flagstring .eq. "pi^tz....=") then
            read(inputstring,"(i1)") out_sel%pitz
            pitz_option=1
          else if(flagstring .eq. "pi^xy....=") then
            read(inputstring,"(i1)") out_sel%pixy
            pixy_option=1
          else if(flagstring .eq. "pi^xz....=") then
            read(inputstring,"(i1)") out_sel%pixz
            pixz_option=1
          else if(flagstring .eq. "pi^yz....=") then
            read(inputstring,"(i1)") out_sel%piyz
            piyz_option=1
          else if(flagstring .eq. "pi^xx....=") then
            read(inputstring,"(i1)") out_sel%pixx
            pixx_option=1
          else if(flagstring .eq. "pi^yy....=") then
            read(inputstring,"(i1)") out_sel%piyy
            piyy_option=1
          else if(flagstring .eq. "pi^zz....=") then
            read(inputstring,"(i1)") out_sel%pizz
            pizz_option=1
          else if(flagstring .eq. "gamma....=") then
            read(inputstring,"(i1)") out_sel%v0
            v0_option=1
          else if(flagstring .eq. "bx.......=") then
            read(inputstring,"(i1)") out_sel%bx
            bx_option=1
          else if(flagstring .eq. "by.......=") then
            read(inputstring,"(i1)") out_sel%by
            by_option=1
          else if(flagstring .eq. "bz.......=") then
            read(inputstring,"(i1)") out_sel%bz
            bz_option=1
          else if(flagstring .eq. "ex.......=") then
            read(inputstring,"(i1)") out_sel%ex
            ex_option=1
          else if(flagstring .eq. "ey.......=") then
            read(inputstring,"(i1)") out_sel%ey
            ey_option=1
          else if(flagstring .eq. "ez.......=") then
            read(inputstring,"(i1)") out_sel%ez
            ez_option=1
          else if(flagstring .eq. "glm......=") then
            read(inputstring,"(i1)") out_sel%glm
            glm_option=1
          else if(flagstring .eq. "rc.......=") then
            read(inputstring,"(i1)") out_sel%rc
            rc_option=1
          else if(flagstring .eq. "derivativ=") then
            read(inputstring,"(i1)") derivatives_flag
            derivatives_flag_option=1
          else if(flagstring .eq. "flows....=") then
            read(inputstring,"(i1)") flows_flag
            flows_flag_option=1
          else if(flagstring .eq. "NUC_DATAF=") then
            read(inputstring,"(a120)") nuclei_data_file
            nuclei_data_file_option=1
          else if(flagstring .eq. "outdir...=") then
            read(inputstring,"(a120)") custom_output_directory 
            custom_output_directory_option=1
          else if(flagstring .eq. "PW.......=") then
            read(inputstring,"(f14.0)")  pw
            pw_option=1
          else if(flagstring .eq. "ALG_SOLV.=") then
            read(inputstring,"(i1)")  algo_solver
            algo_solver_option=1
          end if
         end if
        end do
        close(29)

        !now reading command line arguments
        num_arguments=command_argument_count()

        if (num_arguments .gt. 0) then
         do ia = 1, num_arguments
          if(skip .eq. 1) then
           skip=0
           cycle
          end if
          call get_command_argument(ia, arg)

          select case (arg)
           !we already processed this option, so we just read the value again without doing anything
           case("-PARAM_FILE") 
           call get_command_argument(ia+1, arg)
           read(arg,"(a120)") param_file_cmdline
           !no need of the option flag, because we can choose this option only from command line
           skip=1

           case("-INIT_TYPE") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i2)") init_type
            init_type_option=2
            skip=1
           case("-COORD") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") coordinates 
            coordinates_option=2
            skip=1
           case("-VISCOUS")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") viscosity_flag
            viscosity_flag_option=2
            skip=1
           case("-BULK")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") bulkvis_flag 
            bulkvis_flag_option=2
            skip=1
           case("-MHD")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") mhd_flag 
            mhd_flag_option=2
            skip=1
           case("-DIVCLEAN")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") divclean_flag 
            divclean_flag_option=2
            skip=1
           case("-GLM_PARAM")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") glm_alpha 
            glm_alpha_option=2
            skip=1
           case("-DUMP_IN_B")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") dump_initial_B_flag
            dump_initial_B_flag_option=2
            skip=1
           case("-B_DUM_EN")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") edens_for_B_dump
            edens_for_B_dump_option=2
            skip=1
           case("-Bp_ov_Tp")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") B_th_press_ratio
            B_th_press_ratio_option=2
            skip=1
           case("-CUT_TEMP")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") smooth_temp
            smooth_temp_option=2
            skip=1
           case("-NX")
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") nx
            nx_option=2
            skip=1
           case("-NY")
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") ny
            ny_option=2
            skip=1
           case("-NZ") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") nz
            nz_option=2
            skip=1
           case("-XMIN")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") x1
            x1_option=2
            skip=1
           case("-XMAX") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") x2
            x2_option=2
            skip=1
           case("-YMIN") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") y1
            y1_option=2
            skip=1
           case("-YMAX")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") y2
            y2_option=2
            skip=1
           case("-ZMIN")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") z1
            z1_option=2
            skip=1
           case("-ZMAX")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") z2
            z2_option=2
            skip=1
           case("-TSTART")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") tmin
            tmin_option=2
            skip=1
           case("-TSTOP")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") tmax
            tmax_option=2
            skip=1
           case("-TEMP_END")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") temp_end
            temp_end_option=2
            skip=1
           case("-DTLOG") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") dtlog
            dtlog_option=2
            skip=1
           case("-DTOUT") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") dtout
            dtout_option=2
            skip=1
           case("-OUTP_PREC")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") output_precision
            output_precision_option=2
            skip=1
           case("-MAXDT")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") maxdt
            maxdt_option=2
            skip=1
           case("-RESTART") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") restart_type
            restart_type_option=2
            skip=1
           case("-CFL") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") cfl
            cfl_option=2
            skip=1
           case("-REC_ALGO") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a5)") recal
            recal_option=2
            skip=1
           case("-INT_ALGO") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a3)") nrk_string 
            nrk_option=2
            skip=1
           case("-MAXSPEED") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") maxspeed
            maxspeed_option=2
            skip=1
           case("-NUCLEUS") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a5)") nucleus
            nucleus_option=2
            skip=1
           case("-RADS") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") rads
            rads_option=2
            skip=1
           case("-SIGMA_IN")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") sigma_in
            sigma_in_option=2
            skip=1
           case("-B") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") bimpact
            bimpact_option=2
            skip=1
           case("-IENENTR") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") ienentr
            ienentr_option=2
            skip=1
           case("-AH") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") ah
            ah_option=2
            skip=1
           case("-ECENTER")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") ecenter
            ecenter_option=2
            skip=1
           case("-ENEZERO")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") enezero
            enezero_option=2
            skip=1
           case("-RHOCENTER")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") rhocenter
            rhocenter_option=2
            skip=1
           case("-DETA") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") deta
            deta_option=2
            skip=1
           case("-SIGETA")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") sigeta
            sigeta_option=2
            skip=1
           case("-ETA_S")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") eta_over_s
            eta_over_s_option=2
            skip=1
           case("-TAU_PI_C")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") tau_pi_coeff
            tau_pi_coeff_option=2
            skip=1
           case("-TRACE_IMP")
            call get_command_argument(ia+1, arg)
            read(arg,"(a2)") obtained
            obtained_option=2
            skip=1
           case("-EOS") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") eqstate
            eqstate_option=2
            skip=1
           case("-EOS_FILE")
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") eos_tab_name
            eos_tab_name_option=2
            skip=1
           case("-NUM_DER")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") numerically
            numerically_option=2
            skip=1
           case("-NCONF")
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") nconf
            nconf_option=2
            skip=1
           case("-NBCOLL") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") nbcoll
            nbcoll_option=2
            skip=1
           case("-FIXED_B") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") fixed_b_flag
            fixed_b_flag_option=2
            skip=1
           case("-EV_START") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") events_start
            events_start_option=2
            skip=1
           case("-EV_STOP")
            call get_command_argument(ia+1, arg)
            read(arg,"(i6)") events_stop
            events_stop_option=2
            skip=1
           case("-KAPPA")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") kappa
            kappa_option=2
            skip=1
           case("-SIG") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") sig_mc
            sig_mc_option=2
            skip=1
           case("-EPROT") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") eprot
            eprot_option=2
            skip=1
           case("-AVG_IC") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") avg_events
            avg_events_option=2
            skip=1
           case("-MIN_PART") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i3)") min_participants
            min_participants_option=2
            skip=1
           case("-COLLISION") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") kind_of_collision
            kind_of_collision_option=2
            skip=1
           case("-HYP_COMPU") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_freeze_flag
            out_freeze_flag_option=2
            skip=1
           case("-HYP_RHOEL") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") print_rho_comov_flag
            print_rho_comov_flag_option=2
            skip=1
           case("-FREEZKIND") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") freeze_type
            freeze_type_option=2
            skip=1
           case("-FREEZEVAL") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") freeze_value
            freeze_value_option=2
            skip=1
           case("-HYPSURFTI") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") freeze_time_interval
            freeze_time_interval_option=2
            skip=1
           case("-IN_D_FILE") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") input_edf
            input_edf_option=2
            skip=1
           case("-EX_C_FILE") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") ext_cons_file
            ext_cons_file_option=2
            skip=1
           case("-EX_G_FILE") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") ext_glissando_file
            ext_glissando_file_option=2
            skip=1
           case("-ETAM_TILT")
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") etam
            etam_option=2
            skip=1
           case("-UETA_COEF") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") ueta_coeff
            ueta_coeff_option=2
            skip=1
           case("-B_in_FILE") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") input_B_field
            input_B_field_option=2
            skip=1
           case ('-EL_COND')
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") sigma_el_cond
            sigma_el_cond_option=2
            skip=1
           case ('-CHIR_COND')
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") sigma_chiral_cond
            sigma_chiral_cond_option=2
            skip=1
           case ('-MAGFIELDT') 
            read(arg,"(i1)") magfield_type
            magfield_type_option=2
           case("-B_amp_fac") 
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") B_amplify_factor
            B_amplify_factor_option=2
            skip=1
           case("-bound_con") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") boundary_cond
            boundary_cond_option=2
            skip=1
           case("-density") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%density
            density_option=2
            skip=1
           case("-vx")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%vx
            vx_option=2
            skip=1
           case("-vy")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%vy
            vy_option=2
            skip=1
           case("-vz")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%vz
            vz_option=2
            skip=1
           case("-pressure")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pressure
            pressure_option=2
            skip=1
           case("-ene_dens")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%energy_density
            energy_density_option=2
            skip=1
           case("-temper")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%temperature
            temperature_option=2
            skip=1
           case("-entr_dens")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%entropy_density
            entropy_density_option=2
            skip=1
           case("-bulk_visc")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%bulk
            bulk_option=2
            skip=1
           case("-pi^tt")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pitt
            pitt_option=2
            skip=1
           case("-pi^tx")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pitx
            pitx_option=2
            skip=1
           case("-pi^ty") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pity
            pity_option=2
            skip=1
           case("-pi^tz") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pitz
            pitz_option=2
            skip=1
           case("-pi^xy") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pixy
            pixy_option=2
            skip=1
           case("-pi^xz") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pixz
            pixz_option=2
            skip=1
           case("-pi^yz") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%piyz
            piyz_option=2
            skip=1
           case("-pi^xx") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pixx
            pixx_option=2
            skip=1
           case("-pi^yy") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%piyy
            piyy_option=2
            skip=1
           case("-pi^zz") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%pizz
            pizz_option=2
            skip=1
           case("-gamma") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%v0
            v0_option=2
            skip=1
           case("-bx") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%bx
            bx_option=2
            skip=1
           case("-by")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%by
            by_option=2
            skip=1
           case("-bz") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%bz
            bz_option=2
            skip=1
           case("-ex") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%ex
            ex_option=2
            skip=1
           case("-ey")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%ey
            ey_option=2
            skip=1
           case("-ez")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%ez
            ez_option=2
            skip=1
           case("-glm")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%glm
            glm_option=2
            skip=1
           case("-rc")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") out_sel%rc
            rc_option=2
            skip=1
           case("-derivativ")
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") derivatives_flag
            derivatives_flag_option=2
            skip=1
           case("-flows") 
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") flows_flag
            flows_flag_option=2
            skip=1
           case("-NUC_DATAF") 
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") nuclei_data_file
            nuclei_data_file_option=2
            skip=1
           case ('-outdir')
            call get_command_argument(ia+1, arg)
            read(arg,"(a120)") custom_output_directory
            custom_output_directory_option=2
            skip=1
           case ('-PW')
            call get_command_argument(ia+1, arg)
            read(arg,"(f14.0)") pw
            pw_option=2
            skip=1
           case ('-ALG_SOLV')
            call get_command_argument(ia+1, arg)
            read(arg,"(i1)") algo_solver
            algo_solver_option=2
            skip=1

           case default
            if(pe0) write(*,*) 'Fatal: unrecognized command-line option: ', arg
             stop
          end select
         end do
        end if

        !now we make a few checks about the consistency ot the chosen parameters
         tauzero=tmin !DDD AAA left for compatibily reasons, but tauzero should be changed with tmin everywhere inside the code

         if(.not.((output_precision .eq. 4) .or. (output_precision .eq. 8))) then
            write(*,*) 'Sorry, cannot understand the desired precision in the output files...'
            write(*,*) 'Please, check the OUTP_PREC parameters in the parameters file (it can be only 4 or 8).'
            call exit(1)
         end if

         if((init_type .eq. GLAUBER_GEOMETRIC) .and. (coordinates .eq. MINKOWSKI)) then
            write(*,*) "Since you chose to initialize the code with Geometric-Glauber initial conditions,"
            write(*,*) "we force ECHO-QGP to use Bjorken coordinates even if you set up Minkowski coordinates"
            coordinates=BJORKEN
         end if

         if(init_type .eq. OT) then
          if(boundary_cond .ne. 0) then !this test requires periodic boundary conditions
             boundary_cond=0
             boundary_cond_option=3
          end if

          if(mhd_flag .eq. 0) then
            write(*,*) "Since you chose to initialize the code to perform the Orszang-Tang test,"
            write(*,*) "we force ECHO-QGP to run MHD and not simply HD"
            mhd_flag=1
            mhd_flag_option=3
          end if

          if(eqstate .ne. 2) then
            write(*,*) "Since you chose to initialize the code to perform the Orszang-Tang test,"
            write(*,*) "we force ECHO-QGP to use the ideal gas EOS"
            eqstate=1
            eqstate_option=3
            algo_solver=4
            algo_solver_option=3
          end if
         end if

         if((ienentr .ne. 0) .and. (ienentr .ne. 1)) then
           write(*,*) 'Please, check the IENENTR parameter into the parameters file... Exiting now.'
           call exit(1)
         end if

         if(bulkvis_flag .eq. 0) then
            bulkvis=.false.
            out_sel%bulk=0
         else
            bulkvis=.true.
         end if

         if(viscosity_flag .eq. 0) then
            viscosity=.false.
            bulkvis=.false.
         else
            viscosity=.true.
         end if

         if(smooth_temp .gt. 0.) then
            enable_smooth_viscosity=.true.
         else
            enable_smooth_viscosity=.false.
         end if
  
         if(out_freeze_flag .eq. 0) then
            out_freeze=.false.
         else
            out_freeze=.true.
         end if

         if(mhd_flag .eq. 0.) then
            mhd=.false.
            rmhd=.false.
            if(out_sel%glm .eq. 1) then
              out_sel%glm=0
              glm_option=3
            end if  
         else if((mhd_flag .eq. 1) .or. (mhd_flag .eq. 2)) then
            mhd=.true.
            if(mhd_flag .eq. 2) then
              rmhd=.true.
            else
              rmhd=.false.
            end if
            if(solver .ne. 'HLL') then
              write(*,*) 'Sorry, but for MHD runs, only the HLL solver is implemented, so, please, change your parameters file.'
              call exit()
            end if
            if(divclean_flag .eq. 0) then
              divclean=.false.
              if(out_sel%glm .eq. 1) then
                out_sel%glm=0
                glm_option=3
              end if  
            else if(divclean_flag .eq. 1) then
              divclean=.true.
            else
              write(*,*) "Divergence cleaning option unknown: ", divclean_flag
              write(*,*) "Aborting."
              call exit(1)
            end if
            if(dump_initial_B_flag .eq. 1) then
              dump_initial_B=.true.
            else
              dump_initial_B=.false.
            end if
         else
              write(*,*) "MHD option can be only 0 (no MHD), 1 (ideal MHD) or 2 (resistive MHD)"
              write(*,*) "your choice is unclear:", mhd_flag
              write(*,*) "Aborting."
              call exit(1)
         end if

         if(mhd .and. viscosity) then
           write(*,*) "Sorry, but currently ECHO-QGP can deal with viscosity OR magnetic field, it cannot do resistive MHD"
           write(*,*) "Please, change your parameters file. I exit."
           call exit()
         end if

         if(print_rho_comov_flag .eq. 1) then
            if(out_freeze .and. mhd) then
              print_rho_comov=.true.
            else
              print_rho_comov=.false.
              print_rho_comov_flag_option=3
            end if
         end if
 
         !TESTMHDPRELIMINARE 
         !if((init_type .eq. GLAUBER_GEOMETRIC) .and. (coordinates .eq. MINKOWSKI) .and. (.not. mhd)) then
         !  write(*,*) "Since you chose to initialize the code with Geometric-Glauber initial conditions,"
         !  write(*,*) "we force ECHO-QGP to use Bjorken coordinates even if you set up Minkowski coordinates"
         !  coordinates=BJORKEN
         !end if
  
         !check the correctness of the choice of some parameters in the Glauber-MonteCarlo case
         if(init_type .eq. GLAUBER_MONTECARLO) then

           if((avg_events .ne. 0) .and. (kind_of_collision .ne. 1)) then
             write(*,*) "Sorry, but currently it is possible to run simulations with averaged events initial conditions"
             write(*,*) "only for nucleus-nucleus collisions. I quit."
             call exit(1)
           end if

           if(fixed_b_flag .eq. 1) then
             fixed_b=.true.
             nbcoll=1 !since we are using a fixed impact parameter, we have just 1 value for it
             nbcoll_option=3
           else if(fixed_b_flag .eq. 0) then
             fixed_b=.false.
           else
              write(*,*) "Unclear choice between fixed or randomly sampled b impact parameter. Quitting."
              call exit(1)
           end if
           if(events_stop .lt. events_start) then
             write(*,*) 'The EV_STOP parameter should be greater than the EV_START parameter'
             call exit(1)
           end if
           if(nconf .lt. events_stop) then
             write(*,*) 'The NCONF parameters should be greater than the EV_STOP parameter'
             call exit(1)
           end if
         else
           events_stop=events_start         
         end if
      
         if (ah < 0.0 .or. ah > 1.0) then
           if(pe0) then
             print *, "ah is the parameter for the hardness fraction, must be in [0,1]"
             print *, "your choice was:", ah
             print *, "check your parameter file again, please."
             print *, "I am going to quit."
           endif
           call exit(1)
        end if
        if (bimpact < 0.0) then
           if(pe0) then
             print *, "bimpact is the impact parameter in fm"
             print *, "your choice was:", bimpact
             print *, "check your parameter file again, please."
             print *, "I am going to quit."
           end if
           call exit(1)
        endif


        if(obtained .ne. 'no') then
          if(obtained .ne. 'zz') then
            if(pe0) then
              print *, "Please, check the value of the option TRACE_IMP into the parameters file"
              print *, "TRACE_IMP can be no if you wish tho evolve all spatial shear viscous tensor components"
              print *, "and obtain the other ones imposing the orthogonality condition"
              print *, "TRACE_IMP can also be zz if you wish to obtain the pi^zz component using also the null trace condition"
              print *, "No further value is allowed for TRACE_IMP"
              call exit(1)
            end if
          end if
        end if

        if(etam .lt. 0) then
          tilting=.false.
        else
          tilting=.true.
        end if

        if(derivatives_flag .eq. 1) then
           derivatives_out=.true.
        else
           derivatives_out=.false.
        end if
        if(flows_flag .eq. 1) then
           flows_out=.true.
           tp_anal=.true.
        else
           flows_out=.false.
           tp_anal=.false.
        end if

        if(restart_type .ne. 0) then

          if(mhd) then
            print *,"Restarting simulations is curently not possible in mhd simulations. Forcing restart_type=0."
            restart_type=0
            restart_type_option=3
          end if
          !if we dump the data to allow restarting of simulations, then we must use the output frequency for hypersurface
          !computation and output of primitive variables
          if(freeze_time_interval .ne. dtout) then
            if(freeze_time_interval .lt. dtout) then
              dtout=freeze_time_interval
            else
              freeze_time_interval=dtout
            end if

            print *,"Since you are dumping data to allow restarting of simulations, then the output frequencies for"
            print *,"freezeout hypersurface computation and for output printing of primitive variables have to be"
            print *,"the same. I just forced them to be: ",freeze_time_interval
          end if

          !we disable restarting in the case of Montecarlo simulations
          if(init_type .eq. GLAUBER_MONTECARLO) then
            if(restarted_run) then
              print *,"Sorry, but restarting is not possible in the case of Glauber-Montecarlo simulations..."
              stop
            else
              print *,"This is a Glauber-Montecarlo simulation, so restarting is disabled"
              restart_type=0
            end if
          end if

        end if !end if restart_type .ne. 0

      !here we read the data about the nucleus
      !first, we check that the specific file exists
      if(nuclei_data_file .ne. "none") then
       open(unit=28,status='OLD',file=nuclei_data_file, iostat=filerror, form='formatted')
        if (filerror .ne. 0) then
            write(*,*) nuclei_data_file//" cannot be opened, so I'm forced to quit!"
            close(28)
            call exit(1)
        end if
         read(28,'(a5)', IOSTAT=read_status) name_of_nucleus
          if(read_status .lt. 0) then
            write(*,*) "Sorry, but I can't read nucleus parameters... Please, check your parameters file and try again..."
            call exit(1)
         endif
         do while( name_of_nucleus .ne. nucleus )
             read(28,'(a5)', IOSTAT=read_status) name_of_nucleus
             if(read_status .lt. 0) then
             write(*,*) "Sorry, but I cannot find the parameters related to your selected nucleus..."
               write(*,*) "Please, check your parameters file and launch the program again."
               call exit(1)
             endif
         end do

         backspace(28) !I go back to read the whole line with data corresponding to my selected nucleus

         read(28,014) projmass,radius,delta,roze,zelectrons
        close(28)
      end if

      if(nrk_string .eq. "RK2") then
        nrk=RK2
      else if (nrk_string .eq. "RK3") then
        nrk=RK3
      else if (nrk_string .eq. "SSP") then
        if(.not. mhd) then
          write(*,*) "Sorry, but the IMEX schemes (SSP option) are currently implemented only for MHD simulations."
          write(*,*) "We set the integration algorithm to RK3, i.e. third order Runge-Kutta"
          nrk=RK3
          nrk_option=3
        end if 
        nrk=SSP
      else
        write(*,*) "Sorry, but I am not able to understand which integration algorithm you wish to use."
        write(*,*) "I will use the default, i.e. RK2"
        nrk=RK2
        nrk_option=3
      end if

      !check that the MHD solver in the cons2prim subroutine is appropriate for the eos
      if((init_type .eq. GUBSER_VISCOUS) .or. (init_type .eq. GUBSER_IDEAL)) then
        eqstate=1 !with Gubser flow we need a conformal EoS
      end if

      if(mhd) then
        if(eqstate .eq. 2) then !with the ideal gas law we need this algorithm
          algo_solver=4
          algo_solver_option=3
          write(*,*) "Relativistic ideal gas in use, I imposed algo_solver=4"
        end if

        if(algo_solver .eq. 5) then !with this algorithm we can use only the eqstate 1
          eqstate=1 
          eqstate_option=3
        end if

        if(rmhd) then
          if(algo_solver .ne. 7) then
            algo_solver = 7
            algo_solver_option=3
            write(*,*) "You are running a resistive MHD simulation, I imposed algo_solver=7."
          end if
          if(nrk .ne. SSP) then
             nrk=SSP
             algo_solver_option=3
             write(*,*) "You are running resistive MHD simulation, I nrk=SSP (IMEX-SSP time integration scheme)."
          end if
          
          if(sigma_el_cond .ne. 0) then
            ratio_chir_el=sigma_chiral_cond/sigma_el_cond
          else
            ratio_chir_el=0.d0
          end if
        end if !end if rmhd
      end if!end if mhd

      end if!end if pe0

      !for a parallel run now we spread the values of the variables among all other processors
      if(prl) then
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call MPI_Bcast( init_type, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( boundary_cond, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( derivatives_out, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( flows_out, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( tp_anal, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( coordinates, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( viscosity, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( bulkvis, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( mhd, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( rmhd, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( divclean, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( print_rho_comov, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( dump_initial_B, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( glm_alpha, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( B_th_press_ratio, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( edens_for_B_dump, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( pw, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( maxspeed, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( enable_smooth_viscosity, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( smooth_temp, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( tilting, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( etam, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( obtained, 2, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( bimpact, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( ah, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( events_start, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( events_stop, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( fixed_b, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( out_freeze, 1, MPI_Logical, 0, icomm, ierr)
        call MPI_Bcast( output_precision, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( projmass, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( radius, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( delta, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( roze, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( zelectrons, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( name_of_nucleus, 5, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( out_sel%density, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%vx, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%vy, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%vz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pressure, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%energy_density, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%temperature, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%entropy_density, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%bulk, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pitt, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pitx, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pity, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pitz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pixy, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pixz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%piyz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pixx, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%piyy, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%pizz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%v0, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%bx, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%by, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%bz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%ex, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%ey, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%ez, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%glm, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( out_sel%rc, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( ueta_coeff, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( input_edf, 120, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( input_B_field, 120, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( ext_cons_file, 120, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( ext_glissando_file, 120, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( freeze_type, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( freeze_value, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( freeze_time_interval, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( nconf, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( nbcoll, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( kind_of_collision, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( kappa, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( sig_mc, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( eprot, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( avg_events, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( min_participants, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( eos_tab_name, 120, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( eqstate, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( numerically, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( ecenter, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( enezero, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( rhocenter, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( deta, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( sigeta, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( eta_over_s, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( tau_pi_coeff, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( dtlog, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( dtout, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( maxdt, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( restart_type, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( cfl, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( recal, 5, MPI_CHARACTER, 0, icomm, ierr)
        call MPI_Bcast( nrk, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( rads, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( sigma_in, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( bimpact, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( ienentr, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( nx, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( ny, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( nz, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( x1, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( x2, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( y1, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( y2, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( z1, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( z2, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( tauzero, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( tmin, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( tmax, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( temp_end, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast(custom_output_directory,120,MPI_CHARACTER,0,icomm,ierr)
        call MPI_Bcast( sigma_el_cond, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( sigma_chiral_cond, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( ratio_chir_el, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( magfield_type, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Bcast( B_amplify_factor, 1, mpi_realtype, 0, icomm, ierr)
        call MPI_Bcast( algo_solver, 1, mpi_integer, 0, icomm, ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end if
         
      yb=log(rads/am)              
!     NB: sigma should be given in fm^2; 1 mb = 1/10 fm^2
      sigma_in=sigma_in/10.

      gamma_col=rads/(2.*am)
      beta_col=sqrt(1.-(1./gamma_col**2.))
      amassa=projmass

end subroutine common_get_all_parameters

! *****************************************************************************

subroutine bilinear_interpolation(D2input,nxin, nyin, xinput,yinput, nxout,nyout,xoutput,youtput,D2output)
  implicit none
  real(8), allocatable, intent(in), dimension(:,:) :: D2input
  real(8), allocatable, intent(in), dimension(:) :: xinput,yinput
  real(8), allocatable, intent(in), dimension(:) :: xoutput,youtput
  real(8), allocatable, dimension(:,:) :: D2output
  integer, intent(in) :: nxin, nyin, nxout, nyout 
  integer :: xin,yin,xout,yout
  real(8), allocatable, dimension(:,:) :: temparray
  integer :: allocate_result

  allocate(temparray(1:nxout,1:nyin),stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*) "Error, I can't allocate temparray, contained into common.f08"
    write(*,*)  "(source file common.f08)"
    call exit(1)
  end if
 
  do yin=1,nyin
    xout=1 
    do xin=1,nxin-1
      do while(xoutput(xout) .lt. xinput(1)) !values over the border of old grid are set equal to values at the border of old grid
        temparray(xout,yin)=D2input(xin,yin)
        xout=xout+1
        if (xout .gt. nxout) exit
      end do
      do while ((xoutput(xout) .ge. xinput(xin)) .and. (xoutput(xout) .le. xinput(xin+1)))
         !now we perform linear interpolation
         temparray(xout,yin)=D2input(xin,yin)+(D2input(xin+1,yin)-D2input(xin,yin))*(xoutput(xout)-xinput(xin))/&
         &(xinput(xin+1)-xinput(xin))
         xout=xout+1
         if (xout .gt. nxout) exit
      end do
      if (xout .gt. nxout) exit
      if(xin .eq. nxin-1) then
        do while(xout .le. nxout) !values over the border of old grid are set equal to values at the border of old grid
        temparray(xout,yin)=D2input(xin+1,yin)
        xout=xout+1
        end do
      end if
    end do
  end do 

  do xout=1,nxout
    yout=1 
    do yin=1,nyin-1
      do while(youtput(yout) .lt. yinput(1)) !values over the border of old grid are set equal to values at the border of old grid
        D2output(xout,yout)=temparray(xout,yin)
        yout=yout+1
        if (yout .gt. nyout) exit
      end do
      do while ((youtput(yout) .ge. yinput(yin)) .and. (youtput(yout) .le. yinput(yin+1)))
         !now we perform linear interpolation
         D2output(xout,yout)=temparray(xout,yin)+(temparray(xout,yin+1)-temparray(xout,yin))*(youtput(yout)-yinput(yin))/&
         &(yinput(yin+1)-yinput(yin))
         yout=yout+1
         if (yout .gt. nyout) exit
      end do
      if (yout .gt. nyout) exit
      if(yin .eq. nyin-1) then
        do while(yout .le. nyout) !values over the border of old grid are set equal to values at the border of old grid
        D2output(xout,yout)=temparray(xout,yin+1)
        yout=yout+1
        end do
      end if
    end do
  end do 

end subroutine bilinear_interpolation
! *****************************************************************************

end module common

! *****************************************************************************
