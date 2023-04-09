! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *         
! *                                                                           *         
! *  Copyright (C) 2015,2019,2021,2023 The ECHO-QGP team                      *
! *                                                                           *
! *  File: analysis/common.f90                                                *
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

  module constants
  implicit none
      real, parameter :: PIGRECO=3.14159265358979323846
      real, parameter :: HBARC=0.197
      
      INTEGER, PARAMETER :: maxpar=500 !max amount of particle accepted 
      real, parameter :: NOR_DENS=1.0/(2.0*PIGRECO*PIGRECO*HBARC*HBARC*HBARC)
      real, parameter :: TWOPIhbarc_cube = (2.0*PIGRECO*HBARC)**3
  end module constants
 ! !------------------------------------------------------------------  
  module common ! OVERALL
    use  constants
    implicit none

      integer dimension_flag
      integer viscosity_echo
      integer viscosity_spectra
      integer vorticity_flag
      integer freeze_type
      real  freeze_value
! 	directory and file names for the i/o
      character*32 inputdir,file,input,input1,outdir,output      
      integer LID_in, LID_out
!       number of cells in each direction
      integer nx(0:3)
      real dx(0:3)			!differentials for dV bjorken
      real xlim(2,0:3)		  	!lattice boundaries
      
      REAL Tfo, Efo
      ! !-----	ECHO parameters ------
	REAL, DIMENSION(:), ALLOCATABLE :: 	tau,x,y,eta 	!Bjorken coordinates
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: 	vx,vy,veta  	!Bjorken velocities
      ! !   raw grid on momentum space      
       integer  frozen_cells,  produced_particles
       real Energy_int,particle_energy
       
       integer npart
       character*24 name(maxpar) !particle name
       real, dimension(maxpar) :: m, mu
       	integer,dimension(maxpar) :: pdg_number , bf, g 
!       	real, dimension(:), allocatable :: pt, phi, rapidity
!       	real, dimension(:, :), allocatable :: mt
! 	integer npt, nphi, nrap
! 	real ptmin, ptmax, phimin, phimax, rapmin, rapmax, dphi, drap, dpt

	! need them during the setup reading
	integer BIN_pt, BIN_phi, BIN_y
	real mxv_pt,  mnv_pt, mxv_phi, mnv_phi, mxv_y ,mnv_y
	double precision  PTMAX_BOX
	double precision  YMAX_BOX 
	double precision  PHIMAX_BOX 
	
	integer seed_settings, trueseed
 end module common
 ! !------------------------------------------------------------------
 module common_MC !ONLY MC
  use  constants	
    implicit none
       
    integer, dimension(:), allocatable :: particle_index, particle_position
    integer, dimension(:), allocatable :: n_p_out
    
    real, dimension(:, :), allocatable :: part_coord, vol_elem, flu_mom, part_mom
 end module common_MC
 
! !------------------------------------------------------------------  
  module common_thermal !ONLY THERMAL
    use  constants
    implicit none
     	! !-----	PARTICLES parameters ------
      	real, dimension(:), allocatable :: pt, phi, rapidity
      	real, dimension(:), allocatable :: cosphi, senphi, cosphi2, senphi2,sencosphi
      	real, dimension(:, :), allocatable :: mt
	integer npt, nphi, nrap
	real ptmin, ptmax, phimin, phimax, rapmin, rapmax, dphi, drap, dpt
	
! 	INTEGER, PARAMETER :: maxpar=500 !max amount of particle accepted
       	
	real, dimension(maxpar) :: width, isospin  

	integer, dimension(maxpar) :: baryon, strange, bottom,charme, charge, n_decays
! 	
	integer mc_dec, daughters, mc1,  mc2 , mc3  ,mc4  ,mc5
	real ratio
! ! 	integer ID_start, ID_stop	
	real, dimension(:,:,:,:), allocatable :: spectra
		
! 	----- debugging purpose ----
	real, dimension(:,:,:,:), allocatable :: spectra_x, spectra_eta, spectra_tau, spectra_y
	real, dimension(:,:,:,:), allocatable :: s_0_x, s_0_eta, s_0_tau, s_0_y, spectra_0
! 	----- 
	real, dimension(:,:,:,:), allocatable :: spectra_m, spectra_b
! !     ------- VORTICITY HERE ----
	real, dimension(:,:,:,:,:), allocatable :: pola, pola_boost, pola_0, pola_boost_0
	real, dimension(:,:), allocatable :: u_der
	real, dimension(:), allocatable :: T_der
	real, dimension(:), allocatable :: exter ! de_mu(beta_nu)-de_nu(beta_mu)
	!indexes for exter 
	integer, parameter :: dtbx=1
	integer, parameter :: dtby=2
	integer, parameter :: dtbz=3
	integer, parameter :: dxby=4
	integer, parameter :: dxbz=5
	integer, parameter :: dybz=6


  end module common_thermal
! !------------------------------------------------------------------   
 subroutine check_file(Status, filename)
  implicit none
  integer, intent(in) :: Status
  character(len=*), intent(in) :: filename
  if (Status .ne. 0) then
    print *, "******  ERROR  ******"
    print *, "the file named:    " , filename
    print *, "cannot be opened. I am forced to quit!"
    
    print *, "filename test", filename
    
    call exit(1)
  end if
 end subroutine check_file
 ! !------------------------------------------------------------------   
 subroutine do_you_want_to_go_on()
  implicit none 
  character*1 switch   
  print *, "DO YOU WANT TO GO ON? [Y/n] (default is NO)"
  read (*,*) switch
  select case (switch)
    case('n')
      call exit 
    case('N')
      call exit 
    case('y')	  
      continue
    case('Y')
      continue
    ! 	      case ("")
    ! 		continue
    case default 
      call exit(9)
  end select
  return 
  end subroutine do_you_want_to_go_on
  
		
		
