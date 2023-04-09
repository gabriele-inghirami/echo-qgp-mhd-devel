! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2016,2021 The ECHO-QGP team                           * 
! *                                                                           *
! *  File: analysis/io.f90                                                    *
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


module io
  use constants
  implicit none 
  integer, dimension(:), allocatable :: id_list
  integer antibaryons, chempot, ID_start, ID_stop, los, part_list
  contains
! !******************************************************************
! !******************************************************************
! !***********************   INPUT  ROUTINES   **********************
! !******************************************************************  
  subroutine io_read_setup()
! ! !  This routine reads the setup declared in settings.txt
  use common 
!   use common_PaGe
  use common_thermal, only: npt, nphi, nrap, ptmin, ptmax, phimin, phimax, rapmin, rapmax
  use init
  implicit none
  integer ios, i, AllocateStatus

  
  008	format(10x,i1)		!INT 1 
  009	format(10x,i2)		!INT 2
  010	format(10x,i6)		!INT 6
  011	format(10x,i10)		!INT 10 
  
  111	format(10x,f8.4) 	!FLOAT
  113	format(10x,a32)		!STRING 32
  
  open(unit=3,status='old',file='settings.txt',iostat=ios)
  call check_file(ios, 'settings.txt')
  
  read(3,*) ! 2 dummy lines
  read(3,*)
  
  read(3,009) dimension_flag
  read(3,008) viscosity_echo
  read(3,008) viscosity_spectra  
  read(3,008) vorticity_flag  
  select case(viscosity_spectra)
    case(0)
      continue
    case(1)
      if (viscosity_echo /= 1) then
	print *, "************************* error ********************************"
	print *, "***************** I am going to quit ***************************"
	print *, "I am sorry, but you cannot take into account"
	print*, "viscosity corrections to the distribution function"
	print*, "if you do not perform a viscous-hydro simulation"
	print*, "in the file settings.txt, your choice was"	
	print *, "visco_hyd=", viscosity_echo
	print *, "visco_spe=", viscosity_spectra
	call exit	
      endif 
    case default
      call exit
  end select
  read(3,010) nx(0)
  read(3,*)
  
  read(3,010) npt
  read(3,111) ptmin
  read(3,111) ptmax
  read(3,010) nphi
  read(3,111) phimin
  read(3,111) phimax
  read(3,010) nrap
  read(3,111) rapmin
  read(3,111) rapmax
  
  read (3,008) los
  if (los >3) then
    print *, "***************** error ***************************"
    print *, "You choose a weird option for the flag listorseq", los
    print *, "Please choose between 0 = all particle, 1= sequence of particle, 2-list of particle, 3-stable particles "
    call exit(6)
  endif
  read (3,011) ID_start
  read (3,011) ID_stop
  read (3,010) part_list

  if (los /= 2 ) then
    if (part_list /= 0) then
      print *, "***************** error ***************************"
      print *, "You now have"
      print *, ""
      write(*,*) "part_list=",part_list,"                        ! particles in the list ()"
      print *, ""
      if (los == 0) then
	print *, "You choose to produce ALL the particle available in ../eos_data/pdglist.txt"
      else if (los ==1 ) then
	print *, "You choose to produce a sequence of the particle available in ../eos_data/pdglist.txt"
	print *, "from particle number", ID_start, "to particle number", ID_stop
      else
      	print *, "You choose to produce the stable particles: pions, kaons, eta, omega, proton, neutron, Lambda, Sigma, Xsi "
      endif
      print *, "To perform the calculation correctly plase set the particle in the list to be 0, as following:"
      print *, ""
      write(*,'(a11, a60)') "part_list=0","                        ! particles in the list ()"
      print *, ""
      print *, "And be sure that the following line is "
      print *, ""
      write(*,'(a11, a60)')  "antibar..=","                         ! 1=on 0=off"
      print *, ""
      print *, "Please run the program again, I am quitting... Bye!"
      call exit(6)
    endif
  endif
  
  ALLOCATE ( id_list(part_list), STAT = AllocateStatus)
  if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate id_list in read_setup (io-3D.f08) ***"
  do i=1, part_list
    read (3,'(i10)') id_list(i)	
  end do
  
  if (los == 3) then
    DEALLOCATE ( id_list, STAT = AllocateStatus)
    if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate id_list in read_setup (io-3D.f08) ***"
! !     STABLE PARTICLES
    part_list=18
    ALLOCATE ( id_list(part_list), STAT = AllocateStatus)
    if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate id_list in read_setup (io-3D.f08) ***"
    id_list(1)= 211	!pi +
    id_list(2)= 111	!pi 0
    id_list(3)= -211	!pi -
    id_list(4)= 321	!K +
    id_list(5)= -321	!K -
    id_list(6)= 311	!K 0
    id_list(7)= -311	!K 0 bar
    id_list(8)= 221   	!eta
    id_list(9)= 223	!omega
    id_list(10)=2212	!p
    id_list(11)=2112	!n
    id_list(12)=3122	!Lambda
    id_list(13)=3222	!Sigma +
    id_list(14)=3212	!Sigma 0
    id_list(15)=3112	!Sigma -  
    id_list(16)=3322  	!Xsi 0
    id_list(17)=3312    !Xsi -
    id_list(18)=4       !charm quark +
  endif 
  
  read (3,008) antibaryons
  read (3,008) chempot
  
  read(3,*)

  read(3,113) inputdir,outdir 
  LID_out=0
  LID_in=0
  LID_out=index(outdir, ' ')-1
  LID_in=index(inputdir, ' ')-1
  if (LID_out==0) then
    LID_out=64
  endif
  if (LID_in==0) then
    LID_in=64
  endif

  read(3,113) file
  call check_file_folder_existence(outdir, LID_out, inputdir, file, LID_in)
  
  read(3,*) 
  read(3, 011) seed_settings 
  read(3, 111) PTMAX_BOX
  read(3, 111) PHIMAX_BOX
  if (PHIMAX_BOX == 0) then
    PHIMAX_BOX=(2.0*PIGRECO)
  endif 
  if (PHIMAX_BOX == 0.0) then
    PHIMAX_BOX=(2.0*PIGRECO)
  endif 
  read(3, 111) YMAX_BOX
  read(3,*) 
  read(3, 111) mxv_pt	!(pt
  read(3, 111) mnv_pt 	!(pt
  read(3, 111) mxv_phi	!(phi)
  read(3, 111) mnv_phi	!(phi)
  read(3, 111) mxv_y 	!(y) 
  read(3, 111) mnv_y 	!(y) 
  read(3, 010) BIN_pt 	!(pt)
  read(3, 010) BIN_phi 	!(phi)
  read(3, 010) BIN_y 	!(y)  

  
  close(3)
! ! !-----------------------------------------------  
  call readgrid()
! ! !-----------------------------------------------    
    if (nx(3) < 3 .and. dimension_flag> 2) then 
      print *,     "-----------------------------------------------"
      print *, "ERROR: this is a 3+1D simulation, but i have just less than 3 cells along the 3rd direction"
      print *, "I am quitting ... check it again"
      call exit
    else if (nx(3)>1 .and. dimension_flag<3) then
      print *,     "-----------------------------------------------"
      print *, "ERROR: this is a 2+1D simulation, but i have just more than 1 cell along the 3rd direction"
      print *, "I am quitting ... check it again"
      call exit    
    endif

    select case (dimension_flag)
	case(1)
	  print *, "1+1D still not implemented"
	case(2)
	  print*, "-----------------------------------------------"
	  print *, "THIS IS A (2+1)-D SIMULATION"
	  print *, "Since this a (2+1)-D simulation, I set nrap=1 and rapidity=0"
          nrap=1          
          rapmin=0.0
          rapmax=0.0
          BIN_y=1
          call allocate_momenta()
	case(3)
!!!!	if the range is symmetric, than we have automatically y=0 in the middle index
	  print*, "-----------------------------------------------"
	  print *, "THIS IS A (3+1)-D SIMULATION"
          if (nrap -((nrap/2)*2) .ne. 1 ) then
	    nrap=nrap + 1
	    print *, "I changed nrap from ", nrap-1, "to ", nrap, "so that "
	  else 
	    print *, " HINT :"	    
	  endif
	  print *, "If you have chosen the rapidity range to be symmetric,"
	  print *, "than you have automatically y=0 in the middle index"
	  print *,  "which is: " , nrap/2+1
	  call allocate_momenta()	  
	  
	case default
	  print *, "You can only choose 1, 2 or 3 for the dimension of the simulation"
	  print *, "instead your choice was:  ", dimension_flag
	  print *, "I am going to quit!"
	  call exit
	end select
	

108	format(a30,3x,i3)
109	format(a30,3x,a22)
107	format(a30,f10.3,a3,f10.3,a3,f10.3)
	print*, "-----------------------------------------------"
        print*, "---------   Report of settings.txt  -----------"
	
	write (*,108) 'total amount of time-steps',nx(0)
	write (*,108) 'Viscosity switch ECHO', viscosity_echo
	write (*,108) 'Viscosity switch spectra', viscosity_spectra
	write (*,108) 'Vorticity switch ', vorticity_flag
	write (*,108) 'Points for the transv. p', npt
	write (*,107) 'min/Max value for p_t ', ptmin, '/', ptmax
	write (*,108) 'Points for the polar angle', nphi
	write (*,107) 'min/Max value for phi ', phimin, '/', phimax
	write (*,108) 'Points for the rapidity', nrap
	write (*,107) 'min/Max value for rapidity ', rapmin, '/', rapmax
	write (*,109) 'input directory',inputdir
	write (*,109) 'output directory',outdir
	write (*,109) 'input filename',file 
	print*, "-----------------------------------------------"	
	write (*,*) " seed of the simulation", seed_settings
	write (*,107) "	y/phi/PT BOX for the MC ",  YMAX_BOX, '/',PHIMAX_BOX,'/', PTMAX_BOX
	write (*,107) "min/max for pt", mnv_pt, '/',mxv_pt
	write (*,107) "min/max for phi", mnv_phi, '/',mxv_phi
	write (*,107) "min/max for y", mnv_y, '/',mxv_y
	
	return
  end subroutine io_read_setup
! !******************************************************************

  subroutine readgrid()
    use common, only: nx, dx, xlim, inputdir, tau,x, y,eta
    use init
    implicit none	
    integer lid, i, nout_counter
    integer nxx(0:3)
    integer filerror, AllocateStatus
    real, dimension(:), allocatable :: tempo
    
   
    ALLOCATE ( tempo(10000), STAT = AllocateStatus)
    if (AllocateStatus /= 0)   STOP 	"ERROR:*** Cannot allocate tempo in readgrid ***" 


	  lid=index(inputdir,' ')-1
	  !we read the grid info from ../out/grid_summary.dat
	  open(unit=1,status='old',file=inputdir(1:lid)//'grid_summary.dat',form='formatted', iostat=filerror)
	  call check_file(filerror, inputdir(1:lid)//'grid_summary.dat')
	  read(1,*) nx(1),nx(2),nx(3) !lattice total amount of dumped cells in each  direction
	  read(1,*) dx(1),dx(2),dx(3) !lattice interstitial space
	  read(1,*) (xlim(1,i),i=1,3) !lattice dimensions
	  read(1,*) (xlim(2,i),i=1,3)
	  close(1)

	  lid=index(inputdir,' ')-1
	  open(unit=2,status='old',file=inputdir(1:lid)//'grid.dat',form='formatted', iostat=filerror)
	  call check_file(filerror, inputdir(1:lid)//'grid.dat')
! 
	  open(unit=4,status='old',file=inputdir(1:lid)//'time.dat',form='formatted', iostat=filerror)
	  call check_file(filerror, inputdir(1:lid)//'time.dat')
! 
	  nout_counter=1
	  DO WHILE (nout_counter.lt.10000)
	    read(4,"(I5,f16.9)", end=150) nout_counter, tempo(nout_counter)
	  END DO
150	  CONTINUE 

115       format(a50, 2x, i6)
	  if (nout_counter .ne. nx(0)) then 
	    print*, "-----------------------------------------------"
	    write (*, 115) "the lines in the file time.dat are", nout_counter
	    write (*, 115) "in the settings file, instead you gave me", nx(0)
	    if (nx(0)>nout_counter) then
	      print *,"I need to use the time steps in the settings file. Please check it again"
	      print *,"I have to quit now... BYE!"
	      call exit
	    else
	      if (nx(0)==0) then 
		nx(0)=nout_counter
		write (*, 115) "I am fixing the input file automatically to ", nx(0)
	      else
		write (*, 115) "I will use " , nx(0)
		print *, "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
		call do_you_want_to_go_on()
	      endif
	    endif
	  endif
! 
! 	! !-----------------------------------------------
	  ALLOCATE ( tau(nx(0)), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  ALLOCATE ( x(nx(1)), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0) 	STOP "*** Not enough memory ***"
	  ALLOCATE ( y(nx(2)), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  ALLOCATE ( eta(nx(3)), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0) 	STOP "*** Not enough memory ***"	
! 	! !-----------------------------------------------
	
	do i=1, nx(0)
	  tau(i)=tempo(i)
	end do
	dx(0)=tau(2)-tau(1)
	xlim(1,0)=tau(1)
	xlim(2,0)=tau(nx(0))
	close(4)
! 	  
	read(2,*) nxx(1),nxx(2),nxx(3)
	if (nxx(1).ne.nx(1).or.nxx(2).ne.nx(2).or.nxx(3).ne.nx(3)) then
	    print *, "*** I am forced to quit *** "
	    print *, "ERROR: I have found different values of the lattice size in the files grid.dat and grid summary.dat!!"
	    print *, "grid.dat        ",( nx(i),i=1,3)
	    print *, "grid_summary.dat",(nxx(i),i=1,3)
	    print *, "*** I am forced to quit *** "
	  call exit(8)
	endif
	read(2,*) (  x(i),i=1,nx(1))
	read(2,*) (  y(i),i=1,nx(2))
	read(2,*) (  eta(i),i=1,nx(3))
    deallocate(tempo)
    return
  end subroutine readgrid
! ! !******************************************************************  
! 
  subroutine io_read_particlelist()
	use common
	use common_thermal
	use init
	implicit none
	integer filerror, i, AllocateStatus, ipt, ipart
	integer j, k, check_los
	integer dummy_int
	
	real temperature(2)
	real temp_mu(2,maxpar)
	real weight
	integer check_mu
	real Tfo_MeV
!   
	real, dimension(maxpar) :: R_m, R_mu, R_width, R_isospin  
	integer,dimension(maxpar) :: R_pdg_number , R_bf, R_g 
	integer, dimension(maxpar) :: R_baryon, R_strange, R_bottom, R_charme, R_charge, R_n_decays
! 	
	integer R_mc_dec, R_daughters, R_mc1,  R_mc2 , R_mc3  ,R_mc4  ,R_mc5, R_npart
	real R_ratio

	integer hypercharge
	character*24 R_name(maxpar) !particle name      

	  print*, "-----------------------------------------------"
	  print*, "-------          particles info            --------"

	if (antibaryons .eq. 1) then
	  write(*,'(a30)') "I will produce antibaryons"
	else
	  write(*,'(a30)') "I will NOT produce antibaryons"
	endif
! 
	Tfo_MeV=1000.0*Tfo
	check_mu=0

	open(unit=30,status='old',file='../eos_data/pdglist.txt',form='formatted', iostat=filerror)
	call check_file(filerror, '../eos_data/pdglist.txt')	

! 
! ! skip the first 35 lines
	do i=1, 35
	read (30, *) 
	end do
	i=1
	do while (i.lt.maxpar) 
	  read (30, *, end=166) ! the first line is empty
	  read (30, *, end=166) R_pdg_number(i)
! 	  PRINT *, R_pdg_number(i)
	  READ (30, *, END=166) R_name(i)
! 	  PRINT *, R_name(i)
	  READ(30,*,END=166) R_m(i), R_width(i), R_g(i), R_baryon(i), R_strange(i), R_charme(i), R_bottom(i), &
			   & R_isospin(i), R_charge(i), R_n_decays(i)
! 	  PRINT *,  R_m(i), R_width(i), R_g(i), R_baryon(i), R_strange(i), R_charme(i), R_bottom(i), &
! 			   & R_isospin(i), R_charge(i), R_n_decays(i)
	  do j=1,R_n_decays(i)
	    read(30,*) R_mc_dec, R_daughters, R_ratio, R_mc1,  R_mc2 , R_mc3,  R_mc4, R_mc5
! 	    PRINT *, R_mc_dec, R_daughters, R_ratio, R_mc1,  R_mc2 , R_mc3,  R_mc4, R_mc5
	  end do
	  
	  if (R_baryon(i).NE.0) then
	      R_bf(i)=-1
	    else 
	      R_bf(i)=1
	    endif    	  
	    hypercharge= R_strange(i)+R_charme(i)+R_bottom(i) + R_baryon(i) !top quark out of range
	    R_isospin(i)=2*R_charge(i)-hypercharge
	    i=i+1
	  end do
	  print *, "you have more particle than expected, change maxpar in the common file" 
	  call exit
 166    continue
	R_npart=i-1
! 	WRiTE(*,2001) "INDEX","PDG","NAME","MASS","DEG", "B", "e", "S", "b/f", "mu"
! 	PRINT *, "---------------------------------------------------"
	CLOSE(30)

	if (chempot .eq. 1) then 
	  open(unit=31,status='old',file='../eos_data/chemical_potential.txt',form='formatted', iostat=filerror)	
	    if(0.ne.filerror) then 
	      print *, "*** I am forced to quit *** "
	      print *, "cannot find the file ../eos_data/chemical_potential.txt"
	      print *, "*** I am forced to quit *** "
	      call exit(8)
	    end if
  ! ! 	print *, Tfo
	  R_mu=0.0
	  check_mu=0
	  read (31,*) temperature(2), (temp_mu(2, j), j=1, R_npart)
	  do while (check_mu.lt.1) 
	    read (31,*, end=167) temperature(1), (temp_mu(1, j), j=1, R_npart)
	    if ((temperature(1)-Tfo_MeV)*(Tfo_MeV-temperature(2)) .GE. 0.0) then 
	      check_mu=2
	      weight=(temperature(2)-Tfo_MeV)/(temperature(2)-temperature(1))
	      if (weight<0.0) then 
		print *, "error (weight<0.0)  in io-3D.f08 subroutine read_particlelist ", weight
		call exit (5)
	      endif
	      do j=1, R_npart
		R_mu(j)=temp_mu(2,j)-weight*(temp_mu(2,j)-temp_mu(1,j))
	      end do
	    endif
	    do j=1, R_npart
	      temp_mu(2, j)=temp_mu(1, j)
	    end do
	    temperature(2)=temperature(1)
	  end do
167	  continue
	  CLOSE(31)
	else if (chempot .eq. 0 ) then 
	  R_mu(j)=0.0
	  check_mu=2
	else 
	  print *, "something's wrong, in settings.txt i found chempot=",chempot, " but it can only be 0 or 1"
	  print *, "I am going to quit now"
	  call exit	  
      	endif 

	if (check_mu==0) then
	  print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	  print *, "Hey, something is wrong in routine read_particlelist,"
	  print *, "I cannot interpolate the chemical_potential."
	  print *, "please check again the file ../eos_data/chemical_potential.txt"
	  print *, "and the file ../settings.txt for the Freezeout temperature"
	  print *, "I AM ABORTING THE SIMULATION"
	  call exit (5)
	endif 
! ! 
	select case(los) 
	  case (0) 
	    print *,  "I will produce all the particles in the file ../eos_data/pdglist.txt"
	    m=		R_m
	    mu=		R_mu
	    width=	R_width
	    isospin=	R_isospin  
	    pdg_number=	R_pdg_number
	    bf=		R_bf
	    g=		R_g
	    baryon=	R_baryon
	    strange=	R_strange
	    bottom=	R_bottom
	    charme=	R_charme
	    charge=	R_charge
	    n_decays=	R_n_decays
	    npart=	R_npart
	    name=	R_name
	  case (1) 
	    print *,  "I will produce the sequence of particles in the file ../eos_data/pdglist.txt"
	    print *,  "between the ID", ID_start, "and", ID_stop
	    check_los=0
	    do i=1, R_npart
	      if (R_pdg_number(i) == ID_start) then
		check_los=i
	      else if (R_pdg_number(i) == ID_stop) then
		dummy_int=i
	      endif 
	    end do
! 	    print *, ID_start, ID_stop
	    if (dummy_int .le. check_los) then
	      dummy_int=ID_start
	      ID_start=ID_stop
	      ID_stop=dummy_int
	    endif

! 	    print *, ID_start, ID_stop
! 	    call exit
	    k=1
	    i=1
	    check_los=0
	    do i=1, R_npart
	      if (R_pdg_number(i) == ID_start) then
		check_los=1
	      endif
	      if (check_los>0 .and. check_los<2 ) then
		  m(k)=		R_m(i)
		  mu(k)=	R_mu(i)
		  width(k)=	R_width(i)
		  isospin(k)=	R_isospin(i)
		  pdg_number(k)=R_pdg_number(i)
		  bf(k)=	R_bf(i)
		  g(k)=		R_g(i)
		  baryon(k)=	R_baryon(i)
		  strange(k)=	R_strange(i)
		  bottom(k)=	R_bottom(i)
		  charme(k)=	R_charme(i)
		  charge(k)=	R_charge(i)
		  n_decays(k)=	R_n_decays(i)		
		  name(k)=	R_name(i)
		  mu(k)=	R_mu(i)
		  k=k+1
		  select case (antibaryons)
		    case (0)
		      continue
		    case (1)
		    if (baryon(k-1).NE.0) then
		      m(k)=         m(k-1)
		      mu(k)=	    -mu(k-1)
		      width(k)=	     width(k-1)
		      isospin(k)=   -isospin(k-1)
		      pdg_number(k)=-pdg_number(k-1)
		      bf(k)=	    bf(k-1)	                      
		      g(k)=	    g(k-1)
		      baryon(k)=    -baryon(k-1)
		      strange(k)=   -strange(k-1)
		      bottom(k)=    -bottom(k-1)
		      charme(k)=    -charme(k-1)
		      charge(k)=    -charge(k-1)
		      n_decays(k)=  -n_decays(k-1)
		      name(k)=	     'Anti_'//name(k-1)
		      mu(k)=	    -mu(k-1)
		      k=k+1
		    endif
		    case default
		      print *,  "You can just select 0 or 1 to for the automatic production of antibaryons"
		      print *,  "You selected ", antibaryons
		      call exit (5)
		  end select    
  	      endif
	      if (R_pdg_number(i) == ID_stop) then
		check_los=3
	      endif  
	    end do
	    npart=k-1
	  case(2) 
	    print *,  "I will produce only the particles listed in the settings.txt file"    

	    k=1
	    do i=1, R_npart
	      do j=1, part_list
		if(id_list(j)==R_pdg_number(i)) then
		  m(k)=		R_m(i)
		  mu(k)=		R_mu(i)
		  width(k)=	R_width(i)
		  isospin(k)=	R_isospin(i)
		  pdg_number(k)=	R_pdg_number(i)
		  bf(k)=		R_bf(i)
		  g(k)=		R_g(i)
		  baryon(k)=	R_baryon(i)
		  strange(k)=	R_strange(i)
		  bottom(k)=	R_bottom(i)
		  charme(k)=	R_charme(i)
		  charge(k)=	R_charge(i)
		  n_decays(k)=	R_n_decays(i)		
		  name(k)=	R_name(i)
		  mu(k)=	R_mu(i)
		  k=k+1
		  select case (antibaryons)
		    case (0)
		      continue
		    case (1)
		    if (baryon(k-1).NE.0) then
		      m(k)=         m(k-1)
		      mu(k)=	    -mu(k-1)
		      width(k)=	     width(k-1)
		      isospin(k)=   -isospin(k-1)
		      pdg_number(k)=-pdg_number(k-1)
		      bf(k)=	    bf(k-1)	                      
		      g(k)=	    g(k-1)
		      baryon(k)=    -baryon(k-1)
		      strange(k)=   -strange(k-1)
		      bottom(k)=    -bottom(k-1)
		      charme(k)=    -charme(k-1)
		      charge(k)=    -charge(k-1)
		      n_decays(k)=   n_decays(k-1)
		      name(k)=	     'Anti_'//name(k-1)
		      mu(k)=	    -mu(k-1)
		      k=k+1
		    endif
		    case default
		      print *,  "You can just select 0 or 1 to for the automatic production of antibaryons"
		      print *,  "You selected ", antibaryons
		      call exit (5)
		  end select    
		endif
	      end do
	    end do
	    npart=k-1
	  case(3) 
	    print *,  "I will produce only the stable particles"    

	    k=1
	    do i=1, R_npart
	      do j=1, part_list
		if(id_list(j)==R_pdg_number(i)) then
		  m(k)=		R_m(i)
		  mu(k)=	R_mu(i)
		  width(k)=	R_width(i)
		  isospin(k)=	R_isospin(i)
		  pdg_number(k)=R_pdg_number(i)
		  bf(k)=	R_bf(i)
		  g(k)=		R_g(i)
		  baryon(k)=	R_baryon(i)
		  strange(k)=	R_strange(i)
		  bottom(k)=	R_bottom(i)
		  charme(k)=	R_charme(i)
		  charge(k)=	R_charge(i)
		  n_decays(k)=	R_n_decays(i)		
		  name(k)=	R_name(i)
		  mu(k)=	R_mu(i)
		  k=k+1
		  select case (antibaryons)
		    case (0)
		      continue
		    case (1)
		    if (baryon(k-1).NE.0) then
		      m(k)=         m(k-1)
		      mu(k)=	    -mu(k-1)
		      width(k)=	     width(k-1)
		      isospin(k)=   -isospin(k-1)
		      pdg_number(k)=-pdg_number(k-1)
		      bf(k)=	    bf(k-1)	                      
		      g(k)=	    g(k-1)
		      baryon(k)=    -baryon(k-1)
		      strange(k)=   -strange(k-1)
		      bottom(k)=    -bottom(k-1)
		      charme(k)=    -charme(k-1)
		      charge(k)=    -charge(k-1)
		      name(k)=	     'Anti_'//name(k-1)		      
		      k=k+1
		    endif
		    case default
		      print *,  "You can just select 0 or 1 to for the automatic production of antibaryons"
		      print *,  "You selected ", antibaryons
		      call exit (5)
		  end select    
		endif
	      end do
	    end do
	    npart=k-1	    
	  end select

 2000 FORMAT(I6,I10,2X, A25,F10.6,I5,4I3,3x, F10.6)
 2001 FORMAT(A6,A10,3x,A4,22x,A4,7x, A3, 1x, A1, 2(2x, A1),2x, A3, 2x, A3)
 2002 FORMAT(A10,A7,3x,A4,22x,A4,6x, A3, 1x, A3, 3x, A5, 2x, A3)
 2003 FORMAT(I6,I10,2X, A25,F10.6,I5,I3,3x, F10.6)

       	WRiTE(*,*) "Particle list and relative features"
      	WRiTE(*,2002) "INDEX","PDG","NAME","MASS","DEG", "b/f", "mu"
	do i=1, npart	
	write(*,2003) i,pdg_number(i),name(i),m(i),g(i), bf(i), mu(i)
	end do
	PRINT *, "---------------------------------------------------"
	ALLOCATE ( mt(npart, npt), STAT = AllocateStatus)
	IF (AllocateStatus /= 0) 	STOP "*** Not enough memory ***"
    	do ipt=1, npt
	  pt(ipt)=ptmin+(ipt*1.0-0.5)*dpt
	  do ipart=1, npart
	    mt(ipart,ipt)=sqrt(((pt(ipt))**2) + ((m(ipart))**2))
	  end do
	end do


  return
  end subroutine io_read_particlelist
! !******************************************************************




! !******************************************************************
! !******************************************************************
! !***********************   OUTPUT  ROUTINES   *********************
! !******************************************************************  
  subroutine io_build_particle_filename(idx_part, string, fname, lix)
    use common, only: outdir, LID_out, name, pdg_number
    
    integer, intent (in) :: idx_part
    character*32, intent(inout) :: string
    character*128, intent(out) :: fname
    integer, intent(out) :: lix
    
    character*32 particle  
    character*3 segno
    character*11 pdg_n
    integer lin, lis, lip
    
    lis=index(string, ' ')-1
    
    particle=name(idx_part)
    lip=index(particle, ' ')-1
    
    if (pdg_number(idx_part)<0) then
      segno='_-_'
    else if (pdg_number(idx_part)>0) then
      segno='_+_'
    else 
      segno='_0_'
    endif    
    
    if (abs(pdg_number(idx_part))<10) then !A
      write(pdg_n, '(7I1, I1,A3)')  0,0,0,0,0,0,0,abs(pdg_number(idx_part)), segno(1:3)
    else if (abs(pdg_number(idx_part))<100) then !A
      write(pdg_n, '(6I1, I2,A3)')  0,0,0,0,0,0,abs(pdg_number(idx_part)), segno(1:3)
    else if (abs(pdg_number(idx_part)) .ge. 100 .and. abs(pdg_number(idx_part)) <1000 )  then
      write(pdg_n, '(5I1, I3,A3)')  0,0,0,0,0,  abs(pdg_number(idx_part)), segno(1:3)
    else 
      print *, ' '
      if   (abs(pdg_number(idx_part)) .ge. 1000 .and. abs(pdg_number(idx_part)) <10000 ) then !B
	write(pdg_n, '(4I1, I4,A3)')  0,0,0,0,    abs(pdg_number(idx_part)), segno(1:3)
      else if  (abs(pdg_number(idx_part)) .ge. 10000 .and. abs(pdg_number(idx_part)) <100000 ) then
	write(pdg_n, '(3I1, I5,A3)')  0,0,0,      abs(pdg_number(idx_part)), segno(1:3)
      else 
	print *, ' '
	if (abs(pdg_number(idx_part)) .ge. 100000 .and. abs(pdg_number(idx_part)) <1000000 )  then !C
	  write(pdg_n, '(2I1,I6,A3)') 0,0,         abs(pdg_number(idx_part)), segno(1:3)
	else if (abs(pdg_number(idx_part)) .ge. 1000000 .and. abs(pdg_number(idx_part)) <10000000 )  then
	  write(pdg_n, '(I1,I7,A3)')  0,        abs(pdg_number(idx_part)), segno(1:3)
	else 
	  if (abs(pdg_number(idx_part)) .ge. 10000000 .and. abs(pdg_number(idx_part)) <100000000 )  then !D
	    write(pdg_n, '(I8,A3)')  0,        abs(pdg_number(idx_part)), segno(1:3)
	  else   
	    print *, 'ooops... there was an error in the print_spectra routine.'
	    print *, 'I am sorry, I have to quit. Please contact the developer.'	
! 	    remember that pdg_n is 8+3 here. if you have a bigger pdg number remember to change the strings too
	    call exit (2)
	  endif !D
	endif  !C
      endif	 !B
    endif 	 !A
        
!      128 = 32 // 1 // 32 // 1 // 8+3 // 32 // 4 (=113 < 128)
    fname=outdir(1:LID_out)//'/'//string(1:lis)//'_'//pdg_n(1:11)//particle(1:lip)//'.txt'
    lix=index(fname, ' ')-1
    
    string='                                '
    return
  end subroutine io_build_particle_filename
! !******************************************************************
subroutine io_print_setup(rep, title)
 use common
 use common_MC
 use common_thermal
 implicit none 
 
 integer, intent (in) :: rep
 character*8, intent(in) :: title
 integer filerror, i
 
 2003 FORMAT(I6,I10,2X, A25,F10.6,I5,I3,3x, F10.6)
 108	format(a30,3x,i3)
 109	format(a30,3x,a22)
 107	format(a30,f10.3,a3,f10.3)
 
 select case (rep)
  case(0) !!THERMAL
    open(unit=21,status='replace',file=outdir(1:LID_out)//'report_thermal.txt', form='formatted', iostat=filerror)
    call check_file(filerror, outdir(1:LID_out)//'report_thermal.txt')
  case(1) !!PAGE
   open(unit=21,status='replace',file=outdir(1:LID_out)//'report_PaGe.txt', form='formatted', iostat=filerror)
   call check_file(filerror, outdir(1:LID_out)//'report_PaGe.txt')
  case(2) !!PAMOGE
   open(unit=21,status='replace',file=outdir(1:LID_out)//'report_PaMoGe.txt', form='formatted', iostat=filerror)   
   call check_file(filerror, outdir(1:LID_out)//'report_PaMoGe.txt')
  case(3) !!MERGE
   open(unit=21,status='replace',file=outdir(1:LID_out)//'report_merge.txt', form='formatted', iostat=filerror)   
   call check_file(filerror, outdir(1:LID_out)//'report_merge.txt')
 end select
 

  write (21,*) "###############################################"
  write (21,*) "------   ", title, "   ------" 
  write (21,*) "-----------------------------------------------"
  write (21,*) "---------   Report of settings.txt  -----------"
  write (21,108) 'dimensionality of the simulation', dimension_flag
  write (21,108) 'total amount of time-steps',nx(0)
  write (21,108) 'Viscosity switch ECHO', viscosity_echo
  write (21,108) 'Viscosity switch SPECTRA', viscosity_spectra
  write (21,108) 'Points for the transv. p', npt
  write (21,107) 'min/Max value for p_t ', ptmin, '/', ptmax
  write (21,108) 'Points for the polar angle', nphi
  write (21,107) 'min/Max value for phi ', phimin, '/', phimax
  write (21,108) 'Points for the rapidity', nrap
  write (21,107) 'min/Max value for rapidity ', rapmin, '/', rapmax
  write (21,109) 'input directory',inputdir
  write (21,109) 'output directory',outdir
  write (21,109) 'input filename',file 
  write (21,*) "-----------------------------------------------"	
  write (21,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write (21,*) " seed setup", seed_settings
  write (21,*) " seed of the simulation", trueseed  
  write (21,*) " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  write (21,107) " y/PT BOX for the MC ",  YMAX_BOX, '/', PTMAX_BOX
  write (21,107) "min/max for pt", mnv_pt, '/',mxv_pt
  write (21,107) "min/max for phi", mnv_phi, '/',mxv_phi
  write (21,107) "min/max for y", mnv_y, '/',mxv_y      
  write (21,*) "-----------------------------------------------"
  write (21,*) "-------   Particle list and relative features     --------"  
  write (21,108) 'chemical potential (on/off)', chempot
  write (21,'(7(a7,2x))') "INDEX","PDG","NAME","MASS","DEG", "b/f", "mu"  
  do i=1, npart	
    write(21,2003) i,pdg_number(i),name(i),m(i),g(i), bf(i), mu(i)
  end do
  close(21)  
 return      
end subroutine io_print_setup  
! !******************************************************************
! !******************************************************************
  subroutine io_print_spectra()
    use common
    use common_thermal
    implicit none

! !     proper name of the outout file, es: 'spectra_v2_'
    character*32 stringa32
    character*8 pdg_n
    character*32 particle
    character*3 segno
    integer lip ! particle name index
    ! 	------- 
    
    integer iphi, irap, ipart, ipt, ivort
    integer ind_rap !check on the y=0 index
    integer fstat ! status flag
! ! 	Vorticity and polarization
    real, dimension(0:3) :: pi0, pi0_0
! 	------- 
    
    ! integer for files and files
    !outdir
    integer lid 		
    ! phi_spectra		 (22)
    character*128 filename
    integer lif
    ! spectra and v2 (vs pt)	 (26)
    character*128 filename_int
    integer lif_i
    ! spectra and v1 (vs y)      (21)
    character*128 filename_ys
    integer lif_ys
    ! DEBUGGING ------
    ! components of spectra and v2 (vs pt)	
    character*128 filename_components		!(11)
    integer lcf    
    ! f0 -- spectra and v2 (vs pt)
    character*128 filename_0	!(10)
    integer l0f 
    ! f0 -- spectra and v1 (vs y)
    character*128 filename_y0	!(12)
    integer l0v 
    !!!!
    ! polarization and Vorticity
    character*128 filename_pol		!(56)
    integer lif_pol     
    ! polarization and Vorticity - boosted
    character*128 fn_polboost		!(57)
    integer lif_pb    
    ! polarization and Vorticity, f=f0 
    character*128 filename_pol_0		!(58)
    integer lif_pol_0     
    ! polarization and Vorticity - boosted, f=f0 
    character*128 fn_polboost_0		!(59)
    integer lif_pb_0
    !!---------
    
    
   
    real, dimension (npart, nrap, npt) :: v1, v2, integrated_spectra
    real, dimension (npart, nrap) :: y_spectra, v1_y, v2_y
    real value_rap, value


! 	------- debugging purpose
    ! components
    real value_x, value_y, value_tau, value_eta
    real, dimension (npart, nrap, npt) :: integrated_spectra_x, integrated_spectra_y, integrated_spectra_tau, integrated_spectra_eta
    ! components f=f0 when df is present
    real value_0, value_0_x, value_0_y, value_0_tau, value_0_eta
    real, dimension (npart, nrap, npt) :: int_s0_x,int_s0_y, int_s0_tau, int_s0_eta
    real, dimension (npart, nrap, npt) :: v2_0, v1_0, int_s0 
    real, dimension (npart, nrap) :: y_spectra_0, v1_0_y,v2_0_y
! 	---------
    
    
    v1=0.0
    v2=0.0
    integrated_spectra=0.0
    y_spectra=0.0
    v1_y=0.0
    

! 	------- debugging purpose
    integrated_spectra_x=0.0
    integrated_spectra_y =0.0
    integrated_spectra_tau=0.0
    integrated_spectra_eta=0.0
    
    int_s0=0.0
    v2_0=0.0
    v1_0=0.0
    v1_0_y=0.0
    v2_0_y=0.0
    y_spectra_0=0.0
    int_s0_x=0.0
    int_s0_y=0.0
    int_s0_tau=0.0
    int_s0_eta=0.0
! 	------- 


    lid=index(outdir,' ')-1
    print*,'***************************************************'
    print *, "generating the spectra file, you can find them in the path:"

    do ipart=1, npart	
      particle=''
      particle=name(ipart)
      lip=index(particle, ' ')-1
    
      stringa32= 'phi_spectra_' 
      call io_build_particle_filename(ipart, stringa32, filename, lif)
      open(unit=22,status='replace',file=filename(1:lif),iostat=fstat)
      call check_file(fstat, filename(1:lif))    
      print *, filename(1:lif)
	
      stringa32= 'spectra_v2_' 
      call io_build_particle_filename(ipart, stringa32, filename_int, lif_i)   
      open(unit=26,status='replace',file=filename_int(1:lif_i),iostat=fstat)
      call check_file(fstat, filename_int(1:lif_i))    
      print *, filename_int(1:lif_i)
      
      stringa32= 'rap_spectra_' 
      call io_build_particle_filename(ipart, stringa32, filename_ys, lif_ys)   
      open(unit=21,status='replace',file=filename_ys(1:lif_ys),iostat=fstat)
      call check_file(fstat, filename_ys(1:lif_ys))
      print *, filename_ys(1:lif_ys)
      

      ! 	------- debugging purpose      
      stringa32='components_' 
      print *, stringa32
      call io_build_particle_filename(ipart, stringa32, filename_components, lcf)          
      open(unit=11,status='replace',file=filename_components(1:lcf),iostat=fstat)
      call check_file(fstat, filename_components(1:lcf))      
      
      if (viscosity_spectra == 1) then
	stringa32='f0_dN_v2_'
	call io_build_particle_filename(ipart, stringa32, filename_0, l0f) 
	open(unit=10,status='replace',file=filename_0(1:l0f),iostat=fstat)
	call check_file(fstat, filename_0(1:l0f))
	stringa32='f0_dNdy_v1_'
	call io_build_particle_filename(ipart, stringa32, filename_y0, l0v) 
	open(unit=12,status='replace',file=filename_y0(1:l0v),iostat=fstat)
	call check_file(fstat, filename_y0(1:l0v))
      endif
!           	------- ------- 
      if (vorticity_flag == 1) then
	stringa32='polar_'
	call io_build_particle_filename(ipart, stringa32, filename_pol, lif_pol) 
	open(unit=56,status='replace',file=filename_pol(1:lif_pol),iostat=fstat)
	call check_file(fstat, filename_pol(1:lif_pol))    
	write (56,*) "## This file contains the polarization of the particle ", particle(1:lip)
	write (56,*) "## "
	write (56,*) m(ipart), g(ipart), nrap, npt, nphi
	write (56,*) rapidity
	write (56,*) pt
	write (56,*) phi
	stringa32='polar_boost_'
	call io_build_particle_filename(ipart, stringa32, fn_polboost, lif_pb) 
	open(unit=57,status='replace',file=fn_polboost(1:lif_pb),iostat=fstat)
	call check_file(fstat, fn_polboost(1:lif_pol))    
	write (57,*) "## This file contains the polarization of the particle ", particle(1:lip)
	write (57,*) "## "
	write (57,*) m(ipart), g(ipart), nrap, npt, nphi
	write (57,*) rapidity
	write (57,*) pt
	write (57,*) phi

	if (viscosity_spectra == 1) then
	  stringa32='f0_polar_'
	  call io_build_particle_filename(ipart, stringa32, filename_pol_0, lif_pol_0) 
          open(unit=58,status='replace',file=filename_pol_0(1:lif_pol_0),iostat=fstat)
	  call check_file(fstat, filename_pol_0(1:lif_pol_0))    
	  write (58,*) "## This file contains the polarization of the particle ", particle(1:lip), &
	  & "for its distribution function NOT corrected with viscosity"
	  write (58,*) "## "
	  write (58,*) m(ipart), g(ipart), nrap, npt, nphi
	  write (58,*) rapidity
	  write (58,*) pt
	  write (58,*) phi
	  stringa32='f0_polar_boost_'
	  call io_build_particle_filename(ipart, stringa32, fn_polboost_0, lif_pb_0) 
	  open(unit=59,status='replace',file=fn_polboost_0(1:lif_pb_0),iostat=fstat)
	  call check_file(fstat, fn_polboost_0(1:lif_pol_0))    
	  write (59,*) "## This file contains the polarization of the particle ", particle(1:lip), &
	  & "for its distribution function NOT corrected with viscosity"
	  write (59,*) "## "
	  write (59,*) m(ipart), g(ipart), nrap, npt, nphi
	  write (59,*) rapidity
	  write (59,*) pt
	  write (59,*) phi
	endif 
      endif
      
      write (21,*) "## This spectra is integrated over phi and pt"
      write (21,*) "## columns: 1-particle rapidity  2-spectra 3-v1 4-v2"
! 
      write (22,*) "## Thermal spectra of", particle(1:lip), 'mass =', m(ipart), ", depends of phi and particle rapidity"
      write (22,*) "## column (1)-pt (2n)-phi (2n+1)-spectra(pt, phi, y),  last-rapidity"

      write (26,*) "## Thermal spectra and elliptic flow of", particle(1:lip), 'mass =', m(ipart) 
      write (26,*) "## column 1-pt 2-integrated spectra 3 -averaged (over phi) spectra 4 - direct flow &
      & 5 - elliptic flow 6 -particle rapidity " 
      
      ind_rap= nrap-nrap/2
      write (26,*) "## INDEX for zero rapidity ", ind_rap, "check y=", rapidity(ind_rap)
           
      do irap=1, nrap
	do ipt=1,npt  
	  do iphi=1, nphi
	      if (vorticity_flag == 1) then
		! polarization in the lab
		write(56,*) (pola(ivort, ipart, irap, ipt, iphi)/spectra(ipart, irap, ipt, iphi), ivort=0,3)  
		!polarization in the cdm
		do ivort=0,3
		  pi0(ivort)=(pola(ivort, ipart, irap, ipt, iphi)-pola_boost(ivort, ipart, irap, ipt, iphi))/spectra(ipart, irap, ipt, iphi)
		end do		
		pi0(0)=pi0(0)/m(ipart)
		write(57,*) (pi0(ivort), ivort=0,3)
		!---------------------
		if (viscosity_spectra ==1 ) then
		  ! polarization in the lab f=f0
		  write(58,*) (pola_0(ivort, ipart, irap, ipt, iphi)/spectra_0(ipart, irap, ipt, iphi), ivort=0,3)  
		  !polarization in the cdm
		  do ivort=0,3
		    pi0_0(ivort)=(pola_0(ivort, ipart, irap, ipt, iphi)- &
		    & pola_boost_0(ivort, ipart, irap, ipt, iphi))&
		    & /spectra_0(ipart, irap, ipt, iphi)
		  end do		
		  pi0_0(0)=pi0_0(0)/m(ipart)
		  write(59,*) (pi0_0(ivort), ivort=0,3)
		endif 
	      endif
	  
	  
	  
	      value=g(ipart)*spectra(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      spectra(ipart, irap, ipt, iphi)=value
	      integrated_spectra(ipart, irap, ipt)=integrated_spectra(ipart, irap, ipt)+(dphi*value)
	      v1(ipart, irap, ipt)=v1(ipart, irap, ipt)+(value*dphi*cosphi(iphi))
	      v2(ipart, irap, ipt)=v2(ipart, irap, ipt)+(value*dphi*cos(2.0*phi(iphi)))
	     
! ! 	      ------- debugging purpose --- components
	      value_x=g(ipart)*spectra_x(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      value_y=g(ipart)*spectra_y(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      value_eta=g(ipart)*spectra_eta(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      value_tau=g(ipart)*spectra_tau(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      integrated_spectra_x(ipart, irap, ipt)=integrated_spectra_x(ipart, irap, ipt)+(dphi*value_x)
	      integrated_spectra_y(ipart, irap, ipt)=integrated_spectra_y(ipart, irap, ipt)+(dphi*value_y)
	      integrated_spectra_eta(ipart, irap, ipt)=integrated_spectra_eta(ipart, irap, ipt)+(dphi*value_eta)
	      integrated_spectra_tau(ipart, irap, ipt)=integrated_spectra_tau(ipart, irap, ipt)+(dphi*value_tau)
! ! 	      -------  --- components

! ! 	      ------- debugging purpose --- non corrected distribution function
	      value_0=g(ipart)*spectra_0(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      value_0_tau=g(ipart)*s_0_tau(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      value_0_x=g(ipart)*s_0_x(ipart, irap, ipt, iphi)/TWOPIhbarc_cube
	      value_0_y=g(ipart)*s_0_y(ipart, irap, ipt, iphi)/TWOPIhbarc_cube	      
	      value_0_eta=g(ipart)*s_0_eta(ipart, irap, ipt, iphi)/TWOPIhbarc_cube	
	      
	      int_s0(ipart, irap, ipt)=int_s0(ipart, irap, ipt)+(dphi*value_0)
	      int_s0_tau(ipart, irap, ipt)=int_s0_tau(ipart, irap, ipt)+(dphi*value_0_tau)
	      int_s0_x(ipart, irap, ipt)=int_s0_x(ipart, irap, ipt)+(dphi*value_0_x)
	      int_s0_y(ipart, irap, ipt)=int_s0_y(ipart, irap, ipt)+(dphi*value_0_y)
	      int_s0_eta(ipart, irap, ipt)=int_s0_eta(ipart, irap, ipt)+(dphi*value_0_eta)
	      
	      v1_0(ipart, irap, ipt)=v1_0(ipart, irap, ipt)+(value_0*dphi*cos(phi(iphi)))
	      v2_0(ipart, irap, ipt)=v2_0(ipart, irap, ipt)+(value_0*dphi*cos(2.0*phi(iphi)))      
! ! 	      -------  
	  end do ! phi   
	  v1_y(ipart,irap)=v1_y(ipart,irap)+pt(ipt)*v1(ipart, irap, ipt)*dpt
	  v2_y(ipart,irap)=v2_y(ipart,irap)+pt(ipt)*v2(ipart, irap, ipt)*dpt
	  y_spectra(ipart, irap)=y_spectra(ipart, irap)+pt(ipt)*integrated_spectra(ipart, irap, ipt)*dpt
	  
	  v2(ipart, irap, ipt)= v2(ipart, irap, ipt)/integrated_spectra(ipart, irap, ipt)
	  v1(ipart, irap, ipt)= v1(ipart, irap, ipt)/integrated_spectra(ipart, irap, ipt)
	  ! --- f=f0
	  v1_0_y(ipart,irap)=v1_0_y(ipart,irap)+pt(ipt)*v1_0(ipart, irap, ipt)*dpt
	  v2_0_y(ipart,irap)=v2_0_y(ipart,irap)+pt(ipt)*v2_0(ipart, irap, ipt)*dpt
	  y_spectra_0(ipart, irap)=y_spectra_0(ipart, irap)+pt(ipt)*int_s0(ipart, irap, ipt)*dpt
	  
	  v2_0(ipart, irap, ipt)= v2_0(ipart, irap, ipt)/int_s0(ipart, irap, ipt)
	  v1_0(ipart, irap, ipt)= v1_0(ipart, irap, ipt)/int_s0(ipart, irap, ipt)
	  !----------
	  
	  write (22,*)  pt(ipt), (phi(iphi), spectra(ipart, irap, ipt, iphi),  iphi=1, nphi), rapidity(irap) 

	  write (26,*) pt(ipt), integrated_spectra(ipart, irap, ipt), &
	  & integrated_spectra(ipart, irap, ipt)/(phimax-phimin), &
	  & v1(ipart, irap, ipt), v2(ipart, irap, ipt),  rapidity(irap)	 	  
	  
	  write(11,*)  pt(ipt), integrated_spectra_tau(ipart, irap, ipt),&
	  & integrated_spectra_x(ipart, irap, ipt), integrated_spectra_y(ipart, irap, ipt),&
	  & integrated_spectra_eta(ipart, irap, ipt), rapidity(irap)	  
	  
	  if (viscosity_spectra == 1) then
	    write (10,*) pt(ipt),  int_s0(ipart, irap, ipt), &
	    & int_s0_tau(ipart, irap, ipt), int_s0_x(ipart, irap, ipt),  &
	    & int_s0_y(ipart, irap, ipt),  int_s0_eta(ipart, irap, ipt), &  
	    & v1_0(ipart, irap, ipt), v2_0(ipart, irap, ipt), &
	    & integrated_spectra(ipart, irap, ipt)-int_s0(ipart, irap, ipt), & 
	    & rapidity(irap)
	  endif   
! ! 	      ------- 
	end do !pt

	write (21,*)  rapidity(irap), y_spectra(ipart, irap), v1_y(ipart,irap),&
	& v2_y(ipart,irap), v1_y(ipart,irap)/y_spectra(ipart, irap),&
	& v2_y(ipart,irap)/y_spectra(ipart, irap)
	if (viscosity_spectra == 1) then
	  write (12,*)  rapidity(irap), y_spectra_0(ipart, irap), &
	  &v1_0_y(ipart,irap), v2_0_y(ipart,irap), &
	  &v1_0_y(ipart,irap)/y_spectra_0(ipart,irap), v2_0_y(ipart,irap)/y_spectra_0(ipart, irap)	
	endif 	

	write (22,*) ' '
	write (10,*) ' '
	write (10,*) ' '
	write (11,*) ' '
	write (11,*) ' '
	write (26,*) ' '
	write (26,*) ' '
      end do ! rapidity	         
    close(10)
    close(11)
    close(12)
    close(21)
    close(22)
    close(26)      
    close(56)
    close(57)
    close(58)
    close(59)
  end do ! npart
  return
  end subroutine io_print_spectra


! !******************************************************************
subroutine io_end_setup(rep)
 use common
 use common_MC
 use common_thermal
 implicit none 
 
 integer, intent (in) :: rep

 integer filerror, i
 
 2003 FORMAT(I6,I10,2X, A25,F10.6,I5,I3,3x, F10.6)
 108	format(a30,3x,i3)
 109	format(a30,3x,a22)
 107	format(a30,f10.3,a3,f10.3)
 
 select case (rep)
  case(0) !!THERMAL
    open(unit=21,status='old',position='append', file=outdir(1:LID_out)//'report_thermal.txt', form='formatted', iostat=filerror)
    call check_file(filerror, outdir(1:LID_out)//'report_thermal.txt')
    write (21, *) '____________________________________________________________________'
    write (21, *) ' '
    !write (21, *) "Energy (GeV)", Energy_int
    write (21, *) "Frozen cells", frozen_cells
  case(1) !!PAGE
   open(unit=21,status='old',position='append',file=outdir(1:LID_out)//'report_PaGe.txt', form='formatted', iostat=filerror)
   call check_file(filerror, outdir(1:LID_out)//'report_PaGe.txt')
    write (21, *) '____________________________________________________________________'
    write (21, *) ' '
    !write (21, *) "Energy (GeV)", Energy_int
    write (21, *) "Frozen cells", frozen_cells
    write (21, *) "Produced Particles", produced_particles
  case(2) !!PAMOGE
   open(unit=21,status='old',position='append',file=outdir(1:LID_out)//'report_PaMoGe.txt', form='formatted', iostat=filerror)   
   call check_file(filerror, outdir(1:LID_out)//'report_PaMoGe.txt')
    write (21, *) '____________________________________________________________________'
    write (21, *) ' '
    !write (21, *) "Energy (GeV)", Energy_int
    write (21, *) "Frozen cells", frozen_cells
    write (21, *) "Produced Particles", produced_particles
    write (21, *) "Particle Energy", particle_energy
  case(3) !!MERGE
   open(unit=21,status='old',position='append',file=outdir(1:LID_out)//'report_merge.txt', form='formatted', iostat=filerror)   
   call check_file(filerror, outdir(1:LID_out)//'report_merge.txt')
    write (21, *) '____________________________________________________________________'
    write (21, *) ' '
    !write (21, *) "Energy (GeV)", Energy_int
    write (21, *) "Frozen cells", frozen_cells
    write (21, *) "Produced Particles", produced_particles
    write (21, *) "Particle Energy", particle_energy
 end select
 close(21)
 return      
end subroutine io_end_setup  
! !******************************************************************
  subroutine check_file_folder_existence (oou, LID_o, iin, fffi, LID_i)
    implicit none
      integer, intent(in):: LID_i, LID_o
      character*32, intent(in):: oou, iin, fffi
      
      integer filerror
      
      open(unit=21,status='old',file=iin(1:LID_i)//'hypersurface.txt', form='formatted', iostat=filerror )
      call check_file(filerror, iin(1:LID_i)//'hypersurface.txt')
      close(21)      
  
      !I create the output directory, if it does not exists already
      call EXECUTE_COMMAND_LINE('mkdir -p '//oou(1:LID_o)//'')
      

    return
  end subroutine check_file_folder_existence
end module io
