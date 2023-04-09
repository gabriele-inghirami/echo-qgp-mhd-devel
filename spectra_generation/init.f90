! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *         
! *                                                                           *         
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *         
! *  File: analysis/init.f90                                                  *
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
! *  Acknowledgments: Alessandro Drago                                        *
! *                                                                           *
! *****************************************************************************


  module init 
    use common, only: npart, m, viscosity_echo, viscosity_spectra
    implicit none 
    
    contains
! !******************************************************************  
    subroutine allocate_momenta()	
      use common_thermal
      implicit none 
      integer AllocateStatus, iphi, ipt, irap, ipart
      integer check
      character*1 switch
      
    	  ALLOCATE ( pt(npt), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  ALLOCATE ( rapidity(nrap), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0) 	STOP "*** Not enough memory ***"
	  ALLOCATE ( phi(nphi), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  ALLOCATE ( cosphi(nphi), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  ALLOCATE ( senphi(nphi), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"

	  dpt=(ptmax-ptmin)/(1.0*npt)
	  dphi=(phimax-phimin)/(1.0*nphi)
	  
	  if (nrap==1) then
	    rapidity(1)=0.0
	    check=1
	    drap = 0.0
	  else 
	    drap=(rapmax-rapmin)/(1.0*(nrap-1))
	    check=0
	    do irap=1, nrap
	      rapidity(irap)=rapmin+(irap-1)*drap	
	      if (rapidity(irap)==0.0) then 
		check=1
	      endif
	    end do
	  endif


	  if (check==0) then
	    print *, "Warning: YOU HAVE CHOSEN AN INTERVAL FOR PARTICLE RAPIDITY THAT &
	     & DOES NOT INCLUDE THE VALUE y=0.0"
	      print *, rapidity
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
		  case default 
		    call exit
		  end select
	  endif
	  
	  do iphi=1, nphi
	    phi(iphi)=phimin+(iphi-1)*dphi	    
	    senphi(iphi)=sin(phi(iphi))
	    cosphi(iphi)=cos(phi(iphi))
	  end do
	  
	  if (viscosity_echo .eq. 1) then
	    if (viscosity_spectra .eq. 1) then 
	      ALLOCATE ( cosphi2(nphi), STAT = AllocateStatus)
	      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	      ALLOCATE ( senphi2(nphi), STAT = AllocateStatus)
	      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	      ALLOCATE ( sencosphi(nphi), STAT = AllocateStatus)
	      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	      do iphi=1, nphi
		senphi2(iphi)=senphi(iphi)*senphi(iphi)
		cosphi2(iphi)=cosphi(iphi)*cosphi(iphi)
		sencosphi(iphi)=senphi(iphi)*cosphi(iphi)
	      end do 
	    endif
	 endif

      return
    end subroutine allocate_momenta  
! ! !******************************************************************
    subroutine init_allocate_spectra()	
      use common
      use common_thermal
      implicit none 
      integer AllocateStatus, iphi, ipt, irap, ipart
      
      Energy_int=0.0

      ALLOCATE ( spectra(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      spectra=0.0
!       ----- debugging purpose
      ALLOCATE ( spectra_x(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( spectra_y(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( spectra_tau(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( spectra_eta(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      spectra_tau=0.0
      spectra_x=0.0
      spectra_y=0.0
      spectra_eta=0.0

      ALLOCATE ( s_0_x(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( s_0_y(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( s_0_tau(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( s_0_eta(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
      ALLOCATE ( spectra_0(npart,nrap,npt,nphi), STAT = AllocateStatus)
      IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"	    
      s_0_tau=0.0
      s_0_x=0.0
      s_0_y=0.0
      s_0_eta=0.0	    
      spectra_0=0.0
	    
      if (vorticity_flag==1) then
	ALLOCATE ( pola(0:3,npart,nrap,npt,nphi), STAT = AllocateStatus)
	IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	ALLOCATE ( pola_boost(0:3,npart,nrap,npt,nphi), STAT = AllocateStatus)
	IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	ALLOCATE ( u_der(0:3,0:3), STAT = AllocateStatus)	
	IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"     
	ALLOCATE ( T_der(0:3), STAT = AllocateStatus)	
	IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"  
	ALLOCATE ( exter(1:6), STAT = AllocateStatus)	
	IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"  
	
	pola=0.0
	pola_boost=0.0
	u_der=0.0
	T_der=0.0
	exter=0.0
	if (viscosity_spectra ==1) then 
	  ALLOCATE ( pola_0(0:3,npart,nrap,npt,nphi), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  ALLOCATE ( pola_boost_0(0:3,npart,nrap,npt,nphi), STAT = AllocateStatus)
	  IF (AllocateStatus /= 0)   	STOP "*** Not enough memory ***"
	  pola_0=0.0
	  pola_boost_0=0.0
	  
	endif 
      endif
      return
    end subroutine init_allocate_spectra    
! !******************************************************************  
    end module init
    

 
 

    
