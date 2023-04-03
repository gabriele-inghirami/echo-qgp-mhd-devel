! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 *         
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015,2016 The ECHO-QGP team                                * 
! *                                                                           *
! *  File: eos.f90                                                            *
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
! *           Luca Del Zanna (luca.delzanna@unifi.it)                         *         
! *                                                                           *
! *  Contributors:                                                            *
! *                                                                           *
! *  Acknowledgments: Andrea Beraudo, Francesco Becattini,                    *
! *                   Vinod Chandra, Arturo De Pace,                          *
! *                   Valentina Rolando                                       *
! *                                                                           *
! *****************************************************************************

module eos
      use parallel, only: pe0
      implicit none
      real(8) :: pressure, energy_density, charge_density
      real(8),allocatable,dimension(:) :: temperature_array, energy_density_array, pressure_array, sound_speed_array 
      integer :: entries !number of tabulated data
      integer :: numerically=1 !if set to 1 inside it enables numerical computation of partial derivatives of the analytic eos
      procedure(), pointer :: eos_energy, eos_sound, get_derived_data, eos_pressure, eos_pressure0, get_derived_temp_for_check,&
                             &eos_energy_with_ders
      real(8) :: num_deriv_increment=1.d-10 !increment of numerical derivative of pr(rh,en)
      real(8) :: s_acc=1.d-10 ! accuracy when computing entropy density
      real(8) :: root_accuracy=1.d-10 ! accuracy of the root finding subroutine
      real(8) :: min_en_dens=0.0d0, max_en_dens=2.d4 ! minimum and maximum possible values for energy density in GeV/fm^3
      integer :: temp_check=-1 ! flag to stop the simulation if temperature drops down below a chosen value
      logical :: called_by_tools=.false. !set to true only when this module is called not by echo-qgp, but from postproc tools
      integer :: eqstate !it determines if we wish to use an analytic equation of state (and which one) or a tabulated one
      real(8) :: przero, enezero=-1 ! minimum values for press and energy dens. at the initalization stage (no effects with tab eos)
      character(60) :: eos_tab_name
      logical,parameter :: check_on_pressure=.true. !enable or disable the check on pressure
      real(8),parameter :: gamma_adindex=4.d0/3.d0 !gamma adiabatic index
      real(8) :: param1_defined_in_common !to be initialized with other variables defined in common.f90, which is compiled later
      !real(8),parameter :: gamma_adindex=2. !gamma adiabatic index

contains

! *****************************************************************************

subroutine find_enezero(pr_input, en_output,errcode)
      implicit none
      real(8), intent(in) :: pr_input
      real(8), intent(out) :: en_output
      real(8) :: prguess, prtrya, enguess=1.d-15
      real(8) :: tol=1.d-15
      integer :: n, nmax=100
      real(8) :: a,b,c,g,en,temp,en_min_bag_eos=0.d0
      integer, intent(out) :: errcode

      errcode=0

      !analityc eos when using a temperature.def file
      if(eqstate .lt. 3) then
        INCLUDE 'eos_data/temperature.def'
        en_min_bag_eos=1.1*b
      end if 

      if(eqstate .eq. 3) then
        call get_energy(entries, en_output, energy_density_array, pressure_array, pr_input,errcode)
      else
       call eos_pressure0(0,enguess,prguess,errcode)
       do while ( prguess .lt. pr_input)
          enguess=enguess*10.
          call eos_pressure0(0,enguess,prguess,errcode)
       end do
      
       !bisection method

       b=enguess
       a=enguess/10.
       
       do n=1,nmax
          c=(a+b)/2.
          call eos_pressure0(0,c,prguess,errcode)
          if(( prguess .eq. pr_input) .or. ( (b-a)/2. .lt. tol )) exit
          call eos_pressure0(0,a,prtrya,errcode)
          if(sign(real(1.,8),prguess-pr_input) .eq. sign(real(1.,8),prtrya-pr_input)) then
             a=c
          else
             b=c
          endif
       end do

       if(n .eq. nmax) then
         write(*,*) "Sorry, but I have to quit, because I'm not able to find the energy value corresponding to",pr_input
         errcode=51
         return
       end if
      
       en_output=c
      end if
      
      if(en_output .lt. en_min_bag_eos) then
         write(*,*) 'The bag constant inside the file defining the temperature relation imposes a minum value for energy density'
         call eos_pressure0(0,en_min_bag_eos,przero,errcode)
         en_output=en_min_bag_eos
      end if
      
      if(pe0) write(*,*) 'The minimum value for energy density (GeV/fm^3) allowed before adjusting it is:',en_output

end subroutine find_enezero

! *****************************************************************************

subroutine find_przero(pr_output, en_input,errcode)
      implicit none
      real(8), intent(out) :: pr_output
      real(8), intent(in) :: en_input
      real(8) :: rh_bitbucket, dprdrh_bitbucket, dprden_bitbucket
      integer, intent(out) :: errcode
      
      rh_bitbucket=0.
      errcode=0

      if(check_on_pressure) then
    
        call eos_pressure(rh_bitbucket,en_input,pr_output,dprdrh_bitbucket,dprden_bitbucket,errcode)

        if(pr_output .le. 0.) then
           if(pe0) write(*,*) 'The minimum value for pressure corresponding to enezero is negative, echo-qgp will not start!!!'
           call exit(1)
        end if
      else
        pr_output=0.
      end if
      

end subroutine find_przero

! *****************************************************************************

subroutine set_eos_pointers(eqstate, eos_energy, eos_pressure, eos_sound)
  implicit none
  integer :: eqstate
  procedure(), pointer :: eos_pressure, eos_energy, eos_sound
  if(eqstate .lt. 3) then
    eos_energy=>eos_energy_ana
    eos_energy_with_ders=>eos_energy_with_ders_ana
    eos_pressure=>eos_pressure_ana
    eos_pressure0=>eos_pressure_ana0
    eos_sound=>eos_sound_ana
    get_derived_data=>get_derived_data_ana
    get_derived_temp_for_check=>get_derived_temp_for_check_ana
  else if (eqstate .eq. 3) then
    eos_energy=>eos_energy_tab
    eos_energy_with_ders=>eos_energy_with_ders_tab
    eos_pressure=>eos_pressure_tab
    eos_pressure0=>eos_pressure_tab0
    eos_sound=>eos_sound_tab
    get_derived_data=>get_derived_data_tab
    get_derived_temp_for_check=>get_derived_temp_for_check_tab
  else if (eqstate .eq. 4) then
    call set_tab_from_an_eos()
    eos_energy=>eos_energy_tab
    eos_energy_with_ders=>eos_energy_with_ders_tab
    eos_pressure=>eos_pressure_tab
    eos_pressure0=>eos_pressure_tab0
    eos_sound=>eos_sound_tab
    get_derived_data=>get_derived_data_tab
    get_derived_temp_for_check=>get_derived_temp_for_check_tab
  end if
end subroutine set_eos_pointers

! *****************************************************************************

subroutine eos_pressure_tab(rh,en,pr,dprdrh,dprden,errcode)

  real(8),intent( in) :: rh,en
  real(8),intent(out) :: pr,dprdrh,dprden
  integer, intent(out) :: errcode

  errcode=0
  
     call get_pressure(entries, en, energy_density_array, pressure_array, temperature_array, pr, dprdrh, dprden,errcode)

end subroutine eos_pressure_tab

! *****************************************************************************

subroutine eos_pressure_tab0(rh_input, en_input, pr_output, errcode)

  real(8),intent(out) :: pr_output
  real(8),intent(in) :: en_input
  real(8) :: rh_input
  integer, intent(out) :: errcode

  errcode=0
  rh_input=0.
  
     call get_pressure0(entries, en_input, energy_density_array, pressure_array, pr_output,errcode)

end subroutine eos_pressure_tab0

! *****************************************************************************

subroutine eos_energy_tab(rh,en,pr,errcode)

  real(8),intent( in) :: rh,pr
  real(8),intent(out) :: en
  integer, intent(out) :: errcode

  errcode=0

     call get_energy(entries, en, energy_density_array, pressure_array, pr,errcode)

end subroutine eos_energy_tab

! *****************************************************************************

subroutine eos_energy_with_ders_tab(rh,en,pr,dendrh,dendpr,errcode)

  real(8),intent( in) :: rh,pr
  real(8),intent(out) :: en,dendrh,dendpr
  integer, intent(out) :: errcode

  errcode=0

     call get_energy_with_ders(entries, en, energy_density_array, pressure_array, pr, dendrh, dendpr, errcode)

end subroutine eos_energy_with_ders_tab

! *****************************************************************************

subroutine eos_sound_tab(rh,en,pr,cs2,errcode)

  real(8),intent( in) :: rh,en,pr
  real(8),intent(out) :: cs2
  integer, intent(out) :: errcode
  errcode=0
  
  call get_sound(entries,en,energy_density_array,cs2,sound_speed_array,errcode)
  
end subroutine eos_sound_tab

! *****************************************************************************

!PLEASE, NOTE THAT THIS MODULE WORKS ONLY IF PRESSURE AND ENERGY ARE ALWAYS GROWING WITH TEMPERATURE

      subroutine read_data_from_tabulated_file(entries, eos_tab_name, temperature_array, energy_density_array, pressure_array,&
      & sound_speed_array, min_en_dens, max_en_dens)
      implicit none
      real(8),allocatable,dimension(:) :: temperature_array, energy_density_array, pressure_array, sound_speed_array
      character(60), intent(in) :: eos_tab_name
      real(8) :: min_en_dens, max_en_dens ! minimum and maximum values of energy density
      integer :: entries !number of tabulated data
      integer :: filerror !to check errors when opening data file
      integer :: arrayerror !to check errors when allocating arrays
      integer :: i !simply a counter
      real(8), parameter :: hc3=.00768351d0
      
      if(called_by_tools) then
         open(1,file='../eos_data/'//eos_tab_name,status='OLD',iostat=filerror)
      else
         open(1,file='eos_data/'//eos_tab_name,status='OLD',iostat=filerror)
      end if
      if (filerror .ne. 0) then
      write(*,*) 'eos_data/'//eos_tab_name, " cannot be opened, so I'm forced to quit!"
      close(1)
      call exit(1)
      end if
      
      read(1,*) entries
      
      !allocate array containing temperature data
      allocate(temperature_array(1:entries),STAT=arrayerror)
      if(arrayerror .ne. 0) then
      write(*,*) "Can't allocate memory for the array containing temperature data... Sorry, I'm forced to quit!"
      close(1)
      call exit(1)
      end if
      
      !allocate array containing energy density data
      allocate(energy_density_array(1:entries),STAT=arrayerror)
      if(arrayerror .ne. 0) then
      write(*,*) "Can't allocate memory for the array containing energy density data... Sorry, I'm forced to quit!"
      close(1)
      call exit(1)
      end if
      
      !allocate array containing pressure data
      allocate(pressure_array(1:entries),STAT=arrayerror)
      if(arrayerror .ne. 0) then
      write(*,*) "Can't allocate memory for the array containing pressure data... Sorry, I'm forced to quit!"
      close(1)
      call exit(1)
      end if
      
      !allocate array containing sound speed data
      allocate(sound_speed_array(1:entries),STAT=arrayerror)
      if(arrayerror .ne. 0)then
      write(*,*) "Can't allocate memory for the array containing sound speed data... Sorry, I'm forced to quit!"
      close(1)
      call exit(1)
      end if
      
      !read data from file
      do i=1,entries
      read(1,*) temperature_array(i), energy_density_array(i), pressure_array(i), sound_speed_array(i)
      energy_density_array(i)=energy_density_array(i)/hc3*temperature_array(i)**4.0
      pressure_array(i)=pressure_array(i)/hc3*temperature_array(i)**4.0
      end do
      
      min_en_dens=energy_density_array(1)
      max_en_dens=energy_density_array(entries)
      
      close(1)     
      
      end subroutine
      
! *****************************************************************************

      subroutine get_pressure0(entries, energy_value, energy_density_array, pressure_array, pressure_value,errcode)
      implicit none
      real(8) :: energy_value, pressure_value
      integer :: entries
      real(8),dimension(entries) :: energy_density_array, pressure_array
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: p1, p2, e1, e2
      integer, intent(out) :: errcode
      errcode=0
      
      !check if we are in the range of values where we are able to interpolate
      if(energy_value .lt. energy_density_array(1)) then
      write(*,*) 'There is a value of energy density too small to interpolate...'
      write(*,*) energy_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=121
      return
      endif
      
      if(energy_value .gt. energy_density_array(entries)) then
      write(*,*) 'There is a value of energy density too big to interpolate...'
      write(*,*) energy_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=122
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
            
      do while (indxl+1 .ne. indxh) 
      if( energy_value .gt. energy_density_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
      
      !now we do linear interpolation using the following formula
      !let e2=energy_density_array(indxh), e1=energy_density_array(indxl) and x = (energy_value - e1)/(e2-e1)
      !let p2=pressure_array(indxh) and p1=pressure_array(indxl)
      !we compute the pressure corrisponding at energy_value as p1+(p2-p1)*x
      !now let t2=temperature_array(indxh) and t1=temperature_array(indxl)
      !we compute the derivative of pressure with respect to temperature simply as (p2-p1)/(t2-t1)
      !and the the derivative of energy with respect to temperature simply as (e2-e1)/(t2-t1)
      p2=pressure_array(indxh)
      p1=pressure_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)
      pressure_value=p1+(p2-p1)*(energy_value-e1)/(e2-e1)
      

      end subroutine
      
! ****************************************************************************

      subroutine get_pressure(entries, energy_value, energy_density_array, pressure_array, temperature_array, pressure_value,&
      & deriv_press_dens, deriv_press_energy, errcode)
      implicit none
      real(8) :: energy_value, pressure_value, deriv_press_temp, deriv_energy_temp, deriv_press_dens, deriv_press_energy
      integer :: entries
      real(8),dimension(entries) :: energy_density_array, pressure_array, temperature_array
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: p1, p2, t1, t2, e1, e2
      integer, intent(out) :: errcode

      errcode=0
      
      !check if we are in the range of values where we are able to interpolate
      if(energy_value .lt. energy_density_array(1)) then
      write(*,*) 'There is a value of energy density too small to interpolate...'
      write(*,*) energy_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=123
      return
      endif
      
      if(energy_value .gt. energy_density_array(entries)) then
      write(*,*) 'There is a value of energy density too big to interpolate...'
      write(*,*) energy_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=124
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
            
      do while (indxl+1 .ne. indxh) 
      if( energy_value .gt. energy_density_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
      
      !now we do linear interpolation using the following formula
      !let e2=energy_density_array(indxh), e1=energy_density_array(indxl) and x = (energy_value - e1)/(e2-e1)
      !let p2=pressure_array(indxh) and p1=pressure_array(indxl)
      !we compute the pressure corrisponding at energy_value as p1+(p2-p1)*x
      !now let t2=temperature_array(indxh) and t1=temperature_array(indxl)
      !we compute the derivative of pressure with respect to temperature simply as (p2-p1)/(t2-t1)
      !and the the derivative of energy with respect to temperature simply as (e2-e1)/(t2-t1)
      p2=pressure_array(indxh)
      p1=pressure_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)
      t2=temperature_array(indxh)
      t1=temperature_array(indxl)
      pressure_value=p1+(p2-p1)*(energy_value-e1)/(e2-e1)
      deriv_press_temp=(p2-p1)/(t2-t1)
      deriv_energy_temp=(e2-e1)/(t2-t1)
      
      deriv_press_dens=0 !just to mantain compatibility with a more general version of this routine
      deriv_press_energy=deriv_press_temp/deriv_energy_temp

      end subroutine
      
! ****************************************************************************
      
      subroutine get_energy(entries, energy_value, energy_density_array, pressure_array, pressure_value,errcode)
      implicit none
      real(8),intent(out):: energy_value
      real(8),intent(in):: pressure_value
      integer :: entries
      real(8),dimension(entries) :: energy_density_array, pressure_array
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: p1, p2, e1, e2
      integer, intent(out) :: errcode

      errcode=0
      
      !check if we are in the range of values where we are able to interpolate
      if(pressure_value .lt. pressure_array(1)) then
      write(*,*) 'There is a value of pressure too small to interpolate...'
      write(*,*) pressure_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=125
      return
      endif
      
      if(pressure_value .gt. pressure_array(entries)) then
      write(*,*) 'There is a value of pressure too big to interpolate...'
      write(*,*) pressure_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=126
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
            
      do while (indxl+1 .ne. indxh) 
      if( pressure_value .gt. pressure_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
            
      !now we do linear interpolation using the following formula
      !let p2=pressure_array(indxh) and p1=pressure_array(indxl) and x = (pressure_value - p1)/(p2-p1)
      !let e2=energy_density_array(indxh), e1=energy_density_array(indxl)
      !we compute the energy_density corrisponding at pressure_value as e1+(e2-e1)*x
      p2=pressure_array(indxh)
      p1=pressure_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)
      energy_value=e1+(e2-e1)*(pressure_value-p1)/(p2-p1)
      
      end subroutine

! ****************************************************************************
      
      subroutine get_energy_with_ders(entries, energy_value, energy_density_array, pressure_array, pressure_value, &
      & deriv_en_dens, deriv_en_press, errcode)
      implicit none
      real(8),intent(out):: energy_value,deriv_en_dens, deriv_en_press
      real(8),intent(in):: pressure_value
      integer :: entries
      real(8),dimension(entries) :: energy_density_array, pressure_array
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: p1, p2, e1, e2
      integer, intent(out) :: errcode

      errcode=0
      
      !check if we are in the range of values where we are able to interpolate
      if(pressure_value .lt. pressure_array(1)) then
      write(*,*) 'There is a value of pressure too small to interpolate...'
      write(*,*) pressure_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=125
      return
      endif
      
      if(pressure_value .gt. pressure_array(entries)) then
      write(*,*) 'There is a value of pressure too big to interpolate...'
      write(*,*) pressure_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=126
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
            
      do while (indxl+1 .ne. indxh) 
      if( pressure_value .gt. pressure_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
            
      !now we do linear interpolation using the following formula
      !let p2=pressure_array(indxh) and p1=pressure_array(indxl) and x = (pressure_value - p1)/(p2-p1)
      !let e2=energy_density_array(indxh), e1=energy_density_array(indxl)
      !we compute the energy_density corrisponding at pressure_value as e1+(e2-e1)*x
      p2=pressure_array(indxh)
      p1=pressure_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)
      energy_value=e1+(e2-e1)*(pressure_value-p1)/(p2-p1)
      deriv_en_press=(e2-e1)/(p2-p1)

      deriv_en_press=0 !just to mantain compatibility with a more general version of this routine
      
      end subroutine

! ****************************************************************************
      
      subroutine get_sound(entries, energy_value, energy_density_array, sound_speed, sound_speed_array,errcode)
      implicit none
      real(8),intent(out):: sound_speed
      real(8),intent(in):: energy_value
      integer :: entries
      real(8),dimension(entries) :: energy_density_array, sound_speed_array
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: s1, s2, e1, e2
      integer, intent(out) :: errcode
      errcode=0
      
      !check if we are in the range of values where we are able to interpolate
      if(energy_value .lt. energy_density_array(1)) then
      write(*,*) 'There is a value of energy density too small to interpolate to find the sound speed...'
      write(*,*) energy_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=131
      return
      endif
      
      if(energy_value .gt. energy_density_array(entries)) then
      write(*,*) 'There is a value of energy density too big to interpolate to find the sound speed...'
      write(*,*) energy_value
      write(*,*) 'Sorry, but I have to quit!'
      errcode=132
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
            
      do while (indxl+1 .ne. indxh) 
      if( energy_value .gt. energy_density_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
            
      !now we do linear interpolation using the same formula as the subroutine above
      s2=sound_speed_array(indxh)
      s1=sound_speed_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)
      sound_speed=s1+(s2-s1)*(energy_value-e1)/(e2-e1)
      
      end subroutine
      
! ****************************************************************************

subroutine get_derived_data_tab(rh, pressure, energy_density, temperature, entropy_density,errcode )
      implicit none
      real(8),intent(in) :: rh, pressure
      ! actually rh is unused, it is kept only for using the same interface as in the analytical case
      real(8), intent(out) :: energy_density, temperature, entropy_density 
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: t1, t2, e1, e2, p1, p2
      real(8) :: interpolation_factor
      integer, intent(out) :: errcode

      errcode=0

      !check if we are in the range of values where we are able to interpolate
      if(pressure .lt. pressure_array(1)) then
      write(*,*) 'There is a value of pressure too small to interpolate...'
      write(*,*) pressure
      write(*,*) 'Sorry, but I have to quit!'
      errcode=135
      return
      endif
      
      if(pressure .gt. pressure_array(entries)) then
      write(*,*) 'There is a value of pressure too big to interpolate...'
      write(*,*) pressure
      write(*,*) 'Sorry, but I have to quit!'
      errcode=136
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
 
      do while (indxl+1 .ne. indxh) 
      if( pressure .gt. pressure_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
      
      p2=pressure_array(indxh)
      p1=pressure_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)      
      t2=temperature_array(indxh)
      t1=temperature_array(indxl)
      interpolation_factor=(dble(pressure)-p1)/(p2-p1)
      temperature=t1+(t2-t1)*interpolation_factor
      energy_density=e1+(e2-e1)*interpolation_factor
      entropy_density=(energy_density+pressure)/temperature
    
      end subroutine get_derived_data_tab
      
      
! *****************************************************************************

subroutine get_derived_temp_for_check_tab(rh, pressure, energy_density, temperature,errcode)
      implicit none
      real(8),intent(in) :: rh, pressure
      ! actually rh is unused, it is kept only for using the same interface as in the analytical case
      real(8), intent(out) :: energy_density, temperature
      integer :: indxl
      integer :: indxh
      integer :: indx
      real(8) :: t1, t2, e1, e2, p1, p2
      real(8) :: interpolation_factor
      integer, intent(out) :: errcode
      
      !check if we are in the range of values where we are able to interpolate
      if(pressure .lt. pressure_array(1)) then
      write(*,*) 'There is a value of pressure too small to interpolate...'
      write(*,*) pressure
      write(*,*) 'Sorry, but I have to quit!'
      errcode=137
      return
      endif
      
      if(pressure .gt. pressure_array(entries)) then
      write(*,*) 'There is a value of pressure too big to interpolate...'
      write(*,*) pressure
      write(*,*) 'Sorry, but I have to quit!'
      errcode=138
      return
      endif
      
      indxl=1
      indxh=entries
      indx=(indxl+indxh)/2
            
      do while (indxl+1 .ne. indxh) 
      if( pressure .gt. pressure_array(indx)) then
      indxl=indx
      indx=(indxl+indxh)/2
      else
      indxh=indx
      indx=(indxl+indxh)/2
      endif
      enddo
      
      p2=pressure_array(indxh)
      p1=pressure_array(indxl)
      e2=energy_density_array(indxh)
      e1=energy_density_array(indxl)      
      t2=temperature_array(indxh)
      t1=temperature_array(indxl)
      interpolation_factor=(dble(pressure)-p1)/(p2-p1)
      temperature=t1+(t2-t1)*interpolation_factor
      energy_density=e1+(e2-e1)*interpolation_factor
    
      end subroutine get_derived_temp_for_check_tab
      
      
! *****************************************************************************
! AAAA this subroutine produces worse derivatives than the simpler form when deriving respect a missing variable
! probably it better to remove this subroutine if unused

subroutine num_derivative_2D_5_stencil(xa,xb,dydxa,dydxb)
  implicit none
  real(8) :: xa, xb
  real(8), intent(out) :: dydxa, dydxb
  real(8) xa0, xb0, y_xa1_p, y_xa1_m, y_xb1_p, y_xb1_m, y_xa2_p, y_xa2_m, y_xb2_p, y_xb2_m, y, y0, h

  xa0=xa
  xb0=xb
  xb=xb0*(1.0+num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xb1_p=y
  xb=xb0*(1.0-num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xb1_m=y
  xb=xb0*(1.0+2.0*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xb2_p=y
  xb=xb0*(1.0-2.0*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xb2_m=y
  xb=xb0
  xa=xa0*(1.0+num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xa1_p=y
  xa=xa0*(1.0-num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xa1_m=y
  xa=xa0*(1.0+2.0*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xa2_p=y
  xa=xa0*(1.0-2.0*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xa2_m=y
  xb=xb0
  xa=xa0
  h=xb0*num_deriv_increment
  dydxb=(-y_xb2_p+8.0*y_xb1_p-8.0*y_xb1_m+y_xb2_m)/(12.0*h)
  h=xa0*num_deriv_increment
  dydxa=(-y_xa2_p+8.0*y_xa1_p-8.0*y_xa1_m+y_xa2_m)/(12.0*h)
 

end subroutine num_derivative_2D_5_stencil

! *****************************************************************************

subroutine num_derivative_2D(xa,xb,dydxa,dydxb)
  implicit none
  real(8) :: xa, xb
  real(8), intent(out) :: dydxa, dydxb
  real(8) xa0, xb0, y_xa1_p, y_xa1_m, y_xb1_p, y_xb1_m, y, y0, h

  xa0=xa
  xb0=xb
  xb=xb0*(1.0+0.5*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xb1_p=y
  xb=xb0*(1.0-0.5*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xb1_m=y
  xb=xb0
  xa=xa0*(1.0+0.5*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xa1_p=y
  xa=xa0*(1.0-0.5*num_deriv_increment)
  call get_eos_userdef(y,xa,xb)
  y_xa1_m=y
  xb=xb0
  xa=xa0
  h=xb0*num_deriv_increment
  if ((y_xb1_p .eq. y_xb1_m) .and. ( h .eq. 0 )) then
  dydxb=0
  else
  dydxb=(y_xb1_p-y_xb1_m)/h
  endif
  h=xa0*num_deriv_increment
  if ((y_xa1_p .eq. y_xa1_m) .and. ( h .eq. 0 )) then
  dydxa=0
  else
  dydxa=(y_xa1_p-y_xa1_m)/h
  endif 

end subroutine num_derivative_2D

! *****************************************************************************

subroutine get_eos_userdef(y,xa,xb)
  implicit none
  real(8) y, xa, xb
  real(8) pr, rh, en
  real(8) a, b, c, d, e
     rh=xa
     en=xb
  if (eqstate .eq. 1) then
        INCLUDE 'eos_data/pressure_vs_rh_en1.def'
        y=pr
  else 
        INCLUDE 'eos_data/pressure_vs_rh_en2.def'
        y=pr
  endif
  
end subroutine get_eos_userdef

! ****************************************************************************

subroutine eos_pressure_ana0(xa,xb,y,errcode)
  implicit none

  real(8),intent( in) :: xa,xb
  real(8),intent(out) :: y
  real(8)             :: rh,en
  real(8) :: pr
  integer :: errcode
  errcode=0
 
  call get_eos_userdef(y,xa,xb)

end subroutine eos_pressure_ana0

! *****************************************************************************

subroutine eos_pressure_ana(xa,xb,y,dydxa,dydxb,errcode)
  implicit none

  real(8),intent( in) :: xa,xb
  real(8),intent(out) :: y,dydxa,dydxb
  real(8)             :: rh,en
  real(8) :: pr,dprdrh,dprden
  integer :: errcode
  errcode=0
 
  call get_eos_userdef(y,xa,xb)

  if (eqstate .eq. 1) then
     rh=xa
     en=xb
     pr=y
     if (numerically .eq. 1) then
        call num_derivative_2D(xa,xb,dydxa,dydxb)
     else
        INCLUDE 'eos_data/part_der_pr_vs_rh_en1.def'
        dydxa=dprdrh
        dydxb=dprden
     endif
  else 
     rh=xa
     en=xb
     pr=y
     if (numerically .eq. 1) then
        call num_derivative_2D(xa,xb,dydxa,dydxb)
     else
        INCLUDE 'eos_data/part_der_pr_vs_rh_en2.def'
        dydxa=dprdrh
        dydxb=dprden
     endif
  end if

end subroutine eos_pressure_ana

! *****************************************************************************

subroutine eos_sound_ana(rh,en,pr,cs2,errcode)

  real(8),intent( in) :: rh,en,pr
  real(8),intent(out) :: cs2
  real(8) ::  dprdrh, dprden
  integer, intent(out) :: errcode
  errcode=0
  
     if (numerically .eq. 1) then
        call num_derivative_2D(rh,en,dprdrh,dprden)
     else
       if (eqstate .eq. 1) then
          INCLUDE 'eos_data/part_der_pr_vs_rh_en1.def'
       else
          INCLUDE 'eos_data/part_der_pr_vs_rh_en2.def'
       end if
     endif

  if(en+pr .eq. 0) then
   cs2=dprden
  else
   cs2=dprden+(rh/(en+pr))*dprdrh
  end if
 
  cs2=max(cs2,0.d0) 

end subroutine eos_sound_ana

! *****************************************************************************

subroutine get_derived_data_ana(rh, pr, en, temp, s,errcode )
  implicit none
  real(8),intent(in) :: rh, pr ! pressure
  real(8),intent(out) :: en, temp, s ! en=energy density, temp=temperature, s=entropy_density
  real(8) :: g, b !variables used as parameters inside the temperature.def file
  integer, intent(out) :: errcode

  errcode=0
  
  call eos_energy_ana(rh, en, pr,errcode)
  INCLUDE 'eos_data/temperature.def' ! here we compute the temperature from the pressure
  s=get_entropy_density(rh,en, pr, temp) ! here we compute the entropy density
end subroutine get_derived_data_ana

! *****************************************************************************

subroutine get_derived_temp_for_check_ana(rh, pr, en, temp,error_code)
  implicit none
  real(8),intent(in) :: rh, pr ! pressure
  real(8),intent(out) :: en, temp ! en=energy density, temp=temperature
  real(8) :: g, b !variables used as parameters inside the temperature.def file
  integer, intent(out) :: error_code

  error_code=0
 
  call eos_energy_ana(rh, en, pr,error_code)
  INCLUDE 'eos_data/temperature.def' ! here we compute the temperature from the pressure
end subroutine get_derived_temp_for_check_ana

! *****************************************************************************

real(8) function get_entropy_density(rh, en, pr, temp)
! based on the relation s= (epsilon + p - mu n)/T
  implicit none
  real(8),intent(in) :: rh, en, pr, temp
! AAAAAA for now we don't use rh (i.e. n) and mu, but we should insert them later
  if(temp .eq. 0) then
   get_entropy_density=0.
  else
   get_entropy_density=((en+pr)/temp)
  end if

end function get_entropy_density

! *****************************************************************************

subroutine set_tab_from_an_eos()
  implicit none
  real(8), parameter :: en_min_value=1d-4
  real(8), parameter :: en_max_value=20.0d0
  real(8) :: pr
  real(8) :: en
  real(8) :: rh
  real(8) :: en_start_cycle
  integer, parameter :: step_for_magnitude=2000
  character(15), parameter :: output_filename='tab_an_eos.dat'
  integer :: i !just a counter
  integer :: filerror
  real(8) a, b, c, d, e ! parameters for the user defined analytical equation of state
  
020  format(f14.12,a2,f14.12,a2,f14.12,a2,f14.12)
  
  eos_tab_name=output_filename
  open(22,file='eos_data/'//eos_tab_name,status='REPLACE',iostat=filerror)
      if (filerror .ne. 0) then
      write(*,*) 'eos_data/'//eos_tab_name, " cannot be opened, so I'm forced to quit!"
      close(1)
      call exit(1)
      end if

  en=en_min_value
  do while (en .le. en_max_value)
     i=1
     en_start_cycle=en
     do while (en .le. en_start_cycle*10)
        if(eqstate .eq. 1) then
          INCLUDE 'eos_data/pressure_vs_rh_en1.def'
        else
          INCLUDE 'eos_data/pressure_vs_rh_en2.def'
        end if
        i=i+1
        en=en_start_cycle*i/2000
      
        write(22,020) 0.0,en,pr,0.0
     end do
  end do
  close(22)
end subroutine set_tab_from_an_eos

! *****************************************************************************

subroutine get_en_dens_from_entropy_dens(s, en, min_en_dens, max_en_dens, s_acc,errcode)
  implicit none
  real(8), intent(in) :: s, min_en_dens, max_en_dens, s_acc ! entropy dens., min and max values allowed for energy dens., accuracy
  real(8), intent(out) :: en
  real(8) :: en_try ! guess value for energy density
  real(8) en_low, en_high ! borders for bisections
  real(8) :: s_try, s_low, s_high !  guess values for entropy density
  real(8) :: pr,dprdrh,dprden,temp
  real(8) :: rh=0
  integer :: iterations
  integer, parameter :: max_iterations=200
  real(8) tolerance ! absolute precision required
  integer, intent(out) :: errcode

  tolerance=s*s_acc
  
  ! we make use of routines already written
  ! we provide an initial guess value using the analytic ideal eos pr=en/3, en=a*T^4, where a=PI^2/90/(hc)^3
  ! thus the formula we use is: en=(s^4/3)*(3/4)^(4/3)*1/(a^(1/3))=0.2809*(s^4/3)

  en_try=0.2809*(s**(4./3.))
  en_low=en_try/20.
  en_high=min(en_try*20.,max_en_dens*0.99999999999999)
  
  call eos_pressure(rh,en_try,pr,dprdrh,dprden,errcode) !first, we get the value of pressure corresponding to en_try
  call get_derived_data(rh, pr, en_try, temp,s_try,errcode) !now we get a guess value for entropy density
  call eos_pressure(rh,en_high,pr,dprdrh,dprden,errcode) !first, we get the value of pressure corresponding to en_try
  call get_derived_data(rh, pr, en_high, temp,s_high,errcode) !now we get a value for entropy density at right border
  

  do iterations=1,max_iterations
   if (abs(s_try-s) .le. tolerance) exit
   if (sign(real(1.,8),(s_high-s)) .eq. sign(real(1.,8),(s_try-s))) then
      en_high=en_try
      s_high=s_try
    else
      en_low=en_try
   end if    
      en_try=(en_high+en_low)/2.0
      call eos_pressure(rh,en_try,pr,dprdrh,dprden,errcode) 
      call get_derived_data(rh, pr, en_try, temp,s_try,errcode) 
  end do
   

  if(iterations .eq. max_iterations) then
   write(*,*) "Too many iterations when determining the initial energy density distribution from entropy density distribution"
   write(*,*) "Sorry, but I give up..."
   errcode=44
   return
  end if  
  
  en=en_try

end subroutine get_en_dens_from_entropy_dens

! *****************************************************************************

subroutine eos_energy_ana(rh, en, pr,errcode)
implicit none
real(8), intent(in) :: pr, rh
real(8), intent(out) :: en
real(8) :: min_en, max_en, en_guess
real(8) :: min_pr, max_pr, y
integer :: i
real(8) :: width
real(8) :: dprden, dprdrh 
integer, parameter :: max_iterations=100
integer,intent(out) :: errcode

errcode=0

if(check_on_pressure) then
  if (pr .le. przero) then
      en=enezero*pr/przero
      return
  end if
end if

if (numerically .eq. 1) then
  width=1.0e-1
! now we find an adequate interval for using bisection method
! we start from an initial guess using the ideal equation of state and we expand the interval until we find a change in the sign
! of eos(rh,en)-pr
  do i=1,max_iterations
    min_en=3.0*pr/(1.0+width)
    max_en=3.0*pr*(1.0+width)
  
    call get_eos_userdef(y,rh,min_en)
    min_pr=y-pr
    call get_eos_userdef(y,rh,max_en)
    max_pr=y-pr
    if(( (min_pr .gt. 0.0) .and. (max_pr .gt. 0.0) ) .or. ( (min_pr .lt. 0.0) .and. (max_pr .lt. 0.0))) then
      width=width*2.0
      continue
    else
      exit
    end if
  end do  
  if(i .eq. max_iterations) then
    write(*,*) 'Sorry, but I wasn''t able to find an interval for applying bisection method in the eos_energy_ana subroutine'
    write(*,*) 'I''m leaving...'
    errcode=88
    return
  end if

  do i=1,max_iterations
    en_guess=0.5*(min_en+max_en)
    call get_eos_userdef(y,rh,en_guess)
    if(((y-pr).eq.0) .or. ((max_en-min_en)/2.0 .lt. root_accuracy)) then
      en=en_guess
      exit
    else
      call get_eos_userdef(min_pr,rh,min_en)
      call get_eos_userdef(max_pr,rh,max_en)
      if ( sign(real(1.,8),y-pr) .eq. sign(real(1.,8),min_pr-pr)) then
        min_en=en_guess
      else
        max_en=en_guess
      endif
    endif
  end do
 
  if(i .eq. max_iterations) then
    write(*,*) 'Sorry, but I wasn''t able to find any root when applying bisection method in the eos_energy_ana subroutine'
    write(*,*) 'I''m leaving...'
    errcode=89
    return
  end if

else
  if(eqstate .eq. 1) then
    INCLUDE 'eos_data/energy_den1.def'
  else
    INCLUDE 'eos_data/energy_den2.def'
  end if
endif
  

end subroutine eos_energy_ana

! *****************************************************************************

subroutine eos_energy_with_ders_ana(rh, en, pr,dendrh,dendpr,errcode)
implicit none
real(8), intent(in) :: pr, rh
real(8), intent(out) :: en
real(8) :: dendrh, dendpr 
integer,intent(out) :: errcode

errcode=0

if(check_on_pressure) then
  if (pr .le. przero) then
      en=enezero*pr/przero
      return
  end if
end if

if (numerically .eq. 1) then
   write(*,*) "Sorry, currently it is not possible to run ECHO-QGP with MHD with system solver 5 and numerical derivatives"
   write(*,*) "of an analytic eos."
   call exit(6)
else
  if(eqstate .eq. 1) then
    INCLUDE 'eos_data/energy_den_with_ders1.def'
  else
    INCLUDE 'eos_data/energy_den_with_ders2.def'
  end if
endif
  

end subroutine eos_energy_with_ders_ana

! *****************************************************************************
subroutine find_pressure_from_temperature(temp_input, press_output)
  implicit none
  real(8), intent(in) :: temp_input
  real(8), intent(out) :: press_output
  real(8) :: rh_bb=1.d0, dprdrh_bb, dprden_bb, energy_bb !bitcuket values, just for respecting the argument list when calling eos_pressure
  integer :: errcode !here for compatibility, but not used at this level
  real(8) :: pr_low, pr_high, temp_try, pr_try, en_try
  real(8), parameter :: tolerance=1.d-10

  en_try=min_en_dens
  if(en_try .eq. 0.) en_try=tolerance/2.
  call eos_pressure(rh_bb,en_try,pr_try,dprdrh_bb,dprden_bb,errcode)
  call get_derived_temp_for_check(rh_bb, pr_try, energy_bb, temp_try,errcode)
  do while(temp_try .le. temp_input)
   call get_derived_temp_for_check(rh_bb, pr_try, energy_bb, temp_try,errcode)
   pr_try=pr_try*1.1
  end do

  pr_high=pr_try
  do while(abs(temp_try-temp_input) .gt. temp_input*tolerance)
     pr_try=(pr_low+pr_high)/2.
     call get_derived_temp_for_check(rh_bb, pr_try, energy_bb, temp_try,errcode)
     if(temp_try .ge. temp_input) then
        pr_high=pr_try
     else
        pr_low=pr_try
     end if
  end do
  press_output=pr_try

end subroutine find_pressure_from_temperature


! *****************************************************************************

subroutine eos_pressure_from_s(entropy,press,error_flag)
  implicit none
  integer i
  real(8), intent(in) :: entropy
  real(8), intent(out) :: press
  real(8) :: s_prev, s_next
  integer error_flag
  
  if(eqstate .ne. 3) then
          write(*,*) "Sorry, but currently the subroutine eos_pressure_from_s in eos.f90 works only with a tabuleted EoS. I quit."
          call exit(3)
  end if
  if(entries .eq. 0) then
          write(*,*) "Sorry, but in the subroutine eos_pressure_from_s in eos.f90 I found 0 entries in the tabulated EoS. I quit."
          call exit(3)
  end if

  error_flag=0
  i=1
  s_next=0.
  do while(i .lt. entries)
     s_prev=s_next
     s_next=(energy_density_array(i+1)+pressure_array(i+1))/temperature_array(i+1)
     if((entropy .ge. s_prev) .and. (entropy .lt. s_next)) then
        press=pressure_array(i)+(pressure_array(i+1)-pressure_array(i))*(entropy-s_prev)/(s_next-s_prev)
        return
     end if
     i=i+1
  end do
  !what follows will be executed only in case of failure
  write(*,*) "Sorry, I was not able to obtain the pressure profile from the initial entropy density distribution, so I quit..."
  error_flag=47

end subroutine eos_pressure_from_s

! *****************************************************************************

end module eos

