! ****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Copyright (C) 2015-2019 The ECHO-QGP team                                * 
! *                                                                           *
! *  File: init.f90                                                           *
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
! *  Authors: Gabriele Inghirami (gabriele.g.inghirami@jyu.fi)                *
! *           Francesco Becattini (becattini@fi.infn.it)                      *
! *           Vinod Chandra (vchandra@iitgn.ac.in)                            *
! *           Andrea Beraudo (beraudo@to.infn.it)                             *
! *           Arturo De Pace (depace@to.infn.it)                              *
! *                                                                           *
! *  Contributors: Luca Del Zanna (delzanna@unifi.it)                         *
! *                                                                           *
! *  Acknowledgments: P. Cea, J. Rizzo, L. Cosmai, G. Denicol                 *
! *                                                                           *
! *****************************************************************************

  module qgpinit
  use eos
  use parallel
  use common
  use viscous
  use glaubermc
  
  implicit none

  real(8)                  :: xx,yy,eta
  real(8),allocatable,dimension(:)  :: thick !thickness vector
  integer, parameter    :: minus=-1, center=0, plus=1 !flag/parameters for changing the simmetry of energy distribution

  contains

! ***************************************************************************

  subroutine print_param()
        implicit none
        real(8) :: s_input, en_dens_output
        integer ferror
        character(len=120) :: outbuffer
        

010 format(a50,2x,i6)
011 format(a50,2(2x,i6))
012 format(a50,3(2x,i6))
013 format(a50,2x,e14.7)
014 format(a50,2(2x,e14.7))
015 format(a50,3(2x,e14.7))
016 format(a51,2x,e14.7)

        if(ipe .eq. 0) then

          write(*,*) '  Using as parameters file: '//param_file
          write(*,*) '  Settings and parameters:'
          write(*,*) 

          select case(init_type)
            case(GLAUBER_GEOMETRIC)
              write(*,"(2x,a4,1x,i2,2x,a44)") 'Test',GLAUBER_GEOMETRIC,'- optical-geometrical Glauber initialization'
            case(SHOCK_TUBE_2D)
              write(*,"(2x,a4,1x,i2,2x,a27)") 'Test',SHOCK_TUBE_2D,'- 2D+1 shock tube init_type'
            case(SHEAR_VISCOUS_1D)
              write(*,"(2x,a4,1x,i2,2x,a33)") 'Test',SHEAR_VISCOUS_1D,'- 1D viscous shear flow init_type'
            case(GLAUBER_MONTECARLO)
              if(avg_events .eq. 0) then
                write(*,"(2x,a4,1x,i2,2x,a37)") 'Test',GLAUBER_MONTECARLO,'- Glauber-Monte Carlo initial profile'
              else if(avg_events .eq. 1) then
                write(*,"(2x,a4,1x,i2,2x,a54)") 'Test',GLAUBER_MONTECARLO,'- run with averaged Glauber-Monte Carlo in. conditions'
              else
                write(*,*) "  Sorry, but it is not clear whether I should run event by event simulations"
                write(*,*) "  or average many initial conditions and perform just 1 simulation. I quit."
                call exit(1)
              end if
            case(GUBSER_IDEAL)
              write(*,"(2x,a4,1x,i2,2x,a31)") 'Test',GUBSER_IDEAL,'- ideal Gubser flow init_type'
            case(GUBSER_VISCOUS)
              write(*,"(2x,a4,1x,i2,2x,a31)") 'Test',GUBSER_VISCOUS,'- viscous Gubser flow init_type'
            case(TABULATED_INIT)
              write(*,"(2x,a4,1x,i2,2x,a53)") 'Test',TABULATED_INIT,'- tabulated initial entropy or energy density profile'
            case(ALFVEN2D)
              write(*,"(2x,a4,1x,i2,2x,a20)") 'Test',ALFVEN2D,'- alfven 2D mhd test'
            case(BLAST3D)
              write(*,"(2x,a4,1x,i2,2x,a23)") 'Test',BLAST3D,'- blast wave 3D mhd test'
            case(MHD1DBJ)
              write(*,"(2x,a4,1x,i2,2x,a28)") 'Test',MHD1DBJ,'- 0D+1 bjorken flow mhd test'
            case(SODMHD)
              write(*,"(2x,a4,1x,i2,2x,a18)") 'Test',SODMHD,'- 1 D Sod mhd test'
            case(LYU)
              write(*,"(2x,a4,1x,i2,2x,a24)") 'Test',LYU,'- 1 D Lyutikov rmhd test'
            case(BERHAD)
              write(*,"(2x,a4,1x,i2,2x,a31)") 'Test',BERHAD,'- 1 D entropy density rmhd test'
            case(ROTOR)
              write(*,"(2x,a4,1x,i2,2x,a36)") 'Test',ROTOR,'- 2D or 3D Rotor rmhd test'
            case(OT)
              write(*,"(2x,a4,1x,i2,2x,a42)") 'Test',OT,'- 2D Orszang-Tang test'
            case(DBG)
              write(*,"(2x,a4,1x,i2,2x,a42)") 'Test',DBG,'- 2D grid with B field on transverse plane'
            case(EXT_CONS)
              write(*,"(2x,a4,1x,i2,2x,a43)") 'Test',EXT_CONS,'- external file with conservative variables'
            case(EXT_GLISSANDO)
              write(*,"(2x,a4,1x,i2,2x,a35)") 'Test',EXT_GLISSANDO,'- external file with Glissando data'
             
          end select 

          write(*,"(i1)",advance='no') coordinates_option
          if(coordinates .eq. MINKOWSKI) then
            write(*,*) ' Using Minkowski coordinates'
          else
            write(*,*)' Using Bjorken coordinates'
          end if

        end if

        if(ipe .eq. 0) then
          write(*,"(i1)",advance='no') boundary_cond_option
        end if

        if(boundary_cond .eq. 0) then
          if(ipe .eq. 0) write(*,*) ' Using periodic boundary conditions'
          ibc=0
        else if (boundary_cond .eq. 2) then
          if(ipe .eq. 0) write(*,*) ' Using outflow 0th order boundary conditions'
          ibc=2
        else 
          if(ipe .eq. 0) then
            write(*,*) ' Sorry, but ECHO-QGP has been tested only with 0 (periodic) or 2 (outflow) boundary conditions'
            write(*,*) ' I am not sure to provide correct results, so I quit.'
          end if
          call exit(2)
        end if

        !-- Redefine BCs if parallel
        if (prl) then
                 if (ipe/=    0 .or. ibc(1,1,1)==0) ibc(:,1,1)=-2
                 if (ipe/=npe-1 .or. ibc(1,2,1)==0) ibc(:,2,1)=-2
        end if


        if(ipe .eq. 0) then
          write(*,"(i1)",advance='no') mhd_flag_option
          if(mhd) then
            if(rmhd) then
              write(*,*) ' This is a resistive MHD simulation'
            else
              write(*,*) ' This is an ideal MHD simulation'
            end if
            write(*,"(i1)",advance='no') divclean_flag_option
            if(divclean) then
              write(*,*) " with Dedner's method for divergence cleaning enabled"
              write(*,"(i1)",advance='no') glm_alpha_option
              write(*,013) "  with a value of the dumping parameter glm_alpha: ",glm_alpha
            else
              write(*,*) " without any method to control the solenoidal condition"
            end if
          else
            write(*,*) ' This is an HD simulation'
          end if
          write(*,"(i1)",advance='no') viscosity_flag_option
          if(viscosity) then
            write(*,*) ' This is a viscous simulation'

            write(*,"(i1)",advance='no') bulkvis_flag_option
            if(bulkvis) then
              write(*,*) ' Bulk viscosity is taken into account'
            else
              write(*,*) ' Bulk viscosity is neglected'
            end if

            write(*,"(i1)",advance='no') obtained_option
            if(obtained .eq. 'no') then
              write(*,*) ' Evolved shear viscous tensor components: xx, yy, zz, xy, xz, zz'
              write(*,*) '  tt, tx, ty and tz components are obtained imposing orthogonality'
            else
              write(*,*) ' Evolved shear viscous tensor components: xx, yy, zz, xy, xz'
              write(*,*) '  tt, tx, ty, tz and zz are obtained imposing orth. and null trace'
            end if

          else
            write(*,*) ' This is an unviscid simulation'
          end if

          if(prl) then
            write(*,*) '  This simulation uses MPI and the number of processors is:',npe
          else
            write(*,*) '  This is a serial run (no MPI)'
          end if
            
          write(*,*) 

          write(*,*) '  Grid parameters:'
          write(*,"(i1)",advance='no') nx_option
          write(*,010) '  number of cells along x direction:              ', nx
          write(*,"(i1)",advance='no') ny_option
          write(*,010) '  number of cells along x direction:              ', ny
          write(*,"(i1)",advance='no') nz_option
          write(*,010) '  number of cells along z (or eta)  direction:    ', nz
          write(*,"(i1)",advance='no') x1_option
          write(*,013) '  x min (fm):                                     ',x1
          write(*,"(i1)",advance='no') x2_option
          write(*,013) '  x max (fm):                                     ',x2
          write(*,"(i1)",advance='no') y1_option
          write(*,013) '  y min (fm):                                     ',y1
          write(*,"(i1)",advance='no') y2_option
          write(*,013) '  y max (fm):                                     ',y2
          write(*,"(i1)",advance='no') z1_option
          write(*,013) '  z min (fm):                                     ',z1
          write(*,"(i1)",advance='no') z2_option
          write(*,013) '  z max (fm):                                     ',z2
          write(*,*)
 
          write(*,*) '  Time parameters:'
          write(*,"(i1)",advance='no') tmin_option
          write(*,013) '  starting time:                                  ', tmin 
          write(*,"(i1)",advance='no') tmax_option
          write(*,013) '  ending time:                                    ', tmax 
          write(*,"(i1)",advance='no') temp_end_option
          write(*,013) '  ending temeperature (MeV):                      ', temp_end*1000.
          write(*,"(i1)",advance='no') maxdt_option
          write(*,013) '  maximum timestep:                               ', maxdt
          write(*,*)
 
          if(viscosity) then
            write(*,*) '  Viscous parameters:'
            write(*,"(i1)",advance='no') eta_over_s_option
            write(*,013) '  eta/s parameter for shear viscosity tensor:     ', eta_over_s
            write(*,"(i1)",advance='no') smooth_temp_option
            if(enable_smooth_viscosity) write(*,013) '  Temperature limit for smoothing viscosity:      ', smooth_temp
            write(*,*)
          end if
          write(*,"(i1)",advance='no') eqstate_option
          write(*,010) '  Equation of state:                              ', eqstate
          if(eqstate .eq. 3) then
            write(*,"(i1)",advance='no') eos_tab_name_option
            write(*,*) ' File with tabulated eos:                           '//adjustl(eos_tab_name)
          end if
          write(*,"(i1)",advance='no') numerically_option
          write(*,010) '  numerical derivatives with analytic eos:        ', numerically
          write(*,*)
 
          if((init_type .eq. GLAUBER_GEOMETRIC) .or. (init_type .eq. GLAUBER_MONTECARLO) .or. (init_type .eq. TABULATED_INIT) &
         & .or. (init_type .eq. DBG)) then
            write(*,*) '  Nucleus parameters:'
            write(*,"(i1)",advance='no') nuclei_data_file_option
            write(*,*) ' file with nuclear parameters:                      '//adjustl(nuclei_data_file)
            write(*,016) '   nucleus mass (GeV):                             ',projmass
            write(*,016) '   nuclear radius (fm):                            ', radius
            write(*,016) '   W.-S. width (fm):                               ',delta
            write(*,"(i1)",advance='no') rads_option
            write(*,013) '  sqrt(s) (GeV):                                  ',rads
            write(*,016) '   beam rapidity Y_b:                              ',yb
            write(*,"(i1)",advance='no') sigma_in_option
            write(*,013) '  cross section (mb):                             ',sigma_in
            write(*,"(i1)",advance='no') bimpact_option
            write(*,013) '  impact parameter (fm):                          ',bimpact
            write(*,"(i1)",advance='no') ah_option
            write(*,013) '  initial hardness parameter:                     ',ah

 
            write(*,"(i1)",advance='no') ecenter_option
            if(ienentr .eq. 0) then
              write(*,013) '  central energy density (GeV/fm^3):              ',ecenter
            else if(ienentr .eq. 1) then
              write(*,013) '  central entropy density (fm^-3):                ',ecenter
              if(eqstate .ne. 3) then
                write(*,*) " Sorry, but currently you can initialize with entropy density only if you use a tabulated EoS..."
                call exit(3)
              end if
            else
              write(*,*) " Sorry, but I don't understand the value of the IENENTR parameter..."
              call exit(3)
            end if
          
            write(*,"(i1)",advance='no') enezero_option
            write(*,013) '  enezero (GeV/fm^3):                             ',enezero
            write(*,016) '  przero (GeV/fm^3):                              ',przero
            write(*,"(i1)",advance='no') rhocenter_option
            write(*,013) '  central density:                                ',rhocenter
            write(*,"(i1)",advance='no') deta_option
            write(*,013) '  rapidity distrib. shift (deta or etaflat par.): ',deta
            write(*,"(i1)",advance='no') sigeta_option
            write(*,013) '  rapidity distrib. width (sigeta parameter):     ',sigeta
            write(*,*)
          end if

          if(init_type .eq. GLAUBER_GEOMETRIC) then
            write(*,"(i1)",advance='no') ueta_coeff_option
            write(*,013) '  ueta A coeff. so that u^eta=A*x:                ', ueta_coeff
          end if
          write(*,*)

          if(init_type .eq. EXT_CONS) then
            write(*,"(i1)",advance='no') ext_cons_file_option
            write(*,013) '  File containing initial conservative variables:  '//adjustl(ext_cons_file)
          end if
          write(*,*)

          if(init_type .eq. EXT_GLISSANDO) then
            write(*,"(i1)",advance='no') ext_glissando_file_option
            write(*,013) '  File containing initial variables from Glissando:  '//adjustl(ext_glissando_file)
          end if
          write(*,*)

          if(mhd) then
            write(*,"(i1)",advance='no') input_B_field_option
            write(*,*) ' File with initial magnetic field configuration:    '//adjustl(input_B_field)
            write(*,"(i1)",advance='no') B_amplify_factor_option
            write(*,013) '  Initial B field amplification factor:            ', B_amplify_factor
            write(*,"(i1)",advance='no') dump_initial_B_flag_option
            if(dump_initial_B) then
              write(*,*) ' Initial B field dumping enabled'
            else
              write(*,*) ' Initial B field dumping disabled'
            end if
            if(init_type .eq. GLAUBER_MONTECARLO) then
              write(*,"(i1)",advance='no') magfield_type_option
              if(magfield_type .eq. 1 ) then
                write(*,*) ' The initial B field will have both classical and chiral sources'
              else if(magfield_type .eq. 2 ) then
                write(*,*) ' The initial B field will have classical origin only'
              else if(magfield_type .eq. 3 ) then
                write(*,*) ' The initial B field will have chiral origin only'
              else
                write(*,*) ' Sorry, but I cannot understand which kind of initial magnetic field you want'
                write(*,*) '  (i.e. classical, chiral or both)'
                call exit(2)
              end if
              if((magfield_type .eq. 1) .or. (magfield_type .eq. 2)) then
                write(*,"(i1)",advance='no') sigma_el_cond_option
                write(*,013) '  Electrical cond. of the medium to compute in. B: ', sigma_el_cond 
              end if
              if((magfield_type .eq. 1) .or. (magfield_type .eq. 3)) then
                write(*,"(i1)",advance='no') sigma_chiral_cond_option
                write(*,013) '  Chiral cond. of the medium to compute initial B: ', sigma_chiral_cond 
              end if
            end if
            write(*,"(i1)",advance='no') edens_for_B_dump_option
            write(*,013) '  En. dens. limit for the B field supp.(GeV/fm^3):', edens_for_B_dump 
            write(*,"(i1)",advance='no') B_th_press_ratio_option
            write(*,013) '  Max ratio between magnetic and thermal press.:  ', B_th_press_ratio 
            if((out_sel%rc .eq. 1) .and. out_freeze) then
              print_rho_comov=.true.
              write(*,*) "Printing electric charge density in comoving frame on the freezeout hypersurface"
            end if
          end if
          write(*,*)

          if(init_type .eq. GLAUBER_MONTECARLO) then!for Glauber-MonteCarlo 
            write(*,*) ' Glauber-MonteCarlo specific parameters:'
            write(*,"(i1)",advance='no') kind_of_collision_option
            select case(kind_of_collision)
              case(1)
                 write(*,*) ' kind of collision: nucleus-nucleus'
              case(2)
                 write(*,*) ' kind of collision: deuton-nucleus'
              case(3)
                 write(*,*) ' kind of collision: proton-nucleus'
            end select
            write(*,"(i1)",advance='no') nconf_option
            write(*,010) '  number of nuclear configurations:               ',nconf
            if(fixed_b) then
              write(*,"(i1)",advance='no') fixed_b_flag_option
              write(*,*) ' Using a fixed impact parameter'
            else
              write(*,"(i1)",advance='no') fixed_b_flag_option
              write(*,*) ' Using a randomly sampled impact parameter'
              write(*,"(i1)",advance='no') nbcoll_option
              write(*,010) '  number of impact parameters per configuration:  ', nbcoll
            end if
            write(*,"(i1)",advance='no') events_start_option
            write(*,010) '  id number of starting event:                    ', events_start
            write(*,"(i1)",advance='no') events_stop_option
            write(*,010) '  id number of ending event:                      ', events_stop
            write(*,"(i1)",advance='no') kappa_option
            write(*,013) '  k model parameters:                             ',kappa
            write(*,"(i1)",advance='no') sig_mc_option
            write(*,013) '  sigma model smearing parameter:                 ',sig_mc
            write(*,"(i1)",advance='no') min_participants_option
            write(*,010) '  min number of participants to accept the event: ',min_participants
            if(kind_of_collision .eq. 3) then
              write(*,"(i1)",advance='no') eprot_option
              write(*,013) '  proton beam energy:                             ',eprot
            end if
          end if
          write(*,*)

          if(out_freeze) then !for freeze-out hypersurface computation
            write(*,"(i1)",advance='no') out_freeze_flag_option
            write(*,*) ' Computing freeze-out hypersurface'
            if(freeze_type .eq. 0) then
              write(*,"(i1)",advance='no') freeze_type_option
              write(*,*) ' hypersurface computation based on temperature'
              write(*,"(i1)",advance='no') freeze_value_option
              write(*,013) '  freezeout threshold (MeV):                      ', freeze_value*1000.
            else
              write(*,"(i1)",advance='no') freeze_type_option
              write(*,*) ' hypersurface computation based on energy density'
              write(*,"(i1)",advance='no') freeze_value_option
              write(*,013) '  freezeout threshold (GeV/f^3):                  ', freeze_value
            end if
            write(*,"(i1)",advance='no') freeze_time_interval_option
            write(*,013) '  time interval between hypersurf. computations:  ', freeze_time_interval
          end if
          if(print_rho_comov) then
            write(*,"(i1)",advance='no') print_rho_comov_flag_option
            write(*,*) ' Computing the electric charge density in the comoving frame on f.o. hypersurface'
          end if
          write(*,*)

 
          write(*,*) ' Other numerical parameters:'
          
          if(mhd) then
            write(*,"(i1)",advance='no') algo_solver_option
            write(*,010) '  Inv. algorithm for obtaining prim. variables:   ',algo_solver
          end if 
          write(*,"(i1)",advance='no') cfl_option
          write(*,013) '  Courant-Fr.-Lew. condition parameter:           ',cfl
          write(*,"(i1)",advance='no') nrk_option
          if(nrk .eq. 2) then
            write(*,*) ' Time integration algorithm:                        RK2'
          else if(nrk .eq. 3) then
            write(*,*) ' Time integration algorithm:                        RK3'
          else 
            write(*,*) ' Time integration algorithm:                        IMEX-SSP3(4,3,3)'
          end if
          write(*,"(i1)",advance='no') recal_option
          write(*,*) ' Reconstruction algorithm:                          ',recal
          write(*,"(i1)",advance='no') maxspeed_option
          write(*,013) '  If the velocity is > c, the speed is reduced to:',maxspeed
          write(*,*)

          write(*,*) ' Output parameters:'
          write(*,"(i1)",advance='no') dtlog_option
          write(*,013) '  interval between log updating:                  ', dtlog
          write(*,"(i1)",advance='no') dtout_option
          write(*,013) '  interval between output printing:               ', dtout
          write(*,"(i1)",advance='no') output_precision_option
          if(output_precision .eq. 8) then
            write(*,*) ' output precision:', 'double - 8 bytes'
          else
            write(*,*) ' output precision:', 'float - 4 bytes'
          end if
          write(*,*)
 
 
          write(*,*) ' Variables printed in the output files:'
          write(*,"(i1)",advance='no') density_option
          if(out_sel%density .eq. 1) write(*,*) ' density'
          write(*,"(i1)",advance='no') vx_option
          if(out_sel%vx .eq. 1) write(*,*) ' vx'
          write(*,"(i1)",advance='no') vy_option
          if(out_sel%vy .eq. 1) write(*,*) ' vy'
          write(*,"(i1)",advance='no') vz_option
          if(out_sel%vz .eq. 1) write(*,*) ' vz'
          write(*,"(i1)",advance='no') pressure_option
          if(out_sel%pressure .eq. 1) write(*,*) ' pressure'
          write(*,"(i1)",advance='no') energy_density_option
          if(out_sel%energy_density .eq. 1) write(*,*) ' energy density'
          write(*,"(i1)",advance='no') temperature_option
          if(out_sel%temperature .eq. 1) write(*,*) ' temperature'
          write(*,"(i1)",advance='no') entropy_density_option
          if(out_sel%entropy_density .eq. 1) write(*,*) ' entropy density'
          if(viscosity) then
            write(*,"(i1)",advance='no') bulk_option
            if(out_sel%bulk .eq. 1) then 
              write(*,*) ' bulk viscosity'
            else
              write(*,*) ' no bulk viscosity'
            end if
            write(*,"(i1)",advance='no') pitt_option
            if(out_sel%pitt .eq. 1) then 
              write(*,*) ' pi^tt'
            else
              write(*,*) ' no pi^tt'
            end if
            write(*,"(i1)",advance='no') pitx_option
            if(out_sel%pitx .eq. 1) then 
              write(*,*) ' pi^tx'
            else
              write(*,*) ' no pi^tx'
            end if
            write(*,"(i1)",advance='no') pity_option
            if(out_sel%pity .eq. 1) then
              write(*,*) ' pi^ty'
            else
              write(*,*) ' no pi^ty'
            end if
            write(*,"(i1)",advance='no') pitz_option
            if(out_sel%pitz .eq. 1) then
              write(*,*) ' pi^tz'
            else
              write(*,*) ' no pi^tz'
            end if
            write(*,"(i1)",advance='no') pixy_option
            if(out_sel%pixy .eq. 1) then
              write(*,*) ' pi^xy'
            else
              write(*,*) ' no pi^xy'
            end if
            write(*,"(i1)",advance='no') pixz_option
            if(out_sel%pixz .eq. 1) then
              write(*,*) ' pi^xz'
            else
              write(*,*) ' no pi^xz'
            end if
            write(*,"(i1)",advance='no') piyz_option
            if(out_sel%piyz .eq. 1) then
              write(*,*) ' pi^yz'
            else
              write(*,*) ' no pi^yz'
            end if
            write(*,"(i1)",advance='no') pixx_option
            if(out_sel%pixx .eq. 1) then
              write(*,*) ' pi^xx'
            else
              write(*,*) ' no pi^xx'
            end if
            write(*,"(i1)",advance='no') piyy_option
            if(out_sel%piyy .eq. 1) then
              write(*,*) ' pi^yy'
            else
              write(*,*) ' no pi^yy'
            end if
            write(*,"(i1)",advance='no') pizz_option
            if(out_sel%pizz .eq. 1) then
              write(*,*) ' pi^zz'
            else
              write(*,*) ' no pi^zz'
            end if
          end if
          write(*,"(i1)",advance='no') v0_option
          if(out_sel%v0 .eq. 1) write(*,*) ' u0 or gamma Lorentz factor'
          if(mhd) then
            write(*,"(i1)",advance='no') bx_option
            if(out_sel%bx .eq. 1) write(*,*) ' bx'
            write(*,"(i1)",advance='no') by_option
            if(out_sel%by .eq. 1) write(*,*) ' by'
            write(*,"(i1)",advance='no') bz_option
            if(out_sel%bz .eq. 1) write(*,*) ' bz'
            write(*,"(i1)",advance='no') ex_option
            if(out_sel%ex .eq. 1) write(*,*) ' ex'
            write(*,"(i1)",advance='no') ey_option
            if(out_sel%ey .eq. 1) write(*,*) ' ey'
            write(*,"(i1)",advance='no') ez_option
            if(out_sel%ez .eq. 1) write(*,*) ' ez'
            write(*,"(i1)",advance='no') glm_option
            if(out_sel%glm .eq. 1) write(*,*) ' glm'
            write(*,"(i1)",advance='no') rc_option
            if(out_sel%rc .eq. 1) write(*,*) ' rc'
          end if
          if(derivatives_out) then
            write(*,"(i1)",advance='no') derivatives_flag_option
            write(*,*) ' derivatives will also be printed into separated output files'
          end if
          if(flows_out) then
            write(*,"(i1)",advance='no') flows_flag_option
            write(*,*) ' flows will also be printed into separated output files'
          end if
          
        end if !it closes (ipe .eq. 0)


      end subroutine print_param

!---------------------------------------------------------------

subroutine print_config_summary
  implicit none
  integer ferror
  character(len=120) :: outbuffer
  010 format(a50,2x,i6)
  011 format(a50,2(2x,i6))
  012 format(a50,3(2x,i6))
  013 format(a50,2x,e14.7)
  014 format(a50,2(2x,e14.7))
  015 format(a50,3(2x,e14.7))
  

          open(unit=25,file="config_summary.dat",status='replace',IOSTAT=ferror)
          if(ferror .ne. 0) then
             write(*,*) "Sorry, but I cannot open the file config_summary.dat"
             write(*,*) "Although this is not an essential file, nevertheless this means that there are problems in writing to HD"
             write(*,*) "and therefore I prefer to quit now and let you solve this issue. Bye!"
             call exit(2)
          end if
          write (outbuffer, '(i2)') init_type
          write(25,"(a10,a2)") "INIT_TYPE=",adjustl(outbuffer)
          write(25,"(a10,i1)") "COORD....=",coordinates
          if(viscosity) then
            write(25,"(a11)") "VISCOUS..=1"
          else
            write(25,"(a11)") "VISCOUS..=0"
          end if
          if(bulkvis) then
            write(25,"(a11)") "BULK.....=1"
          else
            write(25,"(a11)") "BULK.....=0"
          end if
          if(mhd) then
            write(25,"(a11)") "MHD......=1"
          else
            write(25,"(a11)") "MHD......=0"
          end if
          if(divclean) then
            write(25,"(a11)") "DIVCLEAN.=1"
          else
            write(25,"(a11)") "DIVCLEAN.=0"
          end if
          write (outbuffer, '(f14.4)') glm_alpha
          write(25,"(a10,a16)") "GLM_PARAM=",adjustl(outbuffer)
          if(dump_initial_B) then
            write(25,"(a11)") "DUMP_IN_B=1"
          else
            write(25,"(a11)") "DUMP_IN_B=0"
          end if
          write (outbuffer, '(f14.8)') edens_for_B_dump
          write(25,"(a10,a16)") "B_DUM_EN.=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') B_th_press_ratio
          write(25,"(a10,a16)") "Bp_ov_Tp.=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') smooth_temp
          write(25,"(a10,a16)") "CUT_TEMP.=",adjustl(outbuffer)
          write (outbuffer, '(i6)') nx
          write(25,"(a10,a6)") "NX.......=",adjustl(outbuffer)
          write (outbuffer, '(i6)') ny
          write(25,"(a10,a6)") "NY.......=",adjustl(outbuffer)
          write (outbuffer, '(i6)') nz
          write(25,"(a10,a6)") "NZ.......=",adjustl(outbuffer)
          write(*,*) "  config_summary.dat written"
          write (outbuffer, '(f14.8)') x1
          write(25,"(a10,a16)") "XMIN.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') x2
          write(25,"(a10,a16)") "XMAX.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') y1
          write(25,"(a10,a16)") "YMIN.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') y2
          write(25,"(a10,a16)") "YMAX.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') z1
          write(25,"(a10,a16)") "ZMIN.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') z2
          write(25,"(a10,a16)") "ZMAX.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') tmin
          write(25,"(a10,a16)") "TSTART...=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') tmax
          write(25,"(a10,a16)") "TSTOP....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') temp_end
          write(25,"(a10,a16)") "TEMP_END.=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') dtlog
          write(25,"(a10,a16)") "DTLOG....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') dtout
          write(25,"(a10,a16)") "DTOUT....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') output_precision
          write(25,"(a10,a1)") "OUTP_PREC=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') maxdt
          write(25,"(a10,a16)") "MAXDT....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') restart_type
          write(25,"(a10,a1)")  "RESTART..=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') cfl
          write(25,"(a10,a16)") "CFL......=",adjustl(outbuffer)
          write(25,"(a10,a5)") "REC_ALGO.=",adjustl(recal)
          write (outbuffer, '(i1)') algo_solver
          write(25,"(a10,a1)") "ALG_SOLV.=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') maxspeed
          write(25,"(a10,a16)") "MAXSPEED.=",adjustl(outbuffer)
          write (outbuffer, '(a60)') nuclei_data_file
          write(25,"(a10,a60)") "NUC_DATAF=",adjustl(outbuffer)
          write(25,"(a10,a5)") "NUCLEUS..=",adjustl(nucleus)
          write (outbuffer, '(f14.8)') rads
          write(25,"(a10,a16)") "RADS.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') sigma_in*10.
          ! sigma has been transformed into fm^2, but it was given in mb = 1/10 fm^2 -> we multiply by 10 to get mb again
          write(25,"(a10,a16)") "SIGMA_IN.=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') bimpact
          write(25,"(a10,a16)") "B........=",adjustl(outbuffer)
          write (outbuffer, '(i1)') ienentr
          write(25,"(a10,a1)") "IENENTR..=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') ah
          write(25,"(a10,a16)") "AH.......=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') ecenter
          write(25,"(a10,a16)") "ECENTER..=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') enezero
          write(25,"(a10,a16)") "ENEZERO..=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') rhocenter
          write(25,"(a10,a16)") "RHOCENTER=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') deta
          write(25,"(a10,a16)") "DETA.....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') sigeta
          write(25,"(a10,a16)") "SIGETA...=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') eta_over_s
          write(25,"(a10,a16)") "ETA_S....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') tau_pi_coeff
          write(25,"(a10,a16)") "TAU_PI_C.=",adjustl(outbuffer)
          write(25,"(a10,a2)") "TRACE_IMP=",adjustl(obtained)
          write (outbuffer, '(i1)') eqstate
          write(25,"(a10,a1)") "EOS......=",adjustl(outbuffer)
          write (outbuffer, '(a60)') eos_tab_name
          write(25,"(a10,a60)") "EOS_FILE.=",adjustl(outbuffer)
          write (outbuffer, '(i1)') numerically
          write(25,"(a10,a1)") "NUM_DER..=",adjustl(outbuffer)
          write (outbuffer, '(i6)') nconf
          write(25,"(a10,a6)") "NCONF....=",adjustl(outbuffer)
          write (outbuffer, '(i6)') nbcoll
          write(25,"(a10,a6)") "NBCOLL...=",adjustl(outbuffer)
          if(fixed_b) then
            write(25,"(a11)") "FIXED_B..=1"
          else
            write(25,"(a11)") "FIXED_B..=0"
          end if
          write (outbuffer, '(i6)') events_start
          write(25,"(a10,a6)") "EV_START.=",adjustl(outbuffer)
          write (outbuffer, '(i6)') events_stop
          write(25,"(a10,a6)") "EV_STOP..=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') kappa
          write(25,"(a10,a16)") "KAPPA....=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') sig_mc
          write(25,"(a10,a16)") "SIG......=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') eprot
          write(25,"(a10,a16)") "EPROT....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') avg_events
          write(25,"(a10,a1)") "AVG_IC...=",adjustl(outbuffer)
          write (outbuffer, '(i3)') min_participants
          write(25,"(a10,a6)") "MIN_PART.=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') pw
          write(25,"(a10,a16)") "PW.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') kind_of_collision
          write(25,"(a10,a1)") "COLLISION=",adjustl(outbuffer)
          if(out_freeze) then
            write(25,"(a11)") "HYP_COMPU=1"
          else
            write(25,"(a11)") "HYP_COMPU=0"
          end if
          write (outbuffer, '(i1)') freeze_type
          write(25,"(a10,a1)") "FREEZKIND=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') freeze_value
          write(25,"(a10,a16)") "FREEZEVAL=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') freeze_time_interval
          write(25,"(a10,a16)") "HYPSURFTI=",adjustl(outbuffer)
          write (outbuffer, '(a60)') input_edf
          write(25,"(a10,a60)") "IN_D_FILE=",adjustl(outbuffer)
          write (outbuffer, '(a60)') ext_cons_file
          write(25,"(a10,a60)") "EX_C_FILE=",adjustl(outbuffer)
          write (outbuffer, '(a60)') ext_glissando_file
          write(25,"(a10,a60)") "EX_G_FILE=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') etam
          write(25,"(a10,a16)") "ETAM_TILT=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') ueta_coeff
          write(25,"(a10,a16)") "UETA_COEF=",adjustl(outbuffer)
          write (outbuffer, '(a60)') input_B_field
          write(25,"(a10,a60)") "B_in_FILE=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') sigma_el_cond
          write(25,"(a10,a16)") "EL_COND..=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') sigma_chiral_cond
          write(25,"(a10,a16)") "CHIR_COND=",adjustl(outbuffer)
          write (outbuffer, '(i1)') magfield_type
          write(25,"(a10,a1)") "MAGFIELDT=",adjustl(outbuffer)
          write (outbuffer, '(f14.8)') B_amplify_factor
          write(25,"(a10,a16)") "B_amp_fac=",adjustl(outbuffer)
          write (outbuffer, '(i1)') boundary_cond
          write(25,"(a10,a1)") "bound_con=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%density
          write(25,"(a10,a1)") "density..=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%vx
          write(25,"(a10,a1)") "vx.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%vy
          write(25,"(a10,a1)") "vy.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%vz
          write(25,"(a10,a1)") "vz.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pressure
          write(25,"(a10,a1)") "pressure.=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%energy_density
          write(25,"(a10,a1)") "ene_dens.=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%temperature
          write(25,"(a10,a1)") "temper...=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%entropy_density
          write(25,"(a10,a1)") "entr_dens=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%bulk
          write(25,"(a10,a1)") "bulk_visc=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pitt
          write(25,"(a10,a1)") "pi^tt....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pitx
          write(25,"(a10,a1)") "pi^tx....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pity
          write(25,"(a10,a1)") "pi^ty....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pitz
          write(25,"(a10,a1)") "pi^tz....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pixy
          write(25,"(a10,a1)") "pi^xy....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pixz
          write(25,"(a10,a1)") "pi^xz....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%piyz
          write(25,"(a10,a1)") "pi^yz....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pixx
          write(25,"(a10,a1)") "pi^xx....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%piyy
          write(25,"(a10,a1)") "pi^yy....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%pizz
          write(25,"(a10,a1)") "pi^zz....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%v0
          write(25,"(a10,a1)") "gamma....=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%bx
          write(25,"(a10,a1)") "bx.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%by
          write(25,"(a10,a1)") "by.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%bz
          write(25,"(a10,a1)") "bz.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%ex
          write(25,"(a10,a1)") "ex.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%ey
          write(25,"(a10,a1)") "ey.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%ez
          write(25,"(a10,a1)") "ez.......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%glm
          write(25,"(a10,a1)") "glm......=",adjustl(outbuffer)
          write (outbuffer, '(i1)') out_sel%rc
          write(25,"(a10,a1)") "rc.......=",adjustl(outbuffer)
          if(derivatives_out) then
            write(25,"(a11)") "derivativ=1"
          else
            write(25,"(a11)") "derivativ=0"
          end if
          if(flows_out) then
            write(25,"(a11)") "flows....=1"
          else
            write(25,"(a11)") "flows....=0"
          end if
          
          
       close(25)

end subroutine print_config_summary

! ***************************************************************************

    subroutine initio (nu,xini,vv,error_flag)
    implicit none
    integer, intent(in):: nu
    real(8),intent(in)    :: xini(3)
        real(8) :: vv(nu)
    real(8) ymax,tr_weight
    real(8) rho,edens,press,dpr1,dpr2
    real(8) velo(3)
       
        integer error_flag
  
        error_flag=0
 
    xx =xini(1)
    yy =xini(2)
    eta=xini(3)
    
        !call calc_enedens(edens,errcode) !it matches the version of calc_enedens provided by A.DP
        call calc_enedens(edens)
        call calc_chargedens(rho)
        if(init_type .eq. GLAUBER_GEOMETRIC) call calc_vel_long_nz_ueta(velo,xx,yy,eta)

        if(ienentr .eq. 0) then
            call eos_pressure(rho,edens,press,dpr1,dpr2,error_flag)
        else
            call eos_pressure_from_s(edens,press,error_flag)
        end if

        press=press+przero
           
        if(error_flag .gt. 0) then
            write(*,*) "An error occurred into module init, subroutine initio, when computing the pressure profile..."
            errcode=error_flag
            write(*,*) "Error code:",errcode
            write(*,*) "Position: x=",xx,"y=",yy,"z=",eta
            run_crashed=.true.
            return
        end if
        vv(1)=rho
        if(.not. mhd) then
            vv(2)=velo(1)    
            vv(3)=velo(2)   
            vv(4)=velo(3)
        end if
        vv(5)=press

    return
    end subroutine initio
    
! ***************************************************************************
    subroutine calc_enedens(edens)
    implicit none
    real(8) edens,tt1,tt2,H, HP, HM, thplus, thminus,sigeta2_inv

        sigeta2_inv=1./sigeta**2 

        if (nz .eq. 1) then                     ! 2+1 with Glauber (ah*N_coll+(1-ah)*N_wound)
          
          call wound(amassa,tt1,tt2,thplus,thminus)
          edens=ecenter*(ah*bincoll()+(1.-ah)*(tt1+tt2))/((ah*weightbin)+((1-ah)*weightwound))

        else   

          ! 3+1 used both by Hirano and Shenke. 
      !  The reference here is hep-ph 1004.1408v1 (Shenke, Jeon, Gale) Formulas 51,52
          call wound(amassa,tt1,tt2, thplus, thminus)
      H=H_of_eta(eta, deta, sigeta2_inv, thplus, thminus,center)
      HP=H_of_eta(eta, deta, sigeta2_inv, thplus, thminus,plus)
      HM=H_of_eta(eta, deta, sigeta2_inv, thplus, thminus,minus)

          if(tilting) then !initial energy density tilting applied
        edens=ecenter*(ah*bincoll()*H+2.*(1.-ah)*(tt1*HM+tt2*HP))/((ah*weightbin)&
                 &+((1-ah)*weightwound))
            !write(*,*) 'edens:',edens
          else
        edens=H*ecenter*(ah*bincoll()+(1.-ah)*(tt1+tt2))/((ah*weightbin)+((1-ah)*weightwound))
          endif
    endif

        end subroutine calc_enedens
    
! ***************************************************************************

    function H_of_eta(etac, etaflat, sigmaetasq_inv, thplus, thminus, simmetry_flag)
    !  The reference here is hep-ph 1004.1408v1 (Shenke, Jeon, Gale) Formulas 51,52
      implicit none
          real(8), intent(in) :: thplus, thminus
      real(8) H_of_eta
      real(8) etac, etaflat, sigmaetasq_inv, eta0
          integer simmetry_flag

          if(tilting) then !in the case of an initial energy density tilting,eta0=0
            eta0=0.
          else 
            eta0=0.5*log(((thminus+thplus)*gamma_col+(thminus-thplus)*gamma_col*beta_col)/&
                &((thminus+thplus)*gamma_col-(thminus-thplus)*gamma_col*beta_col)) 
          end if 

          if (abs(etac-eta0) .gt. (etaflat*0.5)) then 
             H_of_eta=exp(-1.0* (( abs(etac-eta0)-(0.5*etaflat) )**2) * 0.5* sigmaetasq_inv )
          else
             H_of_eta=1.0
          endif

          if(simmetry_flag .eq. plus) then
            if(etac-eta0 .lt. -etam) then
              H_of_eta=0.
            else if(etac-eta0 .le. etam) then
              H_of_eta=H_of_eta*(etac-eta0+etam)/(2.*etam)
            end if
          else if(simmetry_flag .eq. minus) then
            if(etac-eta0 .gt. etam) then
              H_of_eta=0.
            else if(etac-eta0 .ge. -etam) then
              H_of_eta=H_of_eta*(-etac+eta0+etam)/(2.*etam)
            end if
          end if

    end function H_of_eta

! ***************************************************************************

    subroutine calc_chargedens(rho)
    implicit none
    real(8) rho
    rho=rhocenter
    return
    end subroutine calc_chargedens
    
! ***************************************************************************

    subroutine calc_vel_long_nz_ueta(velo,xx,yy,eta)
    implicit none
    real(8) velo(3)
        real(8) :: vsound = 1/sqrt(3.)-1.e-6

        real(8) :: glg_ueta
        real(8) :: ueta
        real(8), intent(in) :: xx, yy, eta
        integer thickarray_index
        real(8) thick_plus, thick_minus, thick_zero, thick_limit, cut_factor, rplus, rminus

    
      velo(1)=0.
      velo(2)=0.
      velo(3)=0.

          ! here we define u^eta - yb is the beam rapidity, xx and eta the x and eta coordinates, t is the time
          ! ueta_coeff is chosen inside the parameter file (usually param.dat )
          ueta=(1./t)*sinh(yb-abs(eta))*tanh(ueta_coeff*xx)
          if(abs(eta) .ge. yb) ueta=0.
          glg_ueta=sqrt(1.+ueta*ueta*g_cov(3))
          velo(3)=ueta/glg_ueta


    return
    end subroutine calc_vel_long_nz_ueta

! ***************************************************************************

    function bincoll()

!   calculates number of binary collisions

    implicit none
        integer rad_index
        real(8) xshifted
    real(8) radius
    real(8) bincoll

    xshifted=xx+b1*0.5
    radius=sqrt(xshifted**2+yy**2)
    rad_index=nint(radius/dr+1)
    thplus=thick(rad_index)

    xshifted=xx-b1*0.5
    radius=sqrt(xshifted**2+yy**2)
    rad_index=nint(radius/dr+1)
    thminus=thick(rad_index)

    bincoll=sigma_in*thplus*thminus
    return

    end function bincoll
    
! ***************************************************************************

        subroutine wound(amassa,tA,tB, thplus, thminus)

!   calculates numbers of wounded nucleons

        implicit none
        integer rad_index
    real(8)  xshifted
        real(8)  radius
    real(8) amassa,tA,tB
        real(8),intent(out) :: thplus,thminus

        xshifted=xx+b1*0.5
        radius=sqrt(xshifted**2+yy**2)
        rad_index=nint(radius/dr+1)
        thplus=thick(rad_index)

        xshifted=xx-b1*0.5
        radius=sqrt(xshifted**2+yy**2)
        rad_index=nint(radius/dr+1)
        thminus=thick(rad_index)

        tA=sigma_in*thminus/amassa
        tA=(1.-tA)**amassa
        tA=1.-tA
        tA=tA*thplus

        tB=sigma_in*thplus/amassa
        tB=(1.-tB)**amassa
        tB=1.-tB
        tB=tB*thminus
 
        end subroutine wound

! ***************************************************************************

    subroutine thickness(radius,delta,acc)
!-------------------------------------
!   calculates thickness function of a given nucleus (radius)
!   for many values of r=(x,y) and stores into the vector <thick>
!       unless parameters are the same as those stored into the file thick_param.dat:
!       in this case, the subroutine read the values already computed from the file thickness.dat
!-------------------------------------
    implicit none
    real(8) acc
    real(8) rho1,rho2,dz,zmax  &
        ,r,delta_inv,espo,erre,z,r2,sum1,sum2,z1,z2
    real(8) ratio,check,radius,delta
        real(8) radiusold, deltaold, accold,dz1old,drold,bitbucket
        integer nmaxold

        real(8) largestx, largesty,rmax
        
    integer n,nstep,nrun,nc,filerror
        logical data_to_be_created
        integer allocate_result
007     format(10x,e14.8)
008     format(a10,e14.8)
009     format(a10,i14)
010     format(10x,i14)

  

        largestx=max(abs(x1),abs(x2))+bimpact
        largesty=max(abs(y1),abs(y2))
        rmax=sqrt(largestx*largestx+largesty*largesty)
        nmax=int(rmax/dr)+1        
        allocate(thick(1:nmax),STAT=allocate_result)
        if(allocate_result /=0) then
          write(*,*) "Error, I can't allocate thick array inside thickness subroutine, contained into init.f90"
          write(*,*)  "(source file init.f90)"
          call exit(1)
        end if

        if(pe0) then

          open(unit=39,status='OLD',file='thick_params.dat', iostat=filerror, form='formatted')
          if (filerror .eq. 0) then !this code will be executed only if the file already exists
           write(*,*) 'Checking if parameters for thickness computation are the same as in a previous run'
           read(39,007) radiusold
           read(39,007) deltaold
           read(39,007) accold
           read(39,007) dz1old
           read(39,007) drold
           read(39,010) nmaxold
           close(39)
           if((radiusold .eq. radius) .and. (deltaold .eq. delta) .and. (accold .eq. acc) .and. (dz1old .eq. dz1) .and.&
             & (drold .eq. dr) .and. (nmaxold .eq. nmax)) then
             data_to_be_created=.false.
             open(unit=2,status='OLD',file='thick.dat', iostat=filerror, form='formatted')
             if (filerror .ne. 0) then 
                data_to_be_created=.true.
                close(2)
             else
                write(*,*) 'Yes, they are the same, reading thickness values from old file thick.dat'
                do n=1,nmax
                   read(2,*,IOSTAT=filerror) bitbucket, thick(n)!we discard the first value
                   if(filerror .lt. 0) then
                     write(*,*) 'Ooops, file thick.dat endend before reading all values...'
                     write(*,*) "I'll compute thickness values again..."
                     data_to_be_created=.true.
                     close(2)
                   end if
                end do
                close(2)
             end if     
            else !some of the old parameters for thickness computation don't match the current ones
             data_to_be_created=.true.
             write(*,*) 'No, they are not the same:'
             write(*,*) 'Old radius:', radiusold,' - Current radius:', radius
             write(*,*) 'Old delta:', deltaold,' - Current delta:', delta
             write(*,*) 'Old acc:', accold,' - Current acc:', acc
             write(*,*) 'Old dz1old:', dz1old,' - Current radius:', dz1
             write(*,*) 'Old radius:', drold,' - Current radius:', dr
             write(*,*) 'Old nmax:', nmaxold,' - Current nmax:', nmax
             write(*,*) 'Computing thickenss function again with the new parameters.'
            end if
        else
           close(39)
           data_to_be_created=.true.
        end if
           
        if(data_to_be_created) then  

          open(unit=2,status='unknown',file='thick.dat')

          delta_inv=1./delta
          zmax=2.*radius

          check=0.
        
          do nc=1,nmax
             dz=2.*dz1

             r=(nc-1)*dr
             r2=r**2

             ratio=1000.
             nrun=0
             !  calculate thickness function at a given r=(x,y) 
             do while(ratio.gt.acc)
                dz=dz*0.5
                nstep=nint(2.*zmax/dz)
                nrun=nrun+1

                z1=-zmax
                z2=z1+dz
                sum1=0.
                sum2=0.

                do n=1,nstep


                   z=z1
                   erre=sqrt(r2+z**2)
                   espo=exp((erre-radius)*delta_inv)
                   rho1=roze/(espo+1.)

                   z=z2
                   erre=sqrt(r2+z**2)
                   espo=exp((erre-radius)*delta_inv)
                   rho2=roze/(espo+1.)

                   sum1=sum1+rho1
                   sum2=sum2+rho2

                   z1=z2
                   z2=z2+dz
               end do

               ratio=abs(2.*(sum2-sum1)/(sum1+sum2))

            end do

            if(pe0) write(2,*) r,sum1*dz,nrun

            check=check+r*sum1*dz

            thick(nc)=sum1*dz
          end do
        
          write(*,'(a,f10.5)')'check nucleus mass',check*2*3.1415927*dr
          close(2)
          open(unit=39,status='REPLACE',file='thick_params.dat', iostat=filerror, form='formatted')
             if (filerror .ne. 0) then
               write(*,*) "File thick_params.dat cannot be opened, so I'm forced to quit!"
               close(39)
               call exit(1)
             end if
             write(39,008) 'radius:   ',radius
             write(39,008) 'delta:    ',delta
             write(39,008) 'accuracy: ',acc
             write(39,008) 'dz1:      ',dz1
             write(39,008) 'dr:       ',dr
             write(39,009) 'nmax:     ',nmax
          close(39)
         end if !end if data_to_be_created
        end if !end if pe0

        if(prl) then
          call MPI_Bcast(thick, nmax, mpi_realtype, 0, icomm, ierr)
          call MPI_Barrier(MPI_COMM_WORLD,ierr)
        end if

end subroutine thickness

!***********************************************************************************************

subroutine dumpBgmc(ix,iy,iz,edens,Tpress)
implicit none
integer :: ix, iy, iz
real(8), intent(in) :: edens, Tpress
real(8) :: bx, by, bz, Bpress, Bpress_over_Tpress, reduce_factor

if(edens .lt. edens_for_B_dump) then
  Bpress=0.5*(v(ix,iy,iz,kbx)**2+v(ix,iy,iz,kby)**2+g_cov(3)*v(ix,iy,iz,kbz)**2)
  Bpress_over_Tpress=Bpress/Tpress 
  if(Bpress_over_Tpress .gt. B_th_press_ratio) then
     reduce_factor=sqrt(B_th_press_ratio/Bpress_over_Tpress)
     v(ix,iy,iz,kbx:kbz)=v(ix,iy,iz,kbx:kbz)*reduce_factor
  end if
end if

end subroutine dumpBgmc

!***********************************************************************************************

subroutine dumpB(ix,iy,iz,edens,Tpress)
implicit none
integer :: ix, iy, iz
real(8), intent(in) :: edens, Tpress
real(8) :: bx, by, bz, Bpress, Bpress_over_Tpress, reduce_factor

if(edens .lt. edens_for_B_dump) then
  Bpress=0.5*(vall(ix,iy,iz,kbx)**2+vall(ix,iy,iz,kby)**2+g_cov(3)*vall(ix,iy,iz,kbz)**2)
  Bpress_over_Tpress=Bpress/Tpress 
  if(Bpress_over_Tpress .gt. B_th_press_ratio) then
     reduce_factor=sqrt(B_th_press_ratio/Bpress_over_Tpress)
     vall(ix,iy,iz,kbx:kbz)=vall(ix,iy,iz,kbx:kbz)*reduce_factor
  end if
end if

end subroutine dumpB

!***********************************************************************************************

subroutine dumpB_subarray(edens,primitives_array)
implicit none
real(8), intent(in) :: edens
real(8), dimension(1:4) :: primitives_array
real(8) :: bx, by, bz, Bpress, Bpress_over_Tpress, reduce_factor, Tpress

Tpress=primitives_array(1)
bx=primitives_array(2)
by=primitives_array(3)
bz=primitives_array(4)

if(edens .lt. edens_for_B_dump) then
  Bpress=0.5*(bx**2+by**2+g_cov(3)*bz**2)
  Bpress_over_Tpress=Bpress/Tpress 
  if(Bpress_over_Tpress .gt. B_th_press_ratio) then
     reduce_factor=sqrt(B_th_press_ratio/Bpress_over_Tpress)
     primitives_array(2)=bx*reduce_factor
     primitives_array(3)=by*reduce_factor
     primitives_array(4)=bz*reduce_factor
  end if
end if

end subroutine dumpB_subarray


!***********************************************************************************************

subroutine glissando_init(xpos,ypos,coll,num_entries,num_events,rad,cutoff,amix)
implicit none
real(8),dimension(:),allocatable :: xpos,ypos
integer,dimension(:),allocatable :: coll
integer :: i,j,k,l,a,b,num_entries,num_events
integer :: cutoff,icut, jcut !cutoff limit in rad (smearing sigma) units and coordinated cutoff indexes
real(8), intent(in) :: rad, amix
real(8) :: norm_gauss, xdelta2, ydelta2, dens_gauss, sr2, tilting, feta, dx, dy

sr2 = 2.d0*rad*rad
norm_gauss=1.d0/(sr2*pi)

!DDD write(*,*) "1/norm_gauss: ", 1.d0/norm_gauss
!DDD write(*,*) "ybeam:  ", yb

dx=(x2-x1)/nx
dy=(y2-y1)/ny

icut=cutoff*int(rad/dx)
jcut=icut

do l=1,num_entries
  a = ceiling((xpos(l)-x1) / dx) 
  b = ceiling((ypos(l)-y1) / dy)

  !DDD write(*,*) "l: ", l, "  xpos: ", xpos(l), "  ypos: ", ypos(l), "  icut: ", icut, " coll: ", coll(l)

  do j = b - jcut, b + jcut
    if((j .gt. 0) .and. (j .lt. ny)) then
      ydelta2=-((y(j)-ypos(l))**2)/sr2
      do i = a - icut, a + icut
         if((i .gt. 0) .and. (i .lt. nx)) then
           xdelta2=-((x(i)-xpos(l))**2)/sr2
           dens_gauss = norm_gauss *  exp(xdelta2+ydelta2)
           !DDD write(*,*) "i: ", i, "x: ", x(i), " j: ", j, "y: ", y(j), " dens_gauss: ", dens_gauss
           do k=1,nz
              if(coll(l) .gt. 0) then
                 tilting=1.d0 + z(k)/yb
              else
                 tilting=1.d0 - z(k)/yb
              end if
              if(abs(z(k)) .le. deta/2.d0) then
                 feta = 1.0
              else if(abs(z(k)) .lt. yb) then
                   feta = exp(-(((abs(z(k))-deta/2.d0)/sigeta)**2)/2.d0)
              else
                   cycle
              end if
              !DDD write(*,*) "eta: ", z(k), " fEta: ", feta, " tilting: ", tilting, " density: ",&
              !DDD           & dens_gauss*feta*tilting*((1.d0 - amix) + abs(coll(l))*amix)
              vall(i,j,k,kpr)=vall(i,j,k,kpr)+dens_gauss*feta*tilting*((1.d0 - amix) + abs(coll(l))*amix)
           end do !k cycle
         end if
       end do ! i cycle 
     end if
   end do !j cycle
end do !l cycle

xpos=0.d0
ypos=0.d0
coll=0    

end subroutine glissando_init

!***********************************************************************************************

end module qgpinit

!***********************************************************************************************
subroutine init()
!-- Set grid and initial conditions

  use common
  use system_eqgp
!
  use qgpinit
  use viscous

  use out, only: restarted_run
!
  implicit none

  real(8),dimension(:),allocatable :: vv,uu
  integer :: ix,iy,iz,ixl,ixr,iyl,iyr,izl,izr
  real(8) :: rho0,rho1,r,f,rho,eps
  real(8) :: k,tk !this variable is tanh(k) for viscous Gubser flow initialization
  real(8),dimension(3) :: xini
  integer :: allocate_result
  real(8) :: dprdrh, dprden, edens, press, eta0,tr_weight,gamm,weight0
  real(8) :: eta1,xmax,ymax

  real(8) :: rad,esum,edif,beta1,arglog,rprime,tt1,tt2

  integer i !just a counter
 
  integer :: error_flag,filerror,iter_error
  real(8),allocatable, dimension(:,:,:) :: glaubermc_edens_array, glaubermc_edens_event_array
  integer,allocatable, dimension(:) :: displsin, recvcountsin
  
  !this is for initialization with tabulated initial entropy or energy density profile
  real(8),allocatable, dimension(:,:) :: raw_ed_array, ied_ed_array
  real(8),allocatable, dimension(:,:,:) :: ied_pr_array
  real(8) :: rh_bb, dprdrh_bb,dprden_bb
  

  real(8) gubdan_x, gubdan_y ! x and y positions for Gubser's initialization
  real(8) gubdan_en, gubdan_ux, gubdan_uy !energy density and velocities ux, uy for Gubser's initialization
  real(8) gubdan_glf !Lorentz factor for the Gubser's initialization
  real(8) gubdan_pixx, gubdan_piyy, gubdan_pixy, gubdan_pitt, gubdan_pitx, gubdan_pity, gubdan_pizz

  integer statfile
  real(8) xvalue, yvalue, edvalue, oldxvalue
  real(8), allocatable, dimension(:) :: x_array, y_array
  integer xdim, ydim
  logical find_ydim

  real(8) :: eta_alf, va_alf, B0_alf, rho_alf, pr_alf, gamma_alf, h_alf, energy_alf
  !variables for the ALFVEN 2D test
  real(8) :: dl, wl, diag, cst, snt, xrng, yrng, xrt, yrt, cf, vpar, vperp, Bperp
  !variable for the 2D cylindrical blast wave test
  real(8) :: rr

  real(8) :: sigma_gauge !variable for 1D MHD Bjorken flow

  !parameters for Lyutikov test
  real(8), parameter :: sigma_lyu=1.
  real(8), parameter :: B1_lyu=100.
  real(8), parameter :: scale_factor_lyu=1.e-4
  real(8) :: B_lyu
  !parameters for the Beraudo's test
  real(8), parameter :: B_berhad=16.
  real(8), parameter :: scale_factor_berhad=1.e-5

  !parameters for the Rotor test
  real(8), parameter :: rotor_radius=0.1
  real(8) :: rotor_r
  real(8), parameter :: rotor_omega=0.99

  !parameters for the Orszag-Tang test
  real(8), parameter :: vxmax_ot=0.6
  real(8), parameter :: vymax_ot=0.6
  real(8), parameter :: p_ot=9.
  real(8), parameter :: rho_ot=1.

  !variable for the local cons2prim in the EXT_CONS initialization
  real(8) :: s2,sx,sy,sz,et,d,glf,drhdpr,dendpr,pr_eos,df,dpr,en,pr,v2,rh,ww
  integer :: iter,end_index
  real, parameter :: tol=1.e-12
  real(8),dimension(1:5) :: cons_temp
  logical :: first_passage
  real(8),dimension(1:3) :: bitbucket_Efield

  !variables for the BLAST3D test
  real(8) :: r_internal, r_external 

  !variables for the ideal Gubser test
  real(8), parameter :: q_idGub=1.
  real(8), parameter :: T0_idGub=2.0
  real(8) :: k_idGub, temp_idGub

  !variables of the initialization with an external Glissando file
  integer, parameter :: array_dim_Glissando_init=20000000 !number of entries in the arrays that store Glissando i.c
  real(8), allocatable :: xpos_gliss(:), ypos_gliss(:) !positions of collisions in Glissando i.c.
  integer, allocatable :: coll_gliss(:) !number of collisions in Glissando i.c.
  integer :: num_Glissando_entries, num_Glissando_events, nwound_gliss, total_Glissando_events, nbin_gliss, l
  real(8), parameter :: srG=0.4 !smearing radius
  real(8), parameter :: alphaMix = 0.125 !wounded nucleons/binary mixing
  integer, parameter :: smear_cutoff=3 !cutoff in smearing radii
  integer, parameter :: scale_factor=14.0 !normalization factor to get the entropy density (or the energy density)
  real(8) :: bitbu_re1, bitbu_re2, bitbu_re3 !bitbuckets for reals
  integer :: bitbu_int1, bitbu_int2 !bitbucktes for integers
  logical, parameter :: initialize_with_entropy_density=.true. !true: initialization with entropy dens., false with energy dens.
  real(8) :: entropy_density_gliss, energy_density_gliss
  real(8), parameter :: dens_threshold=1.d-7 !density threshold for the conversion


007     format(10x,f10.5)

  error_flag=0
  iter_error=0


    call common_grid(x1,x2,y1,y2,z1,z2)

    !if we are restarting an old run we skip the variable initialization
    if(restarted_run) return

    allocate(vv(nv), uu(nv), STAT=allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Sorry, but I can't allocate memory for vv array inside sub init into init.f90... Quitting..."
      call exit(1)
    end if


    call print_param()

  vv(1:nv)=0.
  uu(1:nv)=0.

   t=tmin
   call system_metric()
   timeold=t
!


    if (init_type .eq. SODMHD) then

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)

           vv(kvx:kvz)=0.
           if(x(ix)<=0.5) then
             vv(krh)=1.
             vv(kpr)=1.
            if(mhd)   vv(kbx)=0.0
             if(mhd)  vv(kby)=0.0
           else
             vv(krh)=0.125
             vv(kpr)=0.1
         if(mhd)      vv(kbx)=0.0
        if(mhd)       vv(kby)=0.0
           end if

            
            

        if(.not. viscosity) then
           call system_prim2cons(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nou)=uu(1:nou) 
        end if
!
! for now we initialize only the non viscous variables
        v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do
!
 else if (init_type .eq. LYU) then
   write(*,*) "sigma is:",sigma_lyu
   write(*,*) "B1_lyu is:",B1_lyu
   
   B_lyu=B1_lyu*1.5/(sigma_lyu)**(2./3.)

   write(*,*) "B_lyu is:",B_lyu

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)

        vv(kvx:kvz)=0.
        vv(krh)=1.
        vv(kby)=sqrt(vv(krh)*sigma_lyu)
        vv(kpr)=(B_lyu*vv(kby)**2)/2.
        vv(kbx)=0.
        vv(kbz)=0.
        if(x(ix)>0.) then
          vv(krh)=vv(krh)*scale_factor_lyu
          vv(kpr)=vv(kpr)*scale_factor_lyu
          vv(kbx:kbz)=0.
        end if


        if(.not. viscosity) then
           call system_prim2cons(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nou)=uu(1:nou) 
        end if
!
! for now we initialize only the non viscous variables
        v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do
!

 else if (init_type .eq. ROTOR) then

   
   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)

        rotor_r=sqrt(xini(1)**2+xini(2)**2+xini(3)**2)
        vv(kvx:kvz)=0.
        vv(krh)=1.
        vv(kpr)=1.
        vv(kbx)=1.
        vv(kby)=0.
        vv(kbz)=0.
        if(rotor_r .le. rotor_radius) then
          vv(krh)=10.
          vv(kvx)=-rotor_omega*xini(2)/rotor_radius
          vv(kvy)=rotor_omega*xini(1)/rotor_radius
        end if


        if(.not. viscosity) then
           call system_prim2cons(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if

           call system_cons2prim(uu,vv,iter_error,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nou)=uu(1:nou) 
        end if
!
! for now we initialize only the non viscous variables
        v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do

 else if (init_type .eq. OT) then

  if(pe0) write(*,*) "Starting Orszang-Tang test"

   do iz=iz1,iz2
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)

        vv(kvx)=-vxmax_ot*sin(xini(2))
        vv(kvy)=vymax_ot*sin(xini(1))
        vv(kvz)=0.
        vv(krh)=rho_ot
        vv(kpr)=p_ot
        vv(kbx)=-sin(xini(2))
        vv(kby)=sin(2*xini(1))
        vv(kbz)=0.

        call system_prim2cons(vv,uu,error_flag)
        if(error_flag .gt. 0) then
          write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
          errcode=error_flag
          write(*,*) "Error code: ", errcode
          write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
          run_crashed=.true.
          return
         end if
         call system_cons2prim(uu,vv,iter_error,error_flag)
         if(error_flag .gt. 0) then
           write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
           errcode=error_flag
           write(*,*) "Error code: ", errcode
           write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
           run_crashed=.true.
           return
         end if
         u(ix,iy,iz,1:nou)=uu(1:nou) 
!
! for now we initialize only the non viscous variables
         v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do
!
 else if (init_type .eq. DBG) then

  if(pe0) write(*,*) "Starting simple debugging initialization"

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)

        vv(kvx:kvz)=0.
        vv(krh)=1.
        vv(kpr)=exp(-5.*(xini(1)**2+xini(2)**2))+0.001
        vv(kby)=0.1
        vv(kbx)=0.1
        vv(kbz)=0.

        call system_prim2cons(vv,uu,error_flag)
        if(error_flag .gt. 0) then
          write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
          errcode=error_flag
          write(*,*) "Error code: ", errcode
          write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
          run_crashed=.true.
          return
         end if
         call system_cons2prim(uu,vv,iter_error,error_flag)
         if(error_flag .gt. 0) then
           write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
           errcode=error_flag
           write(*,*) "Error code: ", errcode
           write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
           run_crashed=.true.
           return
         end if
         u(ix,iy,iz,1:nou)=uu(1:nou) 
!
! for now we initialize only the non viscous variables
         v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do
  if(pe0) write(*,*) "Flatbox initialization ended"
!
 else if (init_type .eq. BERHAD) then

   write(*,*) "B  is:",B_berhad
   
   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)

        vv(kvx:kvz)=0.
        vv(krh)=1.e-7
        vv(kpr)=1000.
        vv(kby)=sqrt(2.*vv(kpr)/B_berhad)
        vv(kbx)=0.
        vv(kbz)=0.
        if(x(ix)>0.) then
          vv(kpr)=vv(kpr)*scale_factor_berhad
          vv(kbx:kbz)=0.
        end if


        if(.not. viscosity) then
           call system_prim2cons(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nou)=uu(1:nou) 
        end if
!
! for now we initialize only the non viscous variables
        v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do
!
else if (init_type .eq. GLAUBER_GEOMETRIC) then
          !DDD
          !-------------------------------
          !     initialize the thickness vector <thick>
          !-------------------------------
          if(pe0) print*,'Calculating thickness function...'
          call thickness(radius,delta,acc)

          if(pe0) then
            print*,'Done!'
            print*,'**************'
          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       binary collisions and wounded nucleons at the origin and b=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------
            xx=0.
            yy=0.
            b1=0.
            weightbin=bincoll()
            weightbin_inv=1./weightbin
            call wound(amassa,tt1,tt2,thplus,thminus)
            weightwound=tt1+tt2
            weightwound_inv=1./weightwound
!       
            if(pe0) then
              open(unit=10,status='unknown',file='wound.dat')
              write(10,007) weightbin,weightwound
              close(10)
            end if
!
! ---   Now that the weight in the origin has been calculated,
! ---   b1 is set to the real(8) impact parameter
!
            b1=bimpact
!

! RHIC/Hirano initialization methods provided by A. De Pace
   if(mhd) then
    if(pe0) then
      write(*,*) "*** Loading initial magnetic field configuration ***"

      open(15,file=input_B_field,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) input_B_field, " cannot be opened, so I'm forced to quit!"
         close(15)
         call exit(1)
      end if

      do iz=1,nz
       do iy=1,ny
        do ix=1,nx
           read(15,*) vall(ix,iy,iz,kbx:kbz), bitbucket_Efield(1:3)
           vall(ix,iy,iz,kbx:kbz)= vall(ix,iy,iz,kbx:kbz)*B_amplify_factor
        end do
       end do
      end do
    
      close(15)
    end if !end of commands to be run by pe0 only

    if(prl) then
      allocate(recvcountsin(0:npe-1), displsin(0:npe-1), stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*) "Error, I can't allocate recvcountsin or displsin array for the MPI_Gatherv funcion (source file init.f90)"
        call exit(1)
      end if

      displsin(0)=0
      do i=0,npe-1
         if (i .lt. ipec) then
            recvcountsin(i)=nx/npe+1
         else
            recvcountsin(i)=nx/npe
         end if
         if(i .gt. 0) displsin(i)=displsin(i-1)+recvcountsin(i-1)
      end do

      do i=kbx,kbz
       do iz=iz1,iz2
        do iy=iy1,iy2
          call MPI_Scatterv(vall(1,iy,iz,i),recvcounts,displs,mpi_realtype,v(ix1,iy1+iy-1,iz1+iz-1,i),mx,&
                  &mpi_realtype,0,icomm,ierr)
        end do
       end do
      end do

      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    else
      v(:,:,:,kbx:kbz)=vall(:,:,:,kbx:kbz)
    end if !end if prl

   end if !end if mhd

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        xini(1:3)=(/x(ix),y(iy),z(iz)/)
        if(mhd) vv(kbx:kbz)=v(ix,iy,iz,kbx:kbz)
        call initio(nov,xini,vv,error_flag)

        if(mhd .and. dump_initial_B) then
          call eos_energy(vv(krh),eps,vv(kpr),errcode)
          call dumpB_subarray(eps,vv(kpr:kbz))
        end if

 
        if(error_flag .gt. 0) then
           write(*,*) "Sorry, but I'm not able to find the pressure corresponding to the minimum energy density enezero..."
           errcode=error_flag
           write(*,*) "Error code is: ", errcode
           run_crashed=.true.
           return
        end if

        if(.not. viscosity) then
           call system_prim2cons(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           !DDD here
           call system_cons2prim(uu,vv,iter_error,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if

           u(ix,iy,iz,1:nou)=uu(1:nou) 
        end if
!
! for now we initialize only the non viscous variables
       if(mhd) then
         v(ix,iy,iz,krh:kbz)=vv(krh:kbz)
       else
         v(ix,iy,iz,krh:kpr)=vv(krh:kpr)
       end if
    end do
   end do
  end do
!
  if(viscosity) then 
   do iz=iz1,iz2
    do iy=iy1,iy2
     do ix=ix1,ix2
           call viscous_initio(ix,iy,iz,nv,v,tmin,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling viscous_initio subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           vv(1:nv)=v(ix,iy,iz,1:nv)
           call system_prim2cons_0(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nv)=uu(1:nv)
     end do
    end do
   end do

  end if

else if (init_type .eq. EXT_GLISSANDO) then

    if(pe0) then 
      if(mhd) then
              write(*,*) "Sorry, but currently it not possible to initialize echo-qgp in MHD mode with Glissando"
              call exit(1)
      end if
      write(*,*) "*** Loading data from ",ext_glissando_file

      open(15,file=ext_glissando_file,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) ext_glissando_file, " cannot be opened, so I'm forced to quit!"
         close(15)
         call exit(1)
      end if


      allocate(xpos_gliss(1:array_dim_Glissando_init), ypos_gliss(1:array_dim_Glissando_init),&
      &coll_gliss(1:array_dim_Glissando_init),STAT=allocate_result)
      if(allocate_result /=0) then
         write(*,*) "Error, I can't allocate the arrays for the initialization with a Glissando external file"
         write(*,*)  "(source file init.f90)"
         call exit(1)
      end if
      total_glissando_events=0
      num_Glissando_events=0
      num_Glissando_entries=0

      do
        read(15,*,iostat=filerror) bitbu_re1, nwound_gliss, bitbu_int1, bitbu_re2, bitbu_int2, bitbu_re3
        if(IS_IOSTAT_END(filerror)) exit
        if(nwound_gliss .gt. 0) then
                if((num_Glissando_entries+nwound_gliss) .gt. array_dim_Glissando_init) then
                  call glissando_init&
                  &(xpos_gliss,ypos_gliss,coll_gliss,num_Glissando_entries,num_Glissando_events,srG,smear_cutoff,alphaMix)
                  num_Glissando_events=0
                  num_Glissando_entries=0
                end if
                num_Glissando_events=num_Glissando_events+1
                total_Glissando_events=total_Glissando_events+1
                l=num_Glissando_entries !inside the do loop we use l because it is shorter to write
                do i=1,nwound_gliss
                   l=l+1
                   read(15,*,iostat=filerror) xpos_gliss(l),ypos_gliss(l),coll_gliss(l)
                   if(IS_IOSTAT_END(filerror)) then
                           write(*,*) "Error, premature end of ", ext_glissando_file
                           close(15)
                           call exit(2)
                   end if
                end do
                num_Glissando_entries=l !we use again a more meaningful name
        end if
        read(15,*,iostat=filerror) nbin_gliss
        if(IS_IOSTAT_END(filerror)) exit
        if(nbin_gliss .gt. 0) then
           do i=1,nbin_gliss
              read(15,*,iostat=filerror) bitbu_re1, bitbu_re2, bitbu_int1
              if(IS_IOSTAT_END(filerror)) then
                 write(*,*) "Error, premature end of ", ext_glissando_file
                 close(15)
                 call exit(2)
              end if
           end do
        end if
      end do
               
      close(15) 
      write(*,*) "Events read: ",total_Glissando_events
      call glissando_init(xpos_gliss,ypos_gliss,coll_gliss,num_Glissando_entries,num_Glissando_events,srG,smear_cutoff,alphaMix)

      do iz=1,nz
       do iy=1,ny
        do ix=1,nx
            vall(ix,iy,iz,krh)=1.0!to be changed when, in the future, also baryon density will be provided by the external file 
            vall(ix,iy,iz,kvx:kvz)=0.
            if(initialize_with_entropy_density) then
               entropy_density_gliss=scale_factor*vall(ix,iy,iz,kpr)/total_Glissando_events 
               if(entropy_density_gliss .gt. dens_threshold) then
                 call eos_pressure_from_s(entropy_density_gliss,vall(ix,iy,iz,kpr),error_flag)
                 if(error_flag .ne. 0) then
                    write(*,*) "Error in converting the entropy density into pressure inside init (init.f90) at point: ",ix,iy&
                                &,iz," for entropy density value: ",entropy_density_gliss,". I quit."
                    call exit(2)
                 end if
               end if
            else !initialization with energy density
               energy_density_gliss=scale_factor*vall(ix,iy,iz,kpr)/total_Glissando_events 
               if(energy_density_gliss .gt. dens_threshold) then
                 call eos_pressure(vall(ix,iy,iz,krh),energy_density_gliss,vall(ix,iy,iz,kpr),bitbu_re1,bitbu_re2,error_flag)
                 if(error_flag .ne. 0) then
                    write(*,*) "Error in converting the energy density into pressure inside init (init.f90) at point: ",ix,iy&
                               &,iz," for entropy density value: ",energy_density_gliss,". I quit."
                    call exit(2)
                 end if
               end if
            end if
            if(vall(ix,iy,iz,kpr)>0) then 
                    write(*,*) ix, iy, iz, vall(ix,iy,iz,kpr)
            end if
            vall(ix,iy,iz,kpr)=vall(ix,iy,iz,kpr)+przero
        end do
       end do
      end do
      deallocate(xpos_gliss,ypos_gliss,coll_gliss)
    end if !end of commands to be run by pe0 only
    

    if(mhd) then !option alredy included future uses, albeit now not active because of a previous check that excludes mhd
      end_index=kbz
      if(rmhd) then
        end_index=kez
      end if
    else
      end_index=kpr
    end if

    if(prl) then
      allocate(recvcountsin(0:npe-1), displsin(0:npe-1), stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*) "Error, I can't allocate recvcountsin or displsin array for the MPI_Gatherv funcion (source file init.f90)"
        call exit(1)
      end if

      displsin(0)=0
      do i=0,npe-1
         if (i .lt. ipec) then
            recvcountsin(i)=nx/npe+1
         else
            recvcountsin(i)=nx/npe
         end if
         if(i .gt. 0) displsin(i)=displsin(i-1)+recvcountsin(i-1)
      end do

      do i=krh,end_index
       do iz=iz1,iz2
        do iy=iy1,iy2
          call MPI_Scatterv(vall(1,iy,iz,i),recvcounts,displs,mpi_realtype,v(ix1,iy1+iy-1,iz1+iz-1,i),mx,&
                  &mpi_realtype,0,icomm,ierr)
        end do
       end do
      end do

      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    else
      v(:,:,:,krh:end_index)=vall(:,:,:,krh:end_index)
    end if

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        vv(krh:end_index)=v(ix,iy,iz,krh:end_index)

        call system_prim2cons(vv,uu,error_flag)
        if(error_flag .gt. 0) then
           write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
           errcode=error_flag
           write(*,*) "Error code: ", errcode
           write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
           run_crashed=.true.
           return
        end if
        call system_cons2prim(uu,vv,iter_error,error_flag)
        if(error_flag .gt. 0) then
           write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
           errcode=error_flag
           write(*,*) "Error code: ", errcode
           write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
           run_crashed=.true.
           return
        end if

        u(ix,iy,iz,1:nou)=uu(1:nou) 
!
    end do
   end do
  end do
!
  if(viscosity) then 
   do iz=iz1,iz2
    do iy=iy1,iy2
     do ix=ix1,ix2
           call viscous_initio(ix,iy,iz,nv,v,tmin,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling viscous_initio subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           vv(1:nv)=v(ix,iy,iz,1:nv)
           call system_prim2cons_0(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nv)=uu(1:nv)
     end do
    end do
   end do

  end if


else if (init_type .eq. EXT_CONS) then
!at the moment, the EXT_CONS initialization works only in the ideal mhd case and this conditions are checkec in common.f90
    if(pe0) then 
      write(*,*) "*** Loading the initial values of conservative variables ***"

      open(15,file=ext_cons_file,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) ext_cons_file, " cannot be opened, so I'm forced to quit!"
         close(15)
         call exit(1)
      end if
 
      if(mhd) then
        write(*,*) "*** Loading also the initial magnetic field configuration ***"
        open(16,file=input_B_field,status='OLD',iostat=filerror)
        if (filerror .ne. 0) then
           write(*,*) input_B_field, " cannot be opened, so I'm forced to quit!"
           close(16)
           call exit(1)
        end if
      end if


      do iz=1,nz
       do iy=1,ny
        do ix=1,nx
            if(mhd) read(16,*) vall(ix,iy,iz,kbx:kbz), bitbucket_Efield(1:3)
            read(15,*) uu(kvx:kpr)
            uu(krh)=1.0 !to be changed when, in the future, also baryon density will be provided by the external file 
            if(uu(kpr) .lt. 1.e-10) then
              vall(ix,iy,iz,krh)=1.0!to be changed when, in the future, also baryon density will be provided by the external file 
              vall(ix,iy,iz,kvx)=0.
              vall(ix,iy,iz,kvy)=0.
              vall(ix,iy,iz,kvz)=0.
              vall(ix,iy,iz,kpr)=przero
            else
              call system_cons2prim0(uu,vv,iter_error,error_flag)
              vall(ix,iy,iz,krh)=1.0!to be changed when, in the future, also baryon density will be provided by the external file 
              if(vv(kpr) .gt. 0.1) then
                vall(ix,iy,iz,kvx:kvz)=vv(kvx:kvz)
              else
                vall(ix,iy,iz,kvx:kvz)=0.
              end if
              vall(ix,iy,iz,kpr)=vv(kpr)+przero
            end if

        end do
       end do
      end do
    
      close(15)
      if(mhd) close(16)
    end if !end of commands to be run by pe0 only

    if(mhd) then
      end_index=kbz
      if(rmhd) then
        end_index=kez
      end if
    else
      end_index=kpr
    end if


    if(prl) then
      allocate(recvcountsin(0:npe-1), displsin(0:npe-1), stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*) "Error, I can't allocate recvcountsin or displsin array for the MPI_Gatherv funcion (source file init.f90)"
        call exit(1)
      end if

      displsin(0)=0
      do i=0,npe-1
         if (i .lt. ipec) then
            recvcountsin(i)=nx/npe+1
         else
            recvcountsin(i)=nx/npe
         end if
         if(i .gt. 0) displsin(i)=displsin(i-1)+recvcountsin(i-1)
      end do

      do i=krh,end_index
       do iz=iz1,iz2
        do iy=iy1,iy2
          call MPI_Scatterv(vall(1,iy,iz,i),recvcounts,displs,mpi_realtype,v(ix1,iy1+iy-1,iz1+iz-1,i),mx,&
                  &mpi_realtype,0,icomm,ierr)
        end do
       end do
      end do

      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    else
      v(:,:,:,krh:end_index)=vall(:,:,:,krh:end_index)
    end if

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
     
        vv(krh:end_index)=v(ix,iy,iz,krh:end_index)

        if(mhd .and. dump_initial_B) then
          call eos_energy(vv(krh),eps,vv(kpr),errcode)
          call dumpB_subarray(eps,vv(kpr:kbz))
        end if

 
        call system_prim2cons(vv,uu,error_flag)
        if(error_flag .gt. 0) then
           write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
           errcode=error_flag
           write(*,*) "Error code: ", errcode
           write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
           run_crashed=.true.
           return
        end if
           !DDD here
        call system_cons2prim(uu,vv,iter_error,error_flag)
        if(error_flag .gt. 0) then
           write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
           errcode=error_flag
           write(*,*) "Error code: ", errcode
           write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
           run_crashed=.true.
           return
        end if

        u(ix,iy,iz,1:nou)=uu(1:nou) 
!
! for now we initialize only the non viscous variables
       if(mhd) then
         v(ix,iy,iz,krh:kbz)=vv(krh:kbz)
       else
         v(ix,iy,iz,krh:kpr)=vv(krh:kpr)
       end if
    end do
   end do
  end do
!
  if(viscosity) then 
   do iz=iz1,iz2
    do iy=iy1,iy2
     do ix=ix1,ix2
           call viscous_initio(ix,iy,iz,nv,v,tmin,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling viscous_initio subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           vv(1:nv)=v(ix,iy,iz,1:nv)
           call system_prim2cons_0(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nv)=uu(1:nv)
     end do
    end do
   end do

  end if

else if ( init_type .eq. SHOCK_TUBE_2D) then

  if(pe0) then
     write(*,*) "*** Performing init_type number ", SHOCK_TUBE_2D
     write(*,*) "*** 2D shock tube ***"
  end if
  
  do iz=iz1,iz2
    do iy=iy1,iy2
      do ix=ix1,ix2

     vv(krh)=1.
     vv(2:nv)=0. 
        
     if(y(iy) .lt. -x(ix)) then
       vv(kpr)=5.401411
     else 
        vv(kpr)=0.337588
     end if
    

      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)

    end do
   end do
  end do

  if(pe0) write(*,*) 'Initialization done...'

else if ( init_type .eq. BLAST3D) then

  r_internal = 0.4
  r_external = 0.8

  if(pe0) then
     write(*,*) "*** Performing init_type number ", BLAST3D
     write(*,*) "*** 3D blast wave ***"
  end if
  
  do iz=iz1,iz2
    do iy=iy1,iy2
      do ix=ix1,ix2

      rr=sqrt(x(ix)**2+y(iy)**2+z(iz)**2)
      vv(kvx:kvz)=0
      if(rr<r_internal) then
        vv(krh)=1.e-2
        vv(kpr)=1.
      else if(rr<r_external) then
        vv(krh)=1.e-4*(rr - r_internal)/(r_external - r_internal) + 1.e-2*(rr - r_external)/(r_internal - r_external);
        vv(kpr)=0.01*(rr - r_internal)/(r_external - r_internal) + 1.*(rr - r_external)/(r_internal - r_external);
      else
        vv(krh)=1.e-4
        vv(kpr)=0.01
      end if
      vv(kbx)=0.02
      vv(kby)=0.01
      vv(kbz)=-0.01

      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)

    end do
   end do
  end do

  if(pe0) write(*,*) 'Initialization done...'

else if ( init_type .eq. MHD1DBJ) then

  if(pe0) then
     write(*,*) "*** Performing init_type number ", MHD1DBJ
     write(*,*) "*** 1D MHD Bjorken flow ***"
  end if
  
  sigma_gauge=10. 
  do iz=iz1,iz2
    do iy=iy1,iy2
      do ix=ix1,ix2

      vv(kvx:kvz)=0
      vv(krh)=1.
      vv(kbx)=2.
      vv(kby)=3.
      vv(kbz)=0.
      vv(kpr)=10.
      vv(kbx:kby)=vv(kbx:kby)*sqrt(sigma_gauge*vv(kpr)*3/(vv(kbx)**2+vv(kby)**2))

      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)

    end do
   end do
  end do

  if(pe0) write(*,*) 'Initialization done...'
else  if (init_type .eq. ALFVEN2D) then

   eta_alf=1.0
   B0_alf=1.
   pr_alf=1.
   rho_alf=1.
   gamma_alf=2.
   vpar=0.
   h_alf=1.+gamma_alf/(gamma_alf-1)*pr_alf/rho_alf
   call eos_energy(rho_alf, energy_alf, pr_alf, errcode)
   if (abs(energy_alf/pr_alf - 3) .lt. 1.e-10) then ! this means we are using the eos p=e/3
        va_alf=sqrt(B0_alf**2/(4.*pr_alf+B0_alf**2*(1.+eta_alf**2))/(0.5*(1.+sqrt(1-((2.*eta_alf*B0_alf**2)/&
              &(4.*pr_alf+B0_alf**2*(1.+eta_alf**2)))**2))))
   else
   !Alfven speed for the ideal gas EOS
        va_alf=sqrt(B0_alf**2/(rho_alf*h_alf+B0_alf**2*(1.+eta_alf**2))/(0.5*(1.+sqrt(1-((2.*eta_alf*B0_alf**2)/&
              &(rho_alf*h_alf+B0_alf**2*(1.+eta_alf**2)))**2))))
   end if

   write(*,*) "Alfven speed va: ", va_alf
   !initialization as from:
   !https://svn.einsteintoolkit.org/cactus/EinsteinInitialData/GRHydro_InitData/trunk/src/GRHydro_AlfvenWaveM.F90
   xrng=x2-x1
   yrng=y2-y1
   diag=sqrt(xrng**2+yrng**2)
   dl=xrng*yrng/diag
   cst=yrng/diag
   snt=xrng/diag
   cf=2.*pi/dl

   do iz=iz1,iz2 
    do iy=iy1,iy2
     do ix=ix1,ix2
    
        vv(krh)=rho_alf
        vv(kpr)=pr_alf        
 
        xini(1:3)=(/x(ix),y(iy),z(iz)/)
        xrt=xini(1)*cst+xini(2)*snt      
        vperp=-va_alf*eta_alf*cos(cf*xrt)
        vv(kvx)=vpar*cst-vperp*snt
        vv(kvy)=vpar*snt+vperp*cst
        vv(kvz)=-va_alf*eta_alf*sin(cf*xrt)
        Bperp=eta_alf*B0_alf*cos(cf*xrt)
        vv(kbx)=B0_alf*cst-Bperp*snt
        vv(kby)=B0_alf*snt+Bperp*cst
        vv(kbz)=B0_alf*eta_alf*sin(cf*xrt)

        if(.not. viscosity) then
           call system_prim2cons(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nou)=uu(1:nou) 
        end if
!
! for now we initialize only the non viscous variables
        v(ix,iy,iz,1:nov)=vv(1:nov)
    end do
   end do
  end do

else if ( init_type .eq. SHOCK_TUBE_2D) then

  if(pe0) then
     write(*,*) "*** Performing init_type number ", SHOCK_TUBE_2D
     write(*,*) "*** 2D shock tube ***"
  end if
  
  do iz=iz1,iz2
    do iy=iy1,iy2
      do ix=ix1,ix2

     vv(krh)=1.
     vv(2:nv)=0. 
        
     if(y(iy) .lt. -x(ix)) then
       vv(kpr)=5.401411
     else 
        vv(kpr)=0.337588
     end if
    

      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)

    end do
   end do
  end do

  if(pe0) write(*,*) 'Initialization done...'

else if ( init_type .eq. SHEAR_VISCOUS_1D) then 

  if(pe0) then
     write(*,*) "*** Performing initialization number ", SHEAR_VISCOUS_1D
     write(*,*) "*** (1D viscous shear flow initialization) ***"
  end if
  
  do iz=iz1,iz2
  do iy=iy1,iy2
  do ix=ix1,ix2

     vv(1:nv)=0. 

       vv(kpr)=0.25
       vv(krh)=1.
       vv(kvy)=0.01*erf(x(ix)/(2.*sqrt(eta_over_s*t)))

      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_cons2prim subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)

    end do
   end do
  end do

  if(pe0)  write(*,*) 'Initialization done...'

 else if(init_type .eq. GLAUBER_MONTECARLO) then
  !Glauber-MC initial conditions
  if(pe0) then
     write(*,*) "*** Performing initialization number ",GLAUBER_MONTECARLO
     write(*,*) "*** Glauber-MonteCarlo initial conditions ***"
  end if
  v(:,:,:,kvx:)=0.!all variables, including, if present B field and viscosities, are initialized to 0.
  v(:,:,:,krh)=1.
  u(:,:,:,:)=0.

  allocate(glaubermc_edens_array(ix1:ix2,iy1:iy2,iz1:iz2), stat=allocate_result)
  if(allocate_result /=0) then
    write(*,*)  "Proc.0 - Error, I can't allocate glaubermc_edens_array"
    write(*,*)  "(source file init.f90)"
    call exit(1)
  end if
  glaubermc_edens_array=0.
  if(avg_events .eq. 1) then
    allocate(glaubermc_edens_event_array(ix1:ix2,iy1:iy2,iz1:iz2),stat=allocate_result)
    if(allocate_result /=0) then
      write(*,*)  "Proc.0 - Error, I can't allocate glaubermc_edens_event_array"
      write(*,*)  "(source file init.f90)"
      call exit(1)
    end if
    glaubermc_edens_event_array=0.
  end if

  !the creation of events it is done only by 1 processor, to avoid complications in writing the files with the collision data 
  if((pe0) .and. (id_of_the_run .eq. events_start)) then
     call generate_events
  end if
 !DDD end if !end section to be exectued only once

  do while(.true.)
     if(avg_events .eq. 0) then
       call generate_energy_density_profile(glaubermc_edens_array) 
       exit
     else !at the moment the only possibilities are 0 and 1, a third option is ruled out when printing the configuration
       call generate_energy_density_profile(glaubermc_edens_event_array) 
       glaubermc_edens_array=glaubermc_edens_array+glaubermc_edens_event_array
       glaubermc_edens_event_array=0.
       if(id_of_the_run .eq. events_stop) then
         if(mhd) id_of_the_run=events_start !we reset the counter, e.g. for the B field
         exit
      end if
      id_of_the_run=id_of_the_run+1
    end if
  end do


 if(avg_events .eq. 1) glaubermc_edens_array=glaubermc_edens_array/(events_stop-events_start+1)

 do iz=iz1,iz2
    do iy=iy1,iy2
       do ix=ix1,ix2
          call eos_pressure(rhocenter,glaubermc_edens_array(ix,iy,iz)+enezero,press,dprdrh,dprden,error_flag) 
          if(error_flag .gt. 0) then
            write(*,*) "An error occurred into the init subroutine when calling eos_pressure subroutine"
            errcode=error_flag
            write(*,*) "Error code: ", errcode
            write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
            run_crashed=.true.
            return
          end if
          v(ix,iy,iz,kpr)=press 
          !the minimum energy density has been already taken into account when calling eos_pressure
       end do
    end do
 end do

 if(mhd .and. (input_B_field .ne. "none")) then

  if(pe0) then
      write(*,*) "*** Loading the initial magnetic field configuration ***"
     
      open(15,file=input_B_field,status='OLD',iostat=filerror)
      if (filerror .ne. 0) then
         write(*,*) input_B_field, " cannot be opened, so I'm forced to quit!"
         close(15)
         call exit(1)
      end if
       do iz=1,nz
        do iy=1,ny
         do ix=1,nx
           read(15,*) vall(ix,iy,iz,kbx:kbz)
         end do
        end do
       end do

      close(15)
  end if !end if pe0

  if(prl) then
    allocate(recvcountsin(0:npe-1), displsin(0:npe-1), stat=allocate_result)
    if(allocate_result /=0) then
      write(*,*) "Error, I can't allocate recvcountsin or displsin array for the MPI_Gatherv funcion (source file init.f90)"
      call exit(1)
    end if

    displsin(0)=0
    do i=0,npe-1
      if (i .lt. ipec) then
         recvcountsin(i)=nx/npe+1
      else
         recvcountsin(i)=nx/npe
      end if
      if(i .gt. 0) displsin(i)=displsin(i-1)+recvcountsin(i-1)
    end do

    do i=kbx,kbz
      do iz=iz1,iz2
        do iy=iy1,iy2
           call MPI_Scatterv(vall(1,iy,iz,i),recvcounts,displs,mpi_realtype,v(ix1,iy1+iy-1,iz1+iz-1,i),mx,&
                &mpi_realtype,0,icomm,ierr)
        end do
       end do
    end do

    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  else
    v(:,:,:,kbx:kbz)=vall(:,:,:,kbx:kbz)
  end if !end if prl...

  if(dump_initial_B) then
   do iz=iz1,iz2
      do iy=iy1,iy2
         do ix=ix1,ix2
            call dumpBgmc(ix,iy,iz,glaubermc_edens_array(ix,iy,iz)+enezero,v(ix,iy,iz,kpr))
         end do
      end do
   end do
  end if!end if dump_initial_B
 end if !end (mhd .and. (input_B_field .ne. "none"))

 if(mhd .and. (input_B_field .eq. "none")) then
    !the input file for the B field is not equale to "none", so we create it
      if(kind_of_collision .eq. 3) then
        write(*,*) "Sorry, but currently in ECHO-QGP with Glauber Monte Carlo in. conditions"
        write(*,*) "only AA and dA collisions are possible if you create your initial B field"
        write(*,*) "at runtime. So, I stop here."
        call exit(1)
      end if

      if(pe0)  write(*,*) "*** Computing the initial magnetic field configuration ***"
      do while(.true.)
        if(avg_events .eq. 0) then
          call generate_B_field_profile
          exit
        else !at the moment the only possibilities are 0 and 1, a third option is ruled out when printing the configuration
          call generate_B_field_profile
          if(id_of_the_run .eq. events_stop) then
            exit
          end if
          id_of_the_run=id_of_the_run+1
        end if
      end do

      if(avg_events .eq. 1) v(:,:,:,kbx:kbz)=v(:,:,:,kbx:kbz)/(events_stop-events_start+1)

      if(dump_initial_B) then
       do iz=iz1,iz2
        do iy=iy1,iy2
         do ix=ix1,ix2
           if(dump_initial_B) call dumpBgmc(ix,iy,iz,glaubermc_edens_array(ix,iy,iz)+enezero,v(ix,iy,iz,kpr))
         end do
        end do
       end do
      end if !end if dump_initial_B
 end if !end if mhd .and....   

  do iz=iz1,iz2
  do iy=iy1,iy2
  do ix=ix1,ix2

    vv(1:nv)=v(ix,iy,iz,1:nv)
    if( .not. viscosity) then
     call system_prim2cons(vv,uu,error_flag)
     if(error_flag .gt. 0) then
       write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
       errcode=error_flag
       write(*,*) "Error code: ", errcode
       write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
       run_crashed=.true.
       return
     end if

     u(ix,iy,iz,1:nu)=uu(1:nu)
    end if

    end do
   end do
  end do

   if(viscosity) then
   do iz=iz1,iz2
    do iy=iy1,iy2
     do ix=ix1,ix2
           !DDD
           !write(*,*) 'init: primitive variables before viscous_initio:', v(ix,iy,iz,1:nv)
           !DDD
           call viscous_initio(ix,iy,iz,nv,v,tmin,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling viscous_initio subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           vv(1:nv)=v(ix,iy,iz,1:nv)
           !DDD
           !write(*,*) 'init: primitive variables after viscous_initio:', v(ix,iy,iz,1:nv)
           !write(*,*) 'init: vv variables after viscous_initio:', vv(1:nv)
           !DDD
           call system_prim2cons_0(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nv)=uu(1:nv)
           ! DDD
           !write(*,*) 'visc init - u:', u(ix,iy,iz,1:nv)
           ! DDD
      end do
     end do
    end do
   end if

  !AAADDD we didin't deallocate displsin and recvcountsin because they are used many times, but CAREFUL to future memory leaks

  if(avg_events .eq. 1) deallocate(glaubermc_edens_event_array)
  deallocate(glaubermc_edens_array)
 
  if(pe0) write(*,*) 'Initialization done...'

else if ( init_type .eq. GUBSER_IDEAL) then 
 if(pe0) then
     write(*,*) "*** Performing initialization number ", GUBSER_IDEAL
     write(*,*) "*** ideal Gubser's flow ***"
 end if

 if(prl) then
    if(pe0) then
       write(*,*) "Sorry, this initialization cannot be executed in parallel with mpi..."
       write(*,*) "Please, recompile the code serially"
    end if
    call exit(1)
  end if
  if(viscosity) then
     write(*,*) "Disabling viscosity"
     viscosity=.false.
     bulkvis=.false.
  end if
  if(eqstate .ne. 1) then
     write(*,*) "Setting equation of state to analytic"
     eqstate=1
     write(*,*) "Please, check also that you are using the proper relation for temperature and energy in temperature.def file"
  end if
  iz=1
  do ix=ix1,ix2
   do iy=iy1,iy2
      r=sqrt(x(ix)**2+y(iy)**2)
      temp_idGub=hbar*T0_idGub/tmin*(4*q_idGub**2*tmin**2)**(1.d0/3.d0)/&
     &(1.d0+2.d0*q_idGub**2*(tmin**2+r**2)+q_idGub**4*(tmin**2-r**2)**2)**(1.d0/3.d0)
      edens=47.5*pi**2/(3.d1*hbar**3)*temp_idGub**4
      vv(krh)=1.
      call eos_pressure(vv(krh),edens,vv(kpr),dprdrh,dprden,error_flag)
      k_idGub=atanh(2*(q_idGub**2)*tmin*r/(1.d0+(q_idGub*tmin)**2+(q_idGub*r)**2))
      vv(kvx)=x(ix)/r*tanh(k_idGub)
      vv(kvy)=y(iy)/r*tanh(k_idGub)
      vv(kvz)=0.d0

      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     !DDDwrite(*,*) v(ix,iy,iz,1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)
   end do
 end do
else if ( init_type .eq. GUBSER_VISCOUS) then 

  if(pe0) then 
     write(*,*) "*** Performing initialization number ", GUBSER_VISCOUS
     write(*,*) "*** viscous Gubser's flow with tabulated initial conditions ***"
  end if

  if(prl) then
     if(pe0) then
        write(*,*) "Sorry, this initialization cannot be executed in parallel with mpi..."
        write(*,*) "Please, recompile the code serially"
     end if
     call exit(1)
  end if

  write(*,*) "Disabling bulk viscosity"
  bulkvis=.false.
  write(*,*) "Setting equation of state to analytic"
  eqstate=1
  write(*,*) "Please, check also that you are using the proper relation for temperature and energy in temperature.def file"
  
  if(nz>1) then
     write(*,*) 'Sorry, this is only a 2D+1 initialization, I have to quit...'
     call exit(1)
  end if
  if(.not. viscosity) then
     write(*,*) 'Sorry, this is an initialization designed only for viscous case... I have to quit...'
     call exit(1)
  end if
  !check if the grid is correct
 if(nx .ne. 401) then
   write(*,*) 'Sorry, you need a grid of 401 cells along x and y'
   call exit(1)
 end if
 if(ny .ne. 401) then
   write(*,*) 'Sorry, you need a grid of 401 cells along x and y'
   call exit(1)
 end if
 if(x1 .ne. -10.025) then
   write(*,*) 'Sorry, x range must be from -10.025 to 10.025'
   call exit(1)
 end if
 if(y1 .ne. -10.025) then
   write(*,*) 'Sorry, y range must be from -10.025 to 10.025'
   call exit(1)
 end if

  
  open(unit=12,status='old',file='Initial_Profile_GubserFlow.dat', iostat=filerror)
    if (filerror .ne. 0) then 
       write(*,*) 'Sorry, but I am not able to locate the file:'
       write(*,*) 'Initial_Profile_GubserFlow.dat'
       write(*,*) 'containing tabulated initial conditions'
       write(*,*) 'Please, check that this file exists and re-run echo-qgp.'
    end if

  
  iz=1
  do ix=ix1,ix2
   do iy=iy1,iy2
     
     vv(krh)=1.
     read(12,*) gubdan_x,gubdan_y,gubdan_en,gubdan_ux,gubdan_uy,gubdan_pixx,gubdan_piyy,gubdan_pixy,gubdan_pitt,gubdan_pitx,&
               &gubdan_pity,gubdan_pizz
     gubdan_glf=1./sqrt(1.+gubdan_ux**2.+gubdan_uy**2.)
     vv(kvx)=gubdan_ux*gubdan_glf
     vv(kvy)=gubdan_uy*gubdan_glf
     call eos_pressure(vv(krh),gubdan_en*hbar,vv(kpr),dprdrh,dprden,error_flag)
     vv(kpixx)=gubdan_pixx*hbar
     vv(kpiyy)=gubdan_piyy*hbar
     vv(kpizz)=gubdan_pizz*hbar
     vv(kpixy)=gubdan_pixy*hbar
     vv(kpibu)=0.
   
      call system_prim2cons(vv,uu,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      call system_cons2prim(uu,vv,iter_error,error_flag)
      if(error_flag .gt. 0) then
         write(*,*) "An error occurred into the init subroutine when calling system_prim2cons subroutine"
         errcode=error_flag
         write(*,*) "Error code: ", errcode
         write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
         run_crashed=.true.
         return
      end if
      
     v(ix,iy,iz,1:nv)=vv(1:nv)
     u(ix,iy,iz,1:nu)=uu(1:nu)

    end do
   end do
  close(12)

  if(pe0) write(*,*) 'Initialization done...'

 else if(init_type .eq. TABULATED_INIT) then
  !initial entropy density profile
  if(pe0) then
     write(*,*) "*** Performing initialization number ",TABULATED_INIT
     write(*,*) "*** initial energy density or entropy density profile ***"
  end if

  if(nz>1) then
    write(*,*) "Sorry, this initialization is currently implemented for 2D simulations only..."
    call exit(1)
  end if
  v(:,:,:,kvx:)=0.
  v(:,:,:,krh)=1.
  u(:,:,:,:)=0.

    if(prl) then
      allocate(recvcountsin(0:npe-1), displsin(0:npe-1), stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*) "Error, I can't allocate recvcountsin or displsin array for the MPI_Gatherv funcion (source file init.f90)"
        call exit(1)
      end if
    end if

    if(pe0) then
      allocate(ied_pr_array(1:nx,1:ny,1), ied_ed_array(1:nx,1:ny), stat=allocate_result)
      if(allocate_result /=0) then
        write(*,*)  "Proc.0 - Error, I can't allocate ied_pr_array or ied_ed_array"
        write(*,*)  "(source file init.f90)"
        call exit(1)
      end if
      ied_pr_array=0.
      ied_ed_array=0.
    end if
   if (prl) then
    displsin(0)=0
    do i=0,npe-1
      if (i .lt. ipec) then
        recvcountsin(i)=nx/npe+1
      else
        recvcountsin(i)=nx/npe
      end if
      if(i .gt. 0) displsin(i)=displsin(i-1)+recvcountsin(i-1)
     end do
   end if
  !DDD end if !end section to be exectued only once

  if(prl) CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  if(pe0) then
    open(15,file=input_edf,status='OLD',iostat=filerror)
    if (filerror .ne. 0) then
      write(*,*) input_edf, " cannot be opened, so I'm forced to quit!"
      close(15)
      call exit(1)
    end if

    iz=1
    xdim=1
    ydim=1
    read(15,*,iostat=statfile) xvalue, yvalue, edvalue
    if(statfile<0) then
      write(*,*) input_edf, "has been correctly provided? It seems almost empty..."
      call exit(1)
    end if
    find_ydim=.true.
    do
      oldxvalue=xvalue
      read(15,*,iostat=statfile) xvalue, yvalue, edvalue
      if(statfile<0) exit
      if(xvalue==oldxvalue) then
        if(find_ydim) ydim=ydim+1
      else
        find_ydim=.false.
        xdim=xdim+1
      end if
    end do
    rewind(15) !we rewind the file, DDD - maybe, we could add an error check for this operation...
    
    allocate(raw_ed_array(1:xdim,1:ydim), x_array(1:xdim), y_array(1:ydim), stat=allocate_result)
    if(allocate_result/=0) then
      write(*,*) "An error occurred when allocating arrays for initialization with entropy density distribution"
      call exit(1)
    end if
    
    !DDD write(*,*) 'x dim:', xdim, 'y dim:', ydim
    do ix=1,xdim
      do iy=1,ydim
        read(15,*) x_array(ix), y_array(iy), raw_ed_array(ix,iy)
        !DDD write(*,*) x_array(ix), y_array(iy), raw_ed_array(ix,iy)
      end do
    end do

    close(15)

    call bilinear_interpolation(raw_ed_array,xdim,ydim,x_array,y_array, nx,ny,x,y,ied_ed_array) 

    do iy=1,ny
      do ix=1,nx
        if(ienentr .eq. 1) then
          call eos_pressure_from_s(ied_ed_array(ix,iy),press,error_flag) 
        else 
          call eos_pressure(rh_bb,ied_ed_array(ix,iy),press,dprdrh_bb,dprden_bb,error_flag)
        end if

        if(error_flag .gt. 0) then
          write(*,*) "An error occurred into the init subroutine when calling eos_pressure subroutine"
          errcode=error_flag
          write(*,*) "Error code: ", errcode
          write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
          run_crashed=.true.
          return
        end if
        ied_pr_array(ix,iy,1)=press+przero
      end do
    end do

  end if

 if(prl) then
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      do iz=iz1,iz2
        do iy=iy1,iy2
          call MPI_Scatterv(ied_pr_array(1,iy,iz),recvcountsin,displsin,mpi_realtype,v(ix1,iy1+iy-1,iz1+iz-1,kpr),mx,&
               &mpi_realtype,0,icomm,ierr)
        end do
      end do
      CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
      !AAA DDD think about deallocation
  else
    v(:,:,:,kpr)=max(ied_pr_array(:,:,:),przero)
  end if

  do iz=iz1,iz2
  do iy=iy1,iy2
  do ix=ix1,ix2

    vv(1:nv)=v(ix,iy,iz,1:nv)
    if( .not. viscosity) then
     call system_prim2cons(vv,uu,error_flag)
     if(error_flag .gt. 0) then
       write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
       errcode=error_flag
       write(*,*) "Error code: ", errcode
       write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
       run_crashed=.true.
       return
     end if

     u(ix,iy,iz,1:nu)=uu(1:nu)
    end if

    end do
   end do
  end do

   if(viscosity) then
   do iz=iz1,iz2
    do iy=iy1,iy2
     do ix=ix1,ix2
           !DDD
           !write(*,*) 'init: primitive variables before viscous_initio:', v(ix,iy,iz,1:nv)
           !DDD
           call viscous_initio(ix,iy,iz,nv,v,tmin,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling viscous_initio subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           vv(1:nv)=v(ix,iy,iz,1:nv)
           !DDD
           !write(*,*) 'init: primitive variables after viscous_initio:', v(ix,iy,iz,1:nv)
           !write(*,*) 'init: vv variables after viscous_initio:', vv(1:nv)
           !DDD
           call system_prim2cons_0(vv,uu,error_flag)
           if(error_flag .gt. 0) then
              write(*,*) "An error occurred into the init subroutine when calling system_prim2cons_0 subroutine"
              errcode=error_flag
              write(*,*) "Error code: ", errcode
              write(*,*) "Position: ix=",ix," iy=",iy," iz=",iz
              run_crashed=.true.
              return
           end if
           u(ix,iy,iz,1:nv)=uu(1:nv)
           ! DDD
           !write(*,*) 'visc init - u:', u(ix,iy,iz,1:nv)
           ! DDD
      end do
     end do
    end do
   end if

  !AAADDD we didin't deallocate displsin and recvcountsin because they are used many times, but CAREFUL to future memory leaks
 
  if(pe0) then
    write(*,*) 'Initialization done...'
  end if

end if
  if(pe0) then
    write(*,*) 'Writing config_summary.dat'
    call print_config_summary()
  end if

  !we set the array containing the variables values at the previuos time step
  vold(:,:,:,krh:kvold_end)=v(:,:,:,krh:kvold_end)
  uold(:,:,:,krh:kvold_end)=u(:,:,:,krh:kvold_end)

end subroutine init
