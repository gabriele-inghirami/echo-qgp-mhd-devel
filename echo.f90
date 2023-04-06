! *****************************************************************************
! *                                                                           *
! *  ECHO-QGP                                                                 * 
! *                                                                           *         
! *  Version: 1.5.0-alpha                                                     *
! *                                                                           *
! *  Copyright (C) 2015 The ECHO-QGP team                                     * 
! *                                                                           *
! *  File: echo.f90                                                           *
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

program echo
!-- ECHO, main file

  use parallel
  use common 
  use evolve
  use out
  use eos
  use system_eqgp, only: system_setup, system_metric
  use hypersurface

  implicit none

  integer :: iout,ilog,iter,irk,i1,i2,irate,hypout
  real(8) :: dt,tlog,time ! time is real execution time, not simulation time
  real(8) :: maximum_dt=-1. !maximum timestep passed as argument (called maxdt)
  integer :: error_flag


!-- Initializing processes

  call parallel_start()
  call common_get_all_parameters()

  
  if(init_type .eq. GLAUBER_MONTECARLO) then
    id_of_the_run=events_start
  else
    events_start=1
    events_stop=1
    id_of_the_run=events_stop
  end if

  if(maximum_dt .gt. 0) maxdt=maximum_dt
  
  dt=maxdt

  if(tmin .gt. tmax) then
    write(*,*) 'Sorry, but it is not possible to choose an ending time smaller than the initial time...'
    call exit(1)
  end if

  call system_setup(coordinates, system_metric, viscosity,nv,nu,nov,nou)
  call parallel_grid_setup(nx,ny,nz,mx,my,mz,ix0,ix1,ix2,iy0,iy1,iy2,iz0,iz1,iz2)

  call common_alloc

!-- Zeroth order bcs
!
  call set_eos_pointers(eqstate, eos_energy, eos_pressure, eos_sound)

  if (( eqstate .eq. 3) .and. (id_of_the_run .eq. events_start)) then
      call read_data_from_tabulated_file(entries, eos_tab_name, temperature_array, energy_density_array, pressure_array,&
      & sound_speed_array, min_en_dens, max_en_dens)
  end if

  call find_przero(przero, enezero,error_flag)
  if(prl) call MPI_Bcast(przero, 1, mpi_realtype, 0, icomm, ierr)
  if(error_flag .gt. 0) then
      write(*,*) "Sorry, but I'm not able to find the pressure corresponding to the minimum energy density enezero..."
      errcode=error_flag
      write(*,*) "Error code is: ", errcode
      run_crashed=.true.
      return
  end if

  call evolve_alloc
  call out_alloc

  if (pe0) then
    print*
    write(*,"(3x,20a)") 'ECHO-QGP STARTED'
    print*
  end if

!-- Loop over different initial conditions
  do while(id_of_the_run .le. events_stop)

    !-- Initial conditions

    iter=0; iout=1; hypout=0
    v=0.
    u=0.
    vnewstep=0.
    unewstep=0.
    uold=0.
    vold=0.
    deriv=0.
    run_crashed=.false.
    errcode=0

    call init()
 

    if(pe0) then
      write(*,'(a11,i4,a7,i4,a4,i4)') 'RUN NUMBER ',id_of_the_run,' - RUN ',id_of_the_run-events_start+1,' OF ',events_stop-&
           &events_start+1
      write(*,*) 'Initializations done...'
    end if
  


    if(.not. restarted_run) then

      t=tmin !already into init, here just for clearness
      timeold=t
      timeinterval=0
      call system_metric()
      call out_grid
      if(pe0) then
        write(*,*) 'Grid computed...'
        if(out_freeze) call get_hysufile_ready()
        write(*,*) 'Files for hypersurface computation written...'
      end if
      call out_vars_summary()
      if(pe0) write(*,*) 'Summary of variables written...'
      call out_vars(iout,0)
      if(pe0) then
         write(*,*) 'Copying the parameter file into the output directory...'
         call out_copy_tools
      end if
      !since we just starde probably it is not useful to dump data now
      !if (restart_type .ne. 0) call out_dump(iout)
      iout=iout+1
      if(out_freeze) hypout=hypout+1
    else
      if(pe0) then
        if(out_freeze) call get_hysufile_ready()
      end if
      call out_restart(iout)
      call system_metric()
      iout=iout+1
      restarted_run=.false. ! to avoid another restarting in case of multiple runs, e.g. with Glauber-MC initial conditions
    end if
   
    tout=t+dtout
    tlog=t+dtlog
    if(t .eq. tmin) then
       thyp=tmin
    else !alternative active when the program is restarted
       thyp=t+freeze_time_interval
    end if


  !-- Main cycle

    call system_clock(i1)

    do

      if (t>=tmax) exit
    
      if (temp_end .gt. 0) call check_temp()
      if (temp_check .eq. 1) exit
      
      if(mhd) then
        vnewstep=v(:,:,:,krh:kbz)
        unewstep=u(:,:,:,krh:kbz)
        if(rmhd) then
          vnewstep=v(:,:,:,krh:kez)
          unewstep=u(:,:,:,krh:kez)
        end if
      else
        vnewstep=v(:,:,:,krh:kpr)
        unewstep=u(:,:,:,krh:kpr)
      end if

      do irk=1,nrk
         call evolve_main(irk,dt)
         if(run_crashed) exit           
      end do
      
      if(run_crashed) exit           
          
      timeold=t
      t=t+dt
      g_cov_3_old=g_cov(3) !we keep track of the old value of g_cov(3) for computing time derivatives
      call system_metric() !we update the metric tensor soon after time updating
      iter=iter+1

      if (t>=tlog) then
        call system_clock(i2,irate)
        time=real(i2-i1,8)/real(irate,8)
        if(pe0) then
          open  (10,file='log')
          if((init_type .eq. GLAUBER_MONTECARLO) .and. (avg_events .eq. 1)) then
            write (10,'(a31,i4,a4,i4)') 'RUN AVERAGING I.C. FROM EVENT ',events_start,' TO ',events_stop
          else
            write (10,'(a11,i4,a4,i4)') 'RUN NUMBER ',id_of_the_run,' - RUN ',id_of_the_run-events_start+1,' OF ',events_stop&
                &-events_start+1
          end if
          write (10,'(3x,a,f15.8)') 'Simulation time:',t
          write (10,'(3x,a,f15.8)') 'Real time (sec):',time
          write (10,'(3x,a,i15  )') 'Iterations     :',iter
          if(time .gt. 0.) write (10,'(3x,a,f15.8)') 'Iterations/sec :',iter/time
          close (10)
       end if
       ilog=ilog+1
       tlog=tlog+dtlog
        if (tlog > tmax) tlog=tmax
      end if

      timeinterval=dt
      vold=vnewstep
      uold=unewstep

      if ((t>tout) .or. ((tout-t) .lt. 1.e-10))then
        tout=tout+dtout
        if((abs(tout-tmax) .lt. 1.e-10) .or. (tout .gt. tmax)) tout=tmax
        call out_vars(iout,1)
        !we dump the data for restarting now only if freezeout hypersurface is not computed
        if((restart_type .ne. 0) .and. (.not. out_freeze)) call out_dump(iout)
        iout=iout+1
      end if

      if(out_freeze) then
        if((t>thyp) .or. ((thyp-t) .lt. 1.e-10))then
          thyp=thyp+freeze_time_interval
          if((abs(thyp-tmax) .lt. 1.e-10) .or. (thyp .gt. tmax)) thyp=tmax
          call out_hyp(hypout)
          if(restart_type .ne. 0) call out_dump(iout-1)
          hypout=hypout+1
        end if
      end if

    end do

    if((.not. run_crashed) .and. (temp_check .eq. 1)) then
      if(pe0) then
        if(freeze_type .eq. 0) then
          print*,'Temperature is under the value established for stopping the simulation'
        else
          print*,'Energy density is under the value established for stopping the simulation'
        end if
        print*,'Printing variables again before leaving'
      end if
      call out_vars(iout,1)
      if(out_freeze) call out_hyp(hypout)
      if(restart_type .ne. 0) call out_dump(iout)
    end if

    call system_clock(i2,irate)
    time=real(i2-i1,8)/real(irate,8)

  !-- Closing processes

    if(run_crashed) then
      if (pe0) then
        print*
        print*,'RUN CRASHED BEFORE FINISHING...'
        print*
        print '(3x,a,f15.8)','Simulation time:',t
        print '(3x,a,f15.8)','Real time (sec):',time
        print '(3x,a,i15  )','Iterations     :',iter
        print '(3x,a,f15.8)','Iterations/sec :',iter/time
        print*
      end if
    else
      if (pe0) then
        print*
        print*,'RUN FINISHED'
        print*
        print '(3x,a,f15.8)','Simulation time:',t
        print '(3x,a,f15.8)','Real time (sec):',time
        print '(3x,a,i15  )','Iterations     :',iter
        print '(3x,a,f15.8)','Iterations/sec :',iter/time
        print*
      end if
    end if
    
    id_of_the_run=id_of_the_run+1 
  end do ! end of loop over runs
  call parallel_end

end program echo

! *****************************************************************************
